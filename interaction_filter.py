#!/usr/bin/env python3
"""
Author: Anna Korda

Filters AD4-redocked poses by a real receptor contact: does the ligand make
an Ionic or Salt bridge interaction with residue 116 (Asp116, the salt-bridge
residue this project's DUDE benchmark validated AD4 on)? Uses LUNA
(https://luna-toolkit.readthedocs.io) for interaction perception -- LUNA
perceives charged/ionizable atoms itself via RDKit's standard pharmacophore
feature definitions, so no manual SDF-charge parsing is needed here.

Same chunked/resumable/parallel pattern used throughout this project: list
all pose SDFs once, split into chunks, each worker runs one LUNA
LocalProject per chunk (nproc=1 inside -- parallelism is at the chunk level,
not inside LUNA), passing compounds get copied flat into that chunk's
staging dir and packed into part_NNNNN.tar.gz. Parts concatenate at the end
(gzip streams concatenate without recompression; extract with `tar -ixzf`).

Usage:
    python interaction_filter.py \\
        --conformation c5 \\
        --poses_dir /users/gpcr/annak/ultra-large/vs_results/ad4_redock_c5 \\
        --receptor_pdb /users/gpcr/annak/ultra-large/receptors/c5.pdb \\
        --out_dir /users/gpcr/annak/ultra-large/vs_results \\
        --workers 64
"""
import argparse
import csv
import os
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Module-level (not deferred into the worker like the rest of this file's
# LUNA usage) because LigandOnlyProject below must be picklable by
# EntryResults.save() -- a class defined inside a function has no resolvable
# module path for pickle to find it by, which is a real error hit locally
# ("Can't pickle local object ..."), not a theoretical concern.
from luna.projects import LocalProject, EntryResults
from luna.mol.entry import MolFileEntry
from luna.interaction.filter import InteractionFilter
from luna.interaction.calc import InteractionCalculator

TARGET_RESNAME = "ASP"
TARGET_RESNUM = 116
PASS_TYPES = {"Ionic", "Salt bridge"}
COPY_BUF = 1024 * 1024 * 4


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def list_pose_files(poses_dir: Path):
    """Recursively find pose SDFs -- AD4 redock output is bucket-nested,
    not flat (same layout cn_finalize_ad4.py/mr_finalize_ad4.py packed)."""
    paths = []
    for root, dirs, files in os.walk(poses_dir):
        dirs.sort()
        for name in sorted(files):
            if name.endswith(".sdf"):
                paths.append(Path(root) / name)
    return paths


def _touches_target_residue(atm_grp):
    return any(c.resname == TARGET_RESNAME and c.id[1] == TARGET_RESNUM
               for c in atm_grp.compounds)


def _check_one_compound(entry, entry_result):
    """Return True if this compound makes an Ionic or Salt bridge contact
    with the target residue, on either side of the interaction."""
    hits = entry_result.interactions_mngr.filter_by_types(PASS_TYPES)
    for it in hits:
        if _touches_target_residue(it.src_grp) or _touches_target_residue(it.trgt_grp):
            return True
    return False


def _make_inter_calc():
    calc = InteractionCalculator(inter_filter=InteractionFilter.new_pli_filter())
    # LUNA's default computes the FULL interaction profile for every nearby
    # atom-group pair -- hydrophobic, weak/halogen/chalcogen bonds, all three
    # stacking types, repulsive, multipolar, ion-dipole, metal coordination,
    # and a fully generic Atom-Atom vdW/proximal/clash check across every
    # atom pair regardless of chemistry -- even though PASS_TYPES only ever
    # needs Ionic and Salt bridge. Salt bridge is defined by LUNA as Ionic +
    # a coincident Hydrogen bond between the same groups (confirmed in
    # calc.py's dependent-interaction code), so Donor/Acceptor stays;
    # everything else here is computed and then thrown away by our own
    # filter regardless, so skip computing it in the first place.
    calc.funcs = {
        ("NegativelyIonizable", "PositivelyIonizable"): [calc.calc_ionic],
        ("NegIonizable", "PosIonizable"): [calc.calc_ionic],
        ("Negative", "Positive"): [calc.calc_ionic],
        ("Donor", "Acceptor"): [calc.calc_hbond],
    }
    return calc


class LigandOnlyProject(LocalProject):
    """LUNA's own LocalProject._process_entry (projects.py) calls
    calc_interactions(atm_grps_mngr.atm_grps) with no nb_atm_grps -- which
    means BOTH sides of the pairwise interaction search default to the
    entire cached neighborhood, receptor residues included, not just the
    ligand. Confirmed by profiling a real local run: is_valid_pair alone
    was ~30% of total runtime, called ~11,000+ times for a single ligand,
    because LUNA was needlessly generating and filtering receptor-residue
    to receptor-residue candidate pairs we never use -- we only ever check
    the ligand's own contacts. This override (copied from LUNA's own
    source, changed only where marked) passes just the ligand's own atom
    groups as the search target, keeping the full neighborhood only as the
    (already cached) thing being searched against. There's no `opts` hook
    for this in the public LocalProject API, hence the subclass instead of
    a config change."""

    def _process_entry(self, entry):
        start = time.time()
        self._log("debug", "Starting entry processing: %s." % entry.to_string())
        try:
            if isinstance(entry, MolFileEntry) is False:
                self._validate_entry_format(entry)

            pkl_file = "%s/chunks/%s.pkl.gz" % (self.working_path, entry.to_string())
            if self.append_mode and os.path.exists(pkl_file):
                return

            pdb_parser, structure, ligand = self._parse_complex(entry)
            add_h = self._decide_hydrogen_addition(pdb_parser.get_header(), entry)

            atm_grps_mngr = self._perceive_chemical_groups(entry, structure[0], ligand,
                                                            add_h, self.cache)
            atm_grps_mngr.entry = entry

            # THE FIX (everything else in this method is copied verbatim
            # from LUNA's own _process_entry): only the ligand's own
            # groups as the search target.
            ligand_grps = [g for g in atm_grps_mngr.atm_grps if ligand in g.compounds]

            calc_func = self.inter_calc.calc_interactions
            interactions_mngr = calc_func(ligand_grps, nb_atm_grps=atm_grps_mngr.atm_grps)
            interactions_mngr.entry = entry

            atm_grps_mngr.merge_hydrophobic_atoms(interactions_mngr)

            if self.binding_mode_filter is not None:
                interactions_mngr.filter_out_by_binding_mode(self.binding_mode_filter)

            ifp = None
            if self.calc_ifp:
                ifp = self._create_ifp(atm_grps_mngr)
            mfp = None
            if self.calc_mfp:
                mfp = self._create_mfp(entry)

            entry_results = EntryResults(entry, atm_grps_mngr, interactions_mngr, ifp, mfp)
            entry_results.save(pkl_file)

            csv_file = "%s/results/interactions/%s.csv" % (self.working_path, entry.to_string())
            interactions_mngr.to_csv(csv_file)

            if self.out_pse:
                from luna.interaction.view import InteractionViewer
                pse_path = self.pse_path or "%s/results/pse/" % self.working_path
                pse_file = "%s/%s.pse" % (pse_path, entry.to_string())
                piv = InteractionViewer(add_directional_arrows=True)
                piv.new_session([(entry, interactions_mngr, entry.pdb_file)], pse_file)

            self._log("debug", "Processing of entry '%s' finished successfully."
                      % entry.to_string())
        except Exception:
            self._log("debug", "Processing of entry '%s' failed." % entry.to_string())
            raise

        proc_time = time.time() - start
        self._log("debug", "Processing of entry '%s' took %.2fs." % (entry.to_string(), proc_time))


def _process_chunk(chunk_paths, chunk_idx, receptor_dir, receptor_id,
                    parts_dir_str, work_dir_str):
    parts_dir = Path(parts_dir_str)
    part_path = parts_dir / f"part_{chunk_idx:05d}.tar.gz"
    report_path = parts_dir / f"report_{chunk_idx:05d}.csv"
    if part_path.exists() and report_path.exists():
        return chunk_idx, None, None

    work_dir = Path(work_dir_str) / f"chunk_{chunk_idx}"
    stage_dir = work_dir / "passing_sdfs"
    stage_dir.mkdir(parents=True, exist_ok=True)

    entries = []
    entry_to_path = {}
    for sdf_path in chunk_paths:
        compound_id = sdf_path.stem
        entry = MolFileEntry.from_mol_file(receptor_id, compound_id, str(sdf_path),
                                            mol_obj_type="rdkit", is_multimol_file=False)
        entries.append(entry)
        entry_to_path[compound_id] = sdf_path

    opts = {
        "entries": entries,
        "working_path": str(work_dir / "luna_proj"),
        "pdb_path": receptor_dir,
        "overwrite_path": True,
        # add_h=False looked safer on paper (trust AD4's own protonation) but is
        # actually broken for external MOL-file ligands with explicit H: LUNA's
        # AtomGroupPerceiver sets keep_hydrog = not add_h, and when keep_hydrog
        # is True, its internal atom count includes H while mol_obj's heavy-atom
        # count doesn't -- guaranteed MoleculeSizeError crash. Confirmed by an
        # actual local run and by reading groups.py directly (not guessed).
        # add_h=True is the only path every LUNA tutorial actually exercises.
        "add_h": True,
        "ph": 7.4,
        "amend_mol": False,      # only affects receptor PDB het-residue perception, not our poses --
                                  # and this receptor has residues with missing sidechain atoms
                                  # (see prepare_receptor.sh's --allow_bad_res comment) that make
                                  # amend_mol's atom-count cross-check crash too, confirmed locally
        "calc_ifp": False,
        "out_pse": False,
        "use_cache": True,       # build the receptor's own chemical-group perception ONCE
                                  # from the first compound in this chunk, reuse it for the
                                  # other 99 -- valid here because every compound in a chunk
                                  # is docked into the same fixed box/pocket, so one cache
                                  # covers them all. This is the real fix for the ~23s/compound
                                  # throughput seen in the first production run: without it,
                                  # LUNA redoes full receptor-neighborhood perception (expand_selection
                                  # is hardcoded True internally, not something opts can disable)
                                  # from scratch for every single compound.
        "nproc": None,           # serial inside -- parallelism is at the chunk level
        "inter_calc": _make_inter_calc(),
        # logging_enabled=True writes a real log file to GPFS per operation --
        # on a filesystem this whole project already proved has real per-file
        # I/O latency problems (the MareNostrum redocking slowdown, earlier
        # this session). Off, matching a known-fast reference LUNA script.
        "logging_enabled": False,
        "verbosity": 0,
    }

    rows = []
    try:
        proj = LigandOnlyProject(**opts)
        proj.run()

        for entry in entries:
            compound_id = entry.mol_id
            sdf_path = entry_to_path[compound_id]
            try:
                res = proj.get_entry_results(entry)
                passed = _check_one_compound(entry, res)
            except Exception as exc:
                rows.append({"compound_id": compound_id, "passed": "error", "error": str(exc)})
                continue

            rows.append({"compound_id": compound_id, "passed": passed, "error": ""})
            if passed:
                shutil.copy2(sdf_path, stage_dir / sdf_path.name)
    except Exception as exc:
        return chunk_idx, f"chunk {chunk_idx} LUNA project failed: {exc}", None

    with open(report_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["compound_id", "passed", "error"])
        w.writeheader()
        w.writerows(rows)

    pass_list = stage_dir.glob("*.sdf")
    tmp_part = parts_dir / f"part_{chunk_idx:05d}.tar.gz.partial-{os.getpid()}"
    import subprocess
    list_file = parts_dir / f"list_{chunk_idx:05d}.txt"
    list_file.write_text("\n".join(p.name for p in stage_dir.glob("*.sdf")) + "\n")
    r = subprocess.run(["tar", "-czf", str(tmp_part), "-C", str(stage_dir), "-T", str(list_file)],
                        capture_output=True, text=True)
    list_file.unlink(missing_ok=True)
    if r.returncode != 0:
        tmp_part.unlink(missing_ok=True)
        return chunk_idx, r.stderr, None
    tmp_part.rename(part_path)

    shutil.rmtree(work_dir, ignore_errors=True)

    n_pass = sum(1 for row in rows if row["passed"] is True)
    n_fail = sum(1 for row in rows if row["passed"] is False)
    n_err = sum(1 for row in rows if row["passed"] == "error")
    return chunk_idx, None, (n_pass, n_fail, n_err)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conformation", required=True)
    ap.add_argument("--poses_dir", required=True, help="extracted ad4_redock_<conf>/ directory")
    ap.add_argument("--receptor_pdb", required=True, help="path to <receptor_id>.pdb")
    ap.add_argument("--out_dir", required=True, help="vs_results/ root -- output goes to <out_dir>/interactions_filtered_<conf>/")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--chunk_size", type=int, default=100,
                     help="compounds per LUNA project run -- keep small, this is compute-heavy "
                          "unlike the pure file-copy chunking elsewhere in this project")
    args = ap.parse_args()

    poses_dir = Path(args.poses_dir)
    receptor_pdb = Path(args.receptor_pdb)
    receptor_dir = str(receptor_pdb.parent)
    receptor_id = receptor_pdb.stem
    if not receptor_pdb.is_file():
        raise SystemExit(f"{receptor_pdb} not found")
    if not poses_dir.is_dir():
        raise SystemExit(f"{poses_dir} not found")

    out_root = Path(args.out_dir) / f"interactions_filtered_{args.conformation}"
    parts_dir = out_root / ".parts"
    work_dir = out_root / ".work"
    parts_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    workers = args.workers or len(os.sched_getaffinity(0))
    log(f"Listing pose SDFs under {poses_dir} ...")
    paths = list_pose_files(poses_dir)
    chunks = [paths[i:i + args.chunk_size] for i in range(0, len(paths), args.chunk_size)]
    n_chunks = len(chunks)

    have = {int(p.name.split("_")[1].split(".")[0])
            for p in parts_dir.glob("part_*.tar.gz")
            if (parts_dir / f"report_{p.name.split('_')[1].split('.')[0]}.csv").exists()}
    todo = [i for i in range(n_chunks) if i not in have]
    log(f"{len(paths):,} poses, {n_chunks} chunks of up to {args.chunk_size}, "
        f"{len(have)} already done, {workers} workers")

    total_pass = total_fail = total_err = 0
    if todo:
        t0 = time.time()
        done_count = 0
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_process_chunk, chunks[i], i, receptor_dir, receptor_id,
                                    str(parts_dir), str(work_dir)): i
                       for i in todo}
            for fut in as_completed(futures):
                idx, err, counts = fut.result()
                if err:
                    print(err, file=sys.stderr)
                    raise SystemExit(f"chunk {idx} failed")
                done_count += 1
                if counts:
                    n_pass, n_fail, n_err = counts
                    total_pass += n_pass
                    total_fail += n_fail
                    total_err += n_err
                rate = done_count / (time.time() - t0) if time.time() > t0 else 0
                eta = (len(todo) - done_count) / rate if rate else float("inf")
                log(f"  chunk {done_count}/{len(todo)} this run done "
                    f"({time.time() - t0:.0f}s elapsed, ETA {eta:.0f}s) "
                    f"-- running totals: pass={total_pass} fail={total_fail} err={total_err}")

    log("all chunks done -- merging per-chunk reports and concatenating tars ...")
    report_path = out_root / f"interaction_filter_report_{args.conformation}.csv"
    with open(report_path, "w", newline="") as out_f:
        w = csv.writer(out_f)
        w.writerow(["compound_id", "passed", "error"])
        grand_pass = grand_fail = grand_err = 0
        for i in range(n_chunks):
            with open(parts_dir / f"report_{i:05d}.csv") as in_f:
                r = csv.DictReader(in_f)
                for row in r:
                    w.writerow([row["compound_id"], row["passed"], row["error"]])
                    if row["passed"] == "True":
                        grand_pass += 1
                    elif row["passed"] == "False":
                        grand_fail += 1
                    else:
                        grand_err += 1

    tar_path = out_root / f"interactions_filtered_{args.conformation}.tar.gz"
    tmp_out = f"{tar_path}.partial"
    with open(tmp_out, "wb") as out_f:
        for i in range(n_chunks):
            part = parts_dir / f"part_{i:05d}.tar.gz"
            with open(part, "rb") as in_f:
                shutil.copyfileobj(in_f, out_f, COPY_BUF)
    os.rename(tmp_out, tar_path)
    shutil.rmtree(parts_dir, ignore_errors=True)
    shutil.rmtree(work_dir, ignore_errors=True)

    size_gb = os.path.getsize(tar_path) / (1024 ** 3)
    log(f"Wrote {tar_path} ({size_gb:.2f} GB)")
    log(f"Wrote {report_path}")
    log(f"TOTAL: pass={grand_pass} fail={grand_fail} error={grand_err} "
        f"out of {grand_pass + grand_fail + grand_err} compounds")
    log(f"extract poses with:  tar -ixzf {tar_path.name}   (-i required: concatenated multi-part gzip)")


if __name__ == "__main__":
    main()
