#!/usr/bin/env python3
"""
Author: Anna Korda

Filters AD4-redocked poses by a real, direct geometric criterion: is the
ligand's charged nitrogen atom within 5.0A of Asp116's CG carbon (the
carbon bonded to both carboxylate oxygens)? Same methodology as this
project's own analyze_asp_contact.py (protonated N via formal charge,
distance to Asp116's CG, 5.0A cutoff) -- pure RDKit, no LUNA project
machinery, since none of LUNA's receptor-wide interaction perception is
needed to answer a single distance question.

Molecules are parsed with sanitize=False and are never modified in any
way -- no kekulization, no valence/aromaticity perception, no coordinate
changes. Only the raw connection table (atoms, bonds, coordinates, formal
charges from M CHG lines) is read.

Same chunked/resumable/parallel/flat-tar pattern as the rest of this
project's large-file-count scripts.

Usage:
    python interaction_filter_fast.py \\
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
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

TARGET_RESNAME = "ASP"
TARGET_RESNUM = 116
TARGET_ATOM_NAME = "CG"
MAX_DIST_N_TO_CG = 5.0  # same cutoff as this project's analyze_asp_contact.py
COPY_BUF = 1024 * 1024 * 4


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def list_pose_files(poses_dir: Path):
    paths = []
    n_dirs = 0
    for root, dirs, files in os.walk(poses_dir):
        dirs.sort()
        n_dirs += 1
        if n_dirs % 200 == 0:
            log(f"  ...still listing, {n_dirs} directories walked, {len(paths):,} poses found so far")
        for name in sorted(files):
            if name.endswith(".sdf"):
                paths.append(Path(root) / name)
    return paths


def get_target_coord(receptor_pdb: Path):
    """Real Asp116 CG coordinate (the carboxylate carbon), any chain.
    Errors out if not found or ambiguous -- never silently guesses."""
    matches = []
    with open(receptor_pdb) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            resname = line[17:20].strip()
            atom_name = line[12:16].strip()
            try:
                resnum = int(line[22:26])
            except ValueError:
                continue
            if resname == TARGET_RESNAME and resnum == TARGET_RESNUM and atom_name == TARGET_ATOM_NAME:
                chain = line[21]
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                matches.append((chain, xyz))

    if not matches:
        raise SystemExit(f"Could not find {TARGET_RESNAME}{TARGET_RESNUM}/{TARGET_ATOM_NAME} in {receptor_pdb}")
    if len(matches) > 1:
        chains = ", ".join(f"chain {c}: {xyz}" for c, xyz in matches)
        raise SystemExit(f"Found {len(matches)} copies of {TARGET_RESNAME}{TARGET_RESNUM}/{TARGET_ATOM_NAME} "
                         f"in {receptor_pdb} ({chains}) -- ambiguous, refusing to guess which one.")
    chain, xyz = matches[0]
    log(f"Found exactly one {TARGET_RESNAME}{TARGET_RESNUM} {TARGET_ATOM_NAME} -- "
        f"chain {chain}, coordinates {xyz}")
    return xyz


def _find_charged_n_distances(sdf_path, target_coord):
    """Parse without sanitization -- molecule is read exactly as written,
    never restructured. Formal charges (M CHG) and coordinates are part of
    the raw connection table, not sanitization output, so this is safe."""
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)
    mol = supplier[0] if len(supplier) else None
    if mol is None:
        return None, "rdkit failed to parse"

    conf = mol.GetConformer()
    distances = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N" or atom.GetFormalCharge() <= 0:
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        d = ((pos.x - target_coord[0]) ** 2 + (pos.y - target_coord[1]) ** 2
             + (pos.z - target_coord[2]) ** 2) ** 0.5
        distances.append(d)

    if not distances:
        return None, "no charged N atom found"
    return distances, ""


def _check_one_compound(sdf_path, target_coord):
    distances, error = _find_charged_n_distances(sdf_path, target_coord)
    if distances is None:
        return False, error
    return min(distances) <= MAX_DIST_N_TO_CG, ""


def _process_chunk(chunk_paths, chunk_idx, target_coord, parts_dir_str):
    parts_dir = Path(parts_dir_str)
    part_path = parts_dir / f"part_{chunk_idx:05d}.tar.gz"
    report_path = parts_dir / f"report_{chunk_idx:05d}.csv"
    if part_path.exists() and report_path.exists():
        return chunk_idx, None, None

    stage_dir = parts_dir / f".stage_{chunk_idx}"
    stage_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for sdf_path in chunk_paths:
        compound_id = sdf_path.stem
        try:
            passed, error = _check_one_compound(sdf_path, target_coord)
        except Exception as exc:
            rows.append({"compound_id": compound_id, "passed": "error", "error": str(exc)})
            continue
        rows.append({"compound_id": compound_id, "passed": passed, "error": error})
        if passed:
            # copy2 preserves the file byte-for-byte -- the molecule that
            # gets kept is the original SDF, never anything RDKit touched.
            shutil.copy2(sdf_path, stage_dir / sdf_path.name)

    with open(report_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["compound_id", "passed", "error"])
        w.writeheader()
        w.writerows(rows)

    list_file = parts_dir / f"list_{chunk_idx:05d}.txt"
    passing = [p.name for p in stage_dir.glob("*.sdf")]
    tmp_part = parts_dir / f"part_{chunk_idx:05d}.tar.gz.partial-{os.getpid()}"
    if passing:
        list_file.write_text("\n".join(passing) + "\n")
        r = subprocess.run(["tar", "-czf", str(tmp_part), "-C", str(stage_dir), "-T", str(list_file)],
                            capture_output=True, text=True)
        list_file.unlink(missing_ok=True)
    else:
        r = subprocess.run(["tar", "-czf", str(tmp_part), "--files-from", os.devnull],
                            capture_output=True, text=True)
    shutil.rmtree(stage_dir, ignore_errors=True)
    if r.returncode != 0:
        tmp_part.unlink(missing_ok=True)
        return chunk_idx, r.stderr, None
    tmp_part.rename(part_path)

    n_pass = sum(1 for row in rows if row["passed"] is True)
    n_fail = sum(1 for row in rows if row["passed"] is False)
    n_err = sum(1 for row in rows if row["passed"] == "error")
    return chunk_idx, None, (n_pass, n_fail, n_err)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conformation", required=True)
    ap.add_argument("--poses_dir", required=True)
    ap.add_argument("--receptor_pdb", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--chunk_size", type=int, default=1000,
                     help="per-compound cost is trivial (pure distance check), so this is "
                          "sized for frequent progress updates, not to amortize per-chunk cost")
    args = ap.parse_args()

    poses_dir = Path(args.poses_dir)
    receptor_pdb = Path(args.receptor_pdb)
    if not receptor_pdb.is_file():
        raise SystemExit(f"{receptor_pdb} not found")
    if not poses_dir.is_dir():
        raise SystemExit(f"{poses_dir} not found")

    target_coord = get_target_coord(receptor_pdb)

    out_root = Path(args.out_dir) / f"interactions_passed_{args.conformation}"
    parts_dir = out_root / ".parts"
    parts_dir.mkdir(parents=True, exist_ok=True)

    workers = args.workers or len(os.sched_getaffinity(0))
    log(f"Listing pose SDFs under {poses_dir} ...")
    t_list0 = time.time()
    paths = list_pose_files(poses_dir)
    log(f"Listing done in {time.time() - t_list0:.0f}s.")
    chunks = [paths[i:i + args.chunk_size] for i in range(0, len(paths), args.chunk_size)]
    n_chunks = len(chunks)
    n_total = len(paths)

    have = {int(p.name.split("_")[1].split(".")[0])
            for p in parts_dir.glob("part_*.tar.gz")
            if (parts_dir / f"report_{p.name.split('_')[1].split('.')[0]}.csv").exists()}
    todo = [i for i in range(n_chunks) if i not in have]
    log(f"{n_total:,} poses, {n_chunks} chunks of up to {args.chunk_size:,}, "
        f"{len(have)} already done, {workers} workers")

    total_pass = total_fail = total_err = 0
    t0 = time.time()
    if todo:
        done_count = 0
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_process_chunk, chunks[i], i, target_coord, str(parts_dir)): i
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
                log(f"  chunk {done_count}/{len(todo)} done ({time.time() - t0:.0f}s elapsed, "
                    f"ETA {eta:.0f}s) -- processed so far: pass={total_pass} fail={total_fail} "
                    f"err={total_err} (of {total_pass + total_fail + total_err:,}/{n_total:,} compounds)")

    log("all chunks done -- merging per-chunk reports and concatenating tars ...")
    report_path = out_root / f"interaction_filter_report_{args.conformation}.csv"
    grand_pass = grand_fail = grand_err = 0
    with open(report_path, "w", newline="") as out_f:
        w = csv.writer(out_f)
        w.writerow(["compound_id", "passed", "error"])
        for i in range(n_chunks):
            with open(parts_dir / f"report_{i:05d}.csv") as in_f:
                for row in csv.DictReader(in_f):
                    w.writerow([row["compound_id"], row["passed"], row["error"]])
                    if row["passed"] == "True":
                        grand_pass += 1
                    elif row["passed"] == "False":
                        grand_fail += 1
                    else:
                        grand_err += 1

    tar_path = out_root / f"interactions_passed_{args.conformation}.tar.gz"
    tmp_out = f"{tar_path}.partial"
    with open(tmp_out, "wb") as out_f:
        for i in range(n_chunks):
            with open(parts_dir / f"part_{i:05d}.tar.gz", "rb") as in_f:
                shutil.copyfileobj(in_f, out_f, COPY_BUF)
    os.rename(tmp_out, tar_path)
    shutil.rmtree(parts_dir, ignore_errors=True)

    elapsed = time.time() - t0
    size_gb = os.path.getsize(tar_path) / (1024 ** 3)

    # Human-readable summary log, inside the output folder itself, separate
    # from the per-compound CSV report.
    summary_path = out_root / f"run_summary_{args.conformation}.txt"
    with open(summary_path, "w") as f:
        f.write(f"Interaction filter run summary -- conformation {args.conformation}\n")
        f.write(f"Finished: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Criterion: charged N within {MAX_DIST_N_TO_CG} A of "
                f"{TARGET_RESNAME}{TARGET_RESNUM} CG {target_coord}\n")
        f.write(f"Receptor: {receptor_pdb}\n")
        f.write(f"Poses dir: {poses_dir}\n\n")
        f.write(f"Starting number of compounds: {n_total:,}\n")
        f.write(f"Passed:  {grand_pass:,}\n")
        f.write(f"Failed:  {grand_fail:,}\n")
        f.write(f"Errored: {grand_err:,}\n")
        f.write(f"\nOutput tar: {tar_path.name} ({size_gb:.2f} GB)\n")
        f.write(f"Full per-compound report: {report_path.name}\n")

    log(f"Wrote {tar_path} ({size_gb:.2f} GB)")
    log(f"Wrote {report_path}")
    log(f"Wrote {summary_path}")
    log(f"TOTAL: pass={grand_pass:,} fail={grand_fail:,} error={grand_err:,} "
        f"out of {n_total:,} compounds ({elapsed:.0f}s)")
    log(f"extract poses with:  tar -ixzf {tar_path.name}   (-i required: concatenated multi-part gzip)")


if __name__ == "__main__":
    main()
