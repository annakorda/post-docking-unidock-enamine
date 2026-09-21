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
        --poses_dir vs_results/ad4_redock_c5 \\
        --receptor_pdb receptors/c5.pdb \\
        --out_dir vs_results \\
        --workers 64
"""
import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Same pattern this project's own run_benchmark.py/extract_top_poses.py use
# to read AD4 scores -- unidock_tools SDF output carries them as ENERGY=
# tags under a <Uni-Dock RESULT> property (benchmark.md). Applied per-pose
# (on that pose's own text block, not the whole file) so a compound's score
# always belongs to whichever specific pose it's paired with.
ENERGY_RE = re.compile(r"ENERGY=\s*(-?[\d.]+)")

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


def _split_sdf_blocks(text):
    """Split raw SDF text into per-pose blocks on the $$$$ delimiter, each
    block keeping its own trailing $$$$ line exactly as it was. Pure text
    slicing -- no chemistry library involved, so no risk of a parser
    silently normalizing anything."""
    blocks = []
    start = 0
    idx = text.find("$$$$")
    while idx != -1:
        end = idx + 4
        if text[end:end + 2] == "\r\n":
            end += 2
        elif text[end:end + 1] == "\n":
            end += 1
        blocks.append(text[start:end])
        start = end
        idx = text.find("$$$$", start)
    if text[start:].strip():
        blocks.append(text[start:])
    return blocks


def _parse_poses(sdf_path, target_coord):
    """Parse without sanitization -- molecules are read exactly as written,
    never restructured. Formal charges (M CHG) and coordinates are part of
    the raw connection table, not sanitization output, so this is safe.

    unidock's real defaults (confirmed via `unidock --help`, not assumed)
    are --num_modes 9 and --energy_range 3 (kcal/mol) -- AD4 redock SDFs
    can contain up to 9 poses per compound, and every pose actually written
    is already guaranteed by unidock itself to be within 3 kcal/mol of the
    best-scored pose. So checking only the first (best-ranked) pose misses
    real, near-optimal alternative poses that could make the contact even
    when the top-ranked-by-score pose doesn't -- this checks every pose.

    Distance and score are both read from the SAME per-pose text block
    (RDKit for geometry, this project's own ENERGY_RE for score) so they're
    always tied to the exact same pose, never mismatched.

    Returns (blocks, per_pose) where blocks is the raw text split (reused
    by the reordering step, avoiding a second read+split of the file) and
    per_pose is a list of (distance_or_None, score_or_None) per block."""
    from rdkit import Chem

    text = sdf_path.read_text()
    blocks = _split_sdf_blocks(text)

    per_pose = []
    for block in blocks:
        mol = Chem.MolFromMolBlock(block, sanitize=False, removeHs=False)
        distance = None
        if mol is not None and mol.GetNumConformers() > 0:
            conf = mol.GetConformer()
            for atom in mol.GetAtoms():
                if atom.GetSymbol() != "N" or atom.GetFormalCharge() <= 0:
                    continue
                pos = conf.GetAtomPosition(atom.GetIdx())
                d = ((pos.x - target_coord[0]) ** 2 + (pos.y - target_coord[1]) ** 2
                     + (pos.z - target_coord[2]) ** 2) ** 0.5
                if distance is None or d < distance:
                    distance = d
        m = ENERGY_RE.search(block)
        score = float(m.group(1)) if m else None
        per_pose.append((distance, score))

    return blocks, per_pose


def _check_one_compound(sdf_path, target_coord):
    """Returns (passed, error, winning_idx, winning_dist, winning_score, blocks).
    Among poses within MAX_DIST_N_TO_CG of the target, the winner is the
    one with the best (most negative) AD4 score -- salt-bridge criterion
    first (which poses even qualify), then best score among those that do.
    blocks is returned too so the reordering step can reuse the same split
    instead of re-reading the file."""
    blocks, per_pose = _parse_poses(sdf_path, target_coord)
    if not blocks:
        return False, "rdkit failed to parse", None, None, None, blocks

    passing = [(i, d, s) for i, (d, s) in enumerate(per_pose)
               if d is not None and d <= MAX_DIST_N_TO_CG]
    if not passing:
        if all(d is None for d, _ in per_pose):
            return False, "no charged N atom found in any pose", None, None, None, blocks
        return False, "", None, None, None, blocks

    scored = [(i, d, s) for i, d, s in passing if s is not None]
    if scored:
        winning_idx, winning_dist, winning_score = min(scored, key=lambda t: t[2])
    else:
        # no ENERGY= tag found on any qualifying pose -- fall back to the
        # closest one rather than fail the compound outright.
        winning_idx, winning_dist, _ = min(passing, key=lambda t: t[1])
        winning_score = None

    return True, "", winning_idx, winning_dist, winning_score, blocks


def _write_with_winning_pose_first(blocks, winning_idx, orig_text, out_path):
    """Move the winning pose's block to the front of the file, leaving
    every other block's bytes untouched and in their original relative
    order. Falls back to a plain write of the original text (no
    reordering) if the blocks don't reproduce the original file exactly --
    never write something we're not certain preserves the original data."""
    if "".join(blocks) != orig_text or not (0 <= winning_idx < len(blocks)):
        out_path.write_text(orig_text)
        return False

    reordered = [blocks[winning_idx]] + blocks[:winning_idx] + blocks[winning_idx + 1:]
    out_path.write_text("".join(reordered))
    return True


def _process_chunk(chunk_paths, chunk_idx, target_coord, parts_dir_str):
    parts_dir = Path(parts_dir_str)
    part_path = parts_dir / f"part_{chunk_idx:05d}.tar.gz"
    report_path = parts_dir / f"report_{chunk_idx:05d}.csv"
    reorder_report_path = parts_dir / f"reordered_{chunk_idx:05d}.csv"
    scores_path = parts_dir / f"scores_{chunk_idx:05d}.csv"
    if (part_path.exists() and report_path.exists() and reorder_report_path.exists()
            and scores_path.exists()):
        return chunk_idx, None, None

    stage_dir = parts_dir / f".stage_{chunk_idx}"
    stage_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    reordered_rows = []
    score_rows = []
    for sdf_path in chunk_paths:
        compound_id = sdf_path.stem
        try:
            passed, error, winning_idx, winning_dist, winning_score, blocks = \
                _check_one_compound(sdf_path, target_coord)
        except Exception as exc:
            rows.append({"compound_id": compound_id, "passed": "error", "error": str(exc)})
            continue
        rows.append({"compound_id": compound_id, "passed": passed, "error": error})
        if passed:
            out_path = stage_dir / sdf_path.name
            if winning_idx == 0:
                # already the first pose in the file -- plain copy, no
                # reordering needed.
                shutil.copy2(sdf_path, out_path)
            else:
                orig_text = sdf_path.read_text()
                did_reorder = _write_with_winning_pose_first(blocks, winning_idx, orig_text, out_path)
                if did_reorder:
                    reordered_rows.append({"compound_id": compound_id,
                                            "winning_pose_index_0based": winning_idx,
                                            "distance_A": f"{winning_dist:.3f}"})

            # Score of the winning (salt-bridge-satisfying) pose specifically
            # -- salt-bridge criterion decides which poses even qualify,
            # then the best score among those qualifying poses wins.
            score_rows.append({"compound_id": compound_id,
                                "winning_pose_score": "" if winning_score is None else f"{winning_score:.4f}"})

    with open(report_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["compound_id", "passed", "error"])
        w.writeheader()
        w.writerows(rows)

    with open(reorder_report_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["compound_id", "winning_pose_index_0based", "distance_A"])
        w.writeheader()
        w.writerows(reordered_rows)

    with open(scores_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["compound_id", "winning_pose_score"])
        w.writeheader()
        w.writerows(score_rows)

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
    return chunk_idx, None, (n_pass, n_fail, n_err, len(reordered_rows))


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

    def _chunk_idx_done(i):
        cid = f"{i:05d}"
        return ((parts_dir / f"part_{cid}.tar.gz").exists()
                and (parts_dir / f"report_{cid}.csv").exists()
                and (parts_dir / f"reordered_{cid}.csv").exists())

    have = {i for i in range(n_chunks) if _chunk_idx_done(i)}
    todo = [i for i in range(n_chunks) if i not in have]
    log(f"{n_total:,} poses, {n_chunks} chunks of up to {args.chunk_size:,}, "
        f"{len(have)} already done, {workers} workers")

    total_pass = total_fail = total_err = total_reordered = 0
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
                    n_pass, n_fail, n_err, n_reordered = counts
                    total_pass += n_pass
                    total_fail += n_fail
                    total_err += n_err
                    total_reordered += n_reordered
                rate = done_count / (time.time() - t0) if time.time() > t0 else 0
                eta = (len(todo) - done_count) / rate if rate else float("inf")
                log(f"  chunk {done_count}/{len(todo)} done ({time.time() - t0:.0f}s elapsed, "
                    f"ETA {eta:.0f}s) -- processed so far: pass={total_pass} fail={total_fail} "
                    f"err={total_err} reordered={total_reordered} "
                    f"(of {total_pass + total_fail + total_err:,}/{n_total:,} compounds)")

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

    reorder_path = out_root / f"reordered_compounds_{args.conformation}.csv"
    grand_reordered = 0
    with open(reorder_path, "w", newline="") as out_f:
        w = csv.writer(out_f)
        w.writerow(["compound_id", "winning_pose_index_0based", "distance_A"])
        for i in range(n_chunks):
            reordered_chunk_path = parts_dir / f"reordered_{i:05d}.csv"
            if not reordered_chunk_path.exists():
                continue
            with open(reordered_chunk_path) as in_f:
                for row in csv.DictReader(in_f):
                    w.writerow([row["compound_id"], row["winning_pose_index_0based"], row["distance_A"]])
                    grand_reordered += 1

    # Ranked hit list -- every passing compound's winning-pose score, sorted
    # best (most negative) first, so it's immediately obvious which
    # compounds are the strongest candidates.
    scored_compounds = []
    n_no_score = 0
    for i in range(n_chunks):
        scores_chunk_path = parts_dir / f"scores_{i:05d}.csv"
        if not scores_chunk_path.exists():
            continue
        with open(scores_chunk_path) as in_f:
            for row in csv.DictReader(in_f):
                if row["winning_pose_score"] == "":
                    n_no_score += 1
                    continue
                scored_compounds.append((row["compound_id"], float(row["winning_pose_score"])))
    scored_compounds.sort(key=lambda pair: pair[1])  # most negative (best) first

    ranked_path = out_root / f"ranked_hits_{args.conformation}.csv"
    with open(ranked_path, "w", newline="") as out_f:
        w = csv.writer(out_f)
        w.writerow(["rank", "compound_id", "winning_pose_score"])
        for rank, (compound_id, score) in enumerate(scored_compounds, start=1):
            w.writerow([rank, compound_id, f"{score:.4f}"])
    if n_no_score:
        log(f"WARNING: {n_no_score} passing compound(s) had no ENERGY= tag on their "
            f"winning pose and are excluded from {ranked_path.name} (still present in "
            f"{report_path.name})")

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
        f.write(f"Reordered (winning pose wasn't already first in its SDF): {grand_reordered:,}\n")
        f.write(f"Passing compounds missing an ENERGY= tag (excluded from ranking): {n_no_score:,}\n")
        f.write(f"\nOutput tar: {tar_path.name} ({size_gb:.2f} GB)\n")
        f.write(f"Full per-compound report: {report_path.name}\n")
        f.write(f"Reordered compounds detail: {reorder_path.name}\n")
        f.write(f"Ranked hit list (best score first): {ranked_path.name}\n")

    log(f"Wrote {tar_path} ({size_gb:.2f} GB)")
    log(f"Wrote {report_path}")
    log(f"Wrote {reorder_path} ({grand_reordered:,} compounds reordered)")
    log(f"Wrote {ranked_path} ({len(scored_compounds):,} ranked compounds)")
    log(f"Wrote {summary_path}")
    log(f"TOTAL: pass={grand_pass:,} fail={grand_fail:,} error={grand_err:,} "
        f"reordered={grand_reordered:,} out of {n_total:,} compounds ({elapsed:.0f}s)")
    log(f"extract poses with:  tar -ixzf {tar_path.name}   (-i required: concatenated multi-part gzip)")


if __name__ == "__main__":
    main()
