#!/usr/bin/env python3
"""
Author: Anna Korda

One-off helper: pulls the top N compounds from ranked_hits_<conf>.csv,
locates each one's real original SDF via scores.csv's own path column
(fast, direct lookup -- no scanning the huge interactions_passed tar or
walking the whole ad4_redock_<conf>/ tree), reproduces the exact same
winning-pose-first reordering interaction_filter_fast.py already applied
(reordered_compounds_<conf>.csv tells us which pose won), and bundles the
result + the receptor PDB into one small tar for visual review.

Molecules are never touched beyond pure text pose-block reordering, same
guarantee as interaction_filter_fast.py itself.

Usage:
    python extract_review_set.py --conformation c5 \\
        --vs_results_dir ~/ultra-large/vs_results \\
        --receptor_pdb ~/ultra-large/receptors/c5.pdb \\
        --n 10 --out_dir ~/ultra-large/test_review_c5
"""
import argparse
import csv
import shutil
import tarfile
from pathlib import Path


def split_sdf_blocks(text):
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conformation", required=True)
    ap.add_argument("--vs_results_dir", required=True)
    ap.add_argument("--receptor_pdb", required=True)
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    conf = args.conformation
    vs_dir = Path(args.vs_results_dir)
    passed_dir = vs_dir / f"interactions_passed_{conf}"
    poses_dir = vs_dir / f"ad4_redock_{conf}"
    scores_csv = poses_dir / "scores.csv"
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(passed_dir / f"ranked_hits_{conf}.csv") as f:
        top_ids = [row["compound_id"] for _, row in zip(range(args.n), csv.DictReader(f))]
    print(f"Top {len(top_ids)} compounds: {top_ids}")

    print(f"Looking up real paths in {scores_csv} ...")
    id_to_path = {}
    with open(scores_csv) as f:
        for row in csv.DictReader(f):
            if row["compound_id"] in top_ids:
                id_to_path[row["compound_id"]] = row["path"]
    missing = set(top_ids) - id_to_path.keys()
    if missing:
        raise SystemExit(f"Could not find these in {scores_csv}: {missing}")

    reordered = {}
    reordered_csv = passed_dir / f"reordered_compounds_{conf}.csv"
    if reordered_csv.exists():
        with open(reordered_csv) as f:
            for row in csv.DictReader(f):
                if row["compound_id"] in top_ids:
                    reordered[row["compound_id"]] = int(row["winning_pose_index_0based"])
    print(f"{len(reordered)} of the top {len(top_ids)} need pose reordering: {reordered}")

    for compound_id in top_ids:
        src = poses_dir / id_to_path[compound_id]
        if not src.is_file():
            raise SystemExit(f"{compound_id}: real path {src} doesn't exist -- scores.csv's "
                             f"path column didn't resolve, check it manually before trusting this.")
        dst = out_dir / f"{compound_id}.sdf"
        if compound_id in reordered:
            text = src.read_text()
            blocks = split_sdf_blocks(text)
            widx = reordered[compound_id]
            if "".join(blocks) != text or not (0 <= widx < len(blocks)):
                print(f"  WARNING: {compound_id} reorder verification failed, copying as-is (unreordered)")
                shutil.copy2(src, dst)
            else:
                dst.write_text("".join([blocks[widx]] + blocks[:widx] + blocks[widx + 1:]))
                print(f"  {compound_id}: reordered pose {widx} to front")
        else:
            shutil.copy2(src, dst)
            print(f"  {compound_id}: copied as-is (pose 1 was already the winner)")

    receptor_pdb = Path(args.receptor_pdb).expanduser()
    shutil.copy2(receptor_pdb, out_dir / receptor_pdb.name)

    tar_path = out_dir.parent / f"{out_dir.name}.tar.gz"
    with tarfile.open(tar_path, "w:gz") as tar:
        for f in list(out_dir.glob("*.sdf")) + [out_dir / receptor_pdb.name]:
            tar.add(f, arcname=f.name)
    print(f"Wrote {tar_path}")


if __name__ == "__main__":
    main()
