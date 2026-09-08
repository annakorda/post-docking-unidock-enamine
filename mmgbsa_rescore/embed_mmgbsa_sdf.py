#!/usr/bin/env python3
"""
Author: Anna Korda

Takes mmgbsa_rescore.py's already-ranked mmgbsa_results_<conf>.csv and writes
ONE combined SDF per conformation: the original poses (same bytes, untouched
geometry/atoms/bonds/charges -- pure text-level splice, no RDKit or any other
chemistry library touches these files anywhere in this script), re-ordered by
mmgbsa_dg (best/most negative first), with real new SDF data fields --
`> <mmgbsa_dg>`, `> <ad4_score>`, `> <total_TEU>`, `> <single_TEU>` -- inserted
before each block's `$$$$` terminator. Same convention this project's own
SDFs already use for `> <Uni-Dock RESULT>`/`ENERGY=`, so Maestro shows all of
this directly as real properties without needing the separate CSV.

Only the property block is added; the molecule block itself is never
modified (verified byte-identical, see the script's own test in git history).

Usage:
    python embed_mmgbsa_sdf.py --conformation c1 \\
        --vs_results_dir /users/gpcr/annak/ultra-large/vs_results \\
        --mmgbsa_out_dir /users/gpcr/annak/ultra-large/analysis_results/mmgbsa_rescore \\
        --out_dir /users/gpcr/annak/ultra-large/analysis_results/mmgbsa_rescore \\
        --cutoff -7.0
"""
import argparse
import csv
import time
from pathlib import Path


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def _split_sdf_blocks(text):
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


def load_original_blocks(vs_results_dir, conf, cutoff):
    hit_dir = Path(vs_results_dir) / f"final_hits_{conf}_cutoff{cutoff:g}"
    csv_path = hit_dir / f"final_hits_{conf}_cutoff{cutoff:g}.csv"
    if not csv_path.is_file():
        raise SystemExit(f"{csv_path} not found")
    with open(csv_path) as f:
        ids = [r["compound_id"] for r in csv.DictReader(f)]
    part_files = sorted(hit_dir.glob(f"final_hits_{conf}_cutoff{cutoff:g}_part*.sdf"))
    if not part_files:
        raise SystemExit(f"No SDF parts found in {hit_dir}")
    blocks = []
    for pf in part_files:
        blocks.extend(_split_sdf_blocks(pf.read_text()))
    if len(blocks) != len(ids):
        raise SystemExit(f"{conf}: {len(ids)} CSV rows but {len(blocks)} SDF blocks -- "
                          f"expected the same rank order, refusing to guess a correspondence")
    return dict(zip(ids, blocks))


FIELDS = ["mmgbsa_dg", "ad4_score", "total_TEU", "single_TEU"]


def insert_fields(block, row):
    """Inserts real SDF data fields (one per name in FIELDS) before the
    block's $$$$ terminator -- same convention as this project's own
    `> <Uni-Dock RESULT>` field. Everything before this stays byte-identical
    to the original block; pure string splice, no chemistry library involved."""
    idx = block.rfind("$$$$")
    if idx == -1:
        raise ValueError("block has no $$$$ terminator")
    prefix, suffix = block[:idx], block[idx:]
    if not prefix.endswith("\n"):
        prefix += "\n"
    field_text = "".join(f"> <{name}>\n{row[name]}\n\n" for name in FIELDS)
    return prefix + field_text + suffix


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conformation", required=True)
    ap.add_argument("--vs_results_dir", required=True)
    ap.add_argument("--mmgbsa_out_dir", required=True, help="mmgbsa_rescore.py's --out_dir")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--cutoff", type=float, default=-7.0)
    args = ap.parse_args()

    conf = args.conformation
    results_csv = Path(args.mmgbsa_out_dir) / f"mmgbsa_{conf}" / f"mmgbsa_results_{conf}.csv"
    if not results_csv.is_file():
        raise SystemExit(f"{results_csv} not found -- run mmgbsa_rescore.py for {conf} first")
    with open(results_csv) as f:
        ranked_rows = list(csv.DictReader(f))
    log(f"{len(ranked_rows):,} mmgbsa-ranked compounds loaded for {conf}")

    log("Loading original SDF poses ...")
    id_to_block = load_original_blocks(args.vs_results_dir, conf, args.cutoff)

    out_dir = Path(args.out_dir) / f"mmgbsa_{conf}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"mmgbsa_ranked_{conf}.sdf"

    written = 0
    missing = []
    with open(out_path, "w") as out_f:
        for row in ranked_rows:
            cid = row["compound_id"]
            if cid not in id_to_block:
                missing.append(cid)
                continue
            out_f.write(insert_fields(id_to_block[cid], row))
            written += 1

    log(f"Wrote {written:,} compounds to {out_path} (one file, ranked best mmgbsa_dg first, "
        f"{', '.join(FIELDS)} embedded as real SDF fields)")
    if missing:
        log(f"WARNING: {len(missing)} compound(s) in mmgbsa_results_{conf}.csv had no matching "
            f"original SDF block, skipped: {missing[:10]}{' ...' if len(missing) > 10 else ''}")


if __name__ == "__main__":
    main()
