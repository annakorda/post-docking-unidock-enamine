#!/usr/bin/env python3
"""
Author: Anna Korda

Final MM-GBSA input list, Tc=0.60:
- clusters: trustworthy (Q3, std_z<1.0) AND good (Q4, raw p<0.05 mean-vs-
  pool z-test) -- medoid (real fingerprint centroid, highest mean intra-
  cluster Tanimoto) + best_scorer per cluster (medoid+best-only, no
  second_best_scorer -- dropped to keep the final list a workable size),
  deduplicated where roles coincide on the same compound.
- singletons: z<=-3.0 (q_singletons/qualifying_singletons.csv), included
  as-is, role=singleton_outlier -- unvalidated (no second member to agree
  with), not blended in as equally trustworthy.

Outputs:
- q5_mmgbsa_<conf>.sdf, one per conformation: real molecule blocks (text-
  spliced from the dedup SDF, not reparsed/rewritten) with CLUSTER_ID,
  ROLE, AD4_SCORE, SOURCE property tags appended.
- q5_mmgbsa_input_list.csv: same info, compound-level, all conformations.
- q5_kept_clusters_summary.csv: cluster-level info for only the 1,599
  kept (trustworthy+good) clusters -- not per-compound.
"""
import csv
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import BulkTanimotoSimilarity

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).parent
NESTED_DIR = HERE.parent
Q4_PATH = NESTED_DIR / "q4" / "q4_per_cluster.csv"
SINGLETONS_PATH = NESTED_DIR / "q_singletons" / "qualifying_singletons.csv"
DEDUP_DIR = NESTED_DIR.parent / "vs_results" / "dedup"  # repo_root/vs_results/dedup, dedup_stereoisomers.py's --out_dir
CONFS = ["c1", "ref1", "c5"]
TC = "0.60"


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


def add_sdf_properties(block, props):
    idx = block.rfind("$$$$")
    body = block[:idx]
    if not body.endswith("\n"):
        body += "\n"
    tag_text = "".join(f">  <{k}>\n{v}\n\n" for k, v in props.items())
    return body + tag_text + "$$$$\n"


def load_good_clusters(conf):
    with open(Q4_PATH) as f:
        return {r["cluster_id"]: r for r in csv.DictReader(f)
                if r["conformation"] == conf and r["tc"] == TC and r["good"] == "True"}


def load_cluster_members(conf, cluster_ids_wanted):
    path = NESTED_DIR / f"nested_cluster_assignment_{conf}_{TC}.csv"
    members = {}
    singleton_cluster_id = {}
    with open(path) as f:
        rows = list(csv.DictReader(f))
    by_cluster = {}
    for r in rows:
        by_cluster.setdefault(r["cluster_id"], []).append(r["compound_id"])
    for clid, ids in by_cluster.items():
        if clid in cluster_ids_wanted:
            members[clid] = ids
        if len(ids) == 1:
            singleton_cluster_id[ids[0]] = clid
    return members, singleton_cluster_id


def load_scores_and_blocks(conf, needed_ids):
    stem = f"final_hits_{conf}_cutoff-7_dedup"
    hit_dir = DEDUP_DIR / stem
    with open(hit_dir / f"{stem}.csv") as f:
        rows = list(csv.DictReader(f))
    ids = [r["compound_id"] for r in rows]
    scores = {r["compound_id"]: float(r["ad4_score"]) for r in rows}
    blocks = _split_sdf_blocks((hit_dir / f"{stem}.sdf").read_text())
    id_to_block = {cid: blk for cid, blk in zip(ids, blocks) if cid in needed_ids}
    return scores, id_to_block


def find_medoid(compound_ids, id_to_block, gen):
    fps, valid_ids = [], []
    for cid in compound_ids:
        mol = Chem.MolFromMolBlock(id_to_block[cid])
        if mol is None:
            continue
        fps.append(gen.GetFingerprint(mol))
        valid_ids.append(cid)
    if len(valid_ids) == 1:
        return valid_ids[0]
    if not valid_ids:
        return None
    best_id, best_mean_sim = None, -1.0
    for i, cid in enumerate(valid_ids):
        others = fps[:i] + fps[i + 1:]
        sims = BulkTanimotoSimilarity(fps[i], others)
        mean_sim = sum(sims) / len(sims)
        if mean_sim > best_mean_sim:
            best_mean_sim = mean_sim
            best_id = cid
    return best_id


def main():
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    all_compound_rows = []
    kept_cluster_rows = []

    for conf in CONFS:
        good_clusters = load_good_clusters(conf)
        members, singleton_cluster_id = load_cluster_members(conf, set(good_clusters))
        needed_ids = {cid for ids in members.values() for cid in ids}
        with open(SINGLETONS_PATH) as f:
            singleton_ids_this_conf = [r["compound_id"] for r in csv.DictReader(f) if r["conformation"] == conf]
        needed_ids |= set(singleton_ids_this_conf)
        scores, id_to_block = load_scores_and_blocks(conf, needed_ids)

        conf_rows = []
        for cluster_id, compound_ids in members.items():
            sorted_by_score = sorted(compound_ids, key=lambda cid: scores[cid])
            medoid = find_medoid(compound_ids, id_to_block, gen)

            roles = {}
            if len(sorted_by_score) >= 1:
                roles.setdefault(sorted_by_score[0], []).append("best_scorer")
            if medoid is not None:
                roles.setdefault(medoid, []).append("medoid")

            for cid, role_list in roles.items():
                conf_rows.append({
                    "conformation": conf, "source": "cluster", "cluster_id": cluster_id,
                    "compound_id": cid, "ad4_score": f"{scores[cid]:.3f}",
                    "cluster_size": len(compound_ids), "roles": "+".join(sorted(role_list)),
                })

            cinfo = good_clusters[cluster_id]
            kept_cluster_rows.append({
                "conformation": conf, "cluster_id": cluster_id, "count": cinfo["count"],
                "mean_z": cinfo["mean_z"], "min_z": cinfo["min_z"], "std_z": cinfo["std_z"],
                "t_stat": cinfo["t_stat"], "p_value": cinfo["p_value"],
            })

        for cid in singleton_ids_this_conf:
            conf_rows.append({
                "conformation": conf, "source": "singleton", "cluster_id": singleton_cluster_id[cid],
                "compound_id": cid, "ad4_score": f"{scores[cid]:.3f}",
                "cluster_size": 1, "roles": "singleton_outlier",
            })

        # per-conformation SDF: real molblocks + property tags
        sdf_path = HERE / f"q5_mmgbsa_{conf}.sdf"
        with open(sdf_path, "w") as f:
            for r in conf_rows:
                block = id_to_block[r["compound_id"]]
                block = add_sdf_properties(block, {
                    "COMPOUND_ID": r["compound_id"], "CONFORMATION": conf,
                    "CLUSTER_ID": r["cluster_id"], "ROLE": r["roles"],
                    "AD4_SCORE": r["ad4_score"], "SOURCE": r["source"],
                })
                f.write(block)
        print(f"{conf}: {len(good_clusters):,} kept clusters + {len(singleton_ids_this_conf):,} singletons "
              f"-> {len(conf_rows):,} compounds -> wrote {sdf_path.name}")

        all_compound_rows.extend(conf_rows)

    csv_path = HERE / "q5_mmgbsa_input_list.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["conformation", "source", "cluster_id", "compound_id",
                                           "ad4_score", "cluster_size", "roles"])
        w.writeheader()
        w.writerows(all_compound_rows)
    print(f"\nWrote {csv_path} ({len(all_compound_rows):,} rows)")

    clusters_path = HERE / "q5_kept_clusters_summary.csv"
    with open(clusters_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["conformation", "cluster_id", "count", "mean_z", "min_z",
                                           "std_z", "t_stat", "p_value"])
        w.writeheader()
        w.writerows(kept_cluster_rows)
    print(f"Wrote {clusters_path} ({len(kept_cluster_rows):,} rows, kept clusters only)")

    print(f"\n{'conf':<6}{'cluster_picks':>15}{'singletons':>12}{'total':>8}")
    for conf in CONFS:
        c = sum(1 for r in all_compound_rows if r["conformation"] == conf and r["source"] == "cluster")
        s = sum(1 for r in all_compound_rows if r["conformation"] == conf and r["source"] == "singleton")
        print(f"{conf:<6}{c:>15,}{s:>12,}{c+s:>8,}")


if __name__ == "__main__":
    main()
