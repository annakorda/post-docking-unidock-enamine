"""
ref1 is prepared completely differently from c1/c5 -- no merged ACE/NME cap
atoms anywhere (confirmed by direct inspection: no CAY/CY/OY/CAT/NT atom names
present at all). Its one real gap (228 -> 324, LYS to LYS, ~96 residues, a
disordered loop -- same ICL3-type truncation idea as c1/c5's gap) already has
plain charged-terminus-style atoms at both flanking residues (LYS228 has HXT
like a C-term, LYS324 has H1/H2 like an N-term) -- gmx pdb2gmx handles charged
termini natively via its own residue templates when run with -ignh (discards
input H's, rebuilds correct ones for whichever terminus it applies). The only
real fix needed is a missing TER record between the two fragments -- without
it pdb2gmx would try to bond LYS228 directly to LYS324 as if they were
sequential, which they are not.

Never touches atom names/coordinates/numbering -- inserts one TER line only.
"""
import sys


def main(infile, outfile, gap_before_resnum):
    lines = open(infile).readlines()
    out = []
    inserted = False
    for i, line in enumerate(lines):
        out.append(line)
        if line.startswith("ATOM") or line.startswith("HETATM"):
            resnum = int(line[22:26])
            next_line = lines[i + 1] if i + 1 < len(lines) else ""
            next_resnum = int(next_line[22:26]) if next_line.startswith(("ATOM", "HETATM")) else None
            if resnum == gap_before_resnum and next_resnum is not None and next_resnum != gap_before_resnum:
                out.append(f"TER   {'':>5}\n")
                inserted = True

    if not inserted:
        raise SystemExit(f"Never found the boundary at resSeq {gap_before_resnum} -- refusing to write, check input")

    with open(outfile, "w") as f:
        f.writelines(out)
    print(f"Inserted TER after resSeq {gap_before_resnum}, wrote {outfile}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
