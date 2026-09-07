"""
Generalizes fix_caps.py to all 4 real cap-merge sites in this receptor, found
by direct inspection of c1.pdb (never assumed):

  - resSeq 24  (ASN): real N-terminus, ACE cap merged in     -- same as fix_caps.py
  - resSeq 245 (SER): end of fragment 1, NME cap merged in   -- NEW: different
                       atom names than the C-term case (HT/HAT1-3, not HNT/HT1-3)
  - resSeq 316 (CYS): start of fragment 2, ACE cap merged in -- NEW
  - resSeq 422 (GLN): real C-terminus, NME cap merged in     -- same as fix_caps.py

245->316 is pdb4amber's one genuine gap (70 residues, ~70A away from Asp116 in
the real ligand pocket -- confirmed by direct distance measurement, this is an
intracellular loop truncation, not something to reconstruct). Maestro already
capped both ends of both fragments correctly for this construct; the only
problem is that GROMACS/AMBER expect ACE/NME as separate residues with their
own standard atom names, not merged into the flanking real residue under a
CHARMM-style naming scheme -- exactly what fix_caps.py already solved for the
2 true termini, extended here to the 2 additional fragment-break caps.

Each ACE site also carries 2 stray atoms literally named "HY" with a blank
segid column -- present at BOTH real ACE sites (ASN24 and CYS316), so this is
a consistent Maestro export quirk, not corruption -- dropped, matching
fix_caps.py's original behavior at the N-terminus.

Never touches the original file -- writes a new one, and inserts one TER
record between the two fragments (they are genuinely 2 separate, disconnected
polypeptide chains once released from their artificial merged-cap naming).
"""
import sys

ACE_RENAME = {"CAY": ("CH3", "C"), "HY1": ("HH31", "H"), "HY2": ("HH32", "H"),
              "HY3": ("HH33", "H"), "CY": ("C", "C"), "OY": ("O", "O"),
              # HAY1-3 variant seen in a later Maestro re-export of the same construct
              "HAY1": ("HH31", "H"), "HAY2": ("HH32", "H"), "HAY3": ("HH33", "H")}
NME_RENAME = {"NT": ("N", "N"), "CAT": ("CH3", "C"),
              "HNT": ("H", "H"), "HT": ("H", "H"),
              "HT1": ("HH31", "H"), "HT2": ("HH32", "H"), "HT3": ("HH33", "H"),
              "HAT1": ("HH31", "H"), "HAT2": ("HH32", "H"), "HAT3": ("HH33", "H")}

# (residue carrying the merged cap atoms, new resSeq to give the split-out cap,
#  which fragment it belongs to -- for TER placement)
ACE_SITES = [(24, 23, 1), (316, 315, 2)]
NME_SITES = [(245, 246, 1), (422, 423, 2)]


def parse_atom_line(line):
    return {"name": line[12:16].strip(), "resnum": int(line[22:26])}


def fmt_atom(serial, name, resname, chain, resnum, x, y, z, occ, temp, element):
    name_field = f" {name:<3}" if len(name) < 4 else name[:4]
    return (f"ATOM  {serial:>5} {name_field} {resname:<3} {chain}{resnum:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{temp:6.2f}          {element:>2}\n")


def main(infile, outfile):
    lines = open(infile).readlines()
    ace_by_site = {r: [] for r, _, _ in ACE_SITES}
    nme_by_site = {r: [] for r, _, _ in NME_SITES}
    dropped = []
    header = []
    chain = "X"
    body = {}  # resnum -> list of (serial-agnostic) formatted-or-raw lines, real atoms only

    for line in lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            # CONECT records reference atom serial numbers, which we're about to
            # renumber wholesale -- keeping them would leave stale/wrong bonds
            # in the output. gmx pdb2gmx doesn't read CONECT anyway (topology
            # comes from its own residue templates + sequence order), so just
            # drop them rather than write something misleading.
            if not line.startswith(("TER", "END", "CONECT")):
                header.append(line)
            continue
        info = parse_atom_line(line)
        resnum = info["resnum"]
        name = info["name"]
        chain = line[21]
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        occ = float(line[54:60]) if line[54:60].strip() else 1.0
        temp = float(line[60:66]) if line[60:66].strip() else 0.0

        if resnum in ace_by_site and name in ACE_RENAME:
            new_name, elem = ACE_RENAME[name]
            ace_by_site[resnum].append((new_name, x, y, z, occ, temp, elem))
        elif resnum in ace_by_site and name == "HY":
            dropped.append((resnum, name))
        elif resnum in nme_by_site and name in NME_RENAME:
            new_name, elem = NME_RENAME[name]
            nme_by_site[resnum].append((new_name, x, y, z, occ, temp, elem))
        else:
            body.setdefault(resnum, []).append(line)

    print(f"Dropped (stray, not part of the standard 6-atom ACE cap): {dropped}")
    for r, new_r, frag in ACE_SITES:
        n = len(ace_by_site[r])
        flag = "" if n == 6 else "  <-- WARNING: expected 6, got a naming variant not yet in ACE_RENAME"
        print(f"ACE at old-resSeq {r} -> new residue {new_r} ({n} atoms), fragment {frag}{flag}")
    for r, new_r, frag in NME_SITES:
        n = len(nme_by_site[r])
        flag = "" if n == 6 else "  <-- WARNING: expected 6, got a naming variant not yet in NME_RENAME"
        print(f"NME at old-resSeq {r} -> new residue {new_r} ({n} atoms), fragment {frag}{flag}")

    ace_newnum = {r: n for r, n, _ in ACE_SITES}
    nme_newnum = {r: n for r, n, _ in NME_SITES}
    # fragment for each of the 4 new cap-residue numbers comes from its own site
    # definition (which residue it caps), NOT from a resnum threshold -- the new
    # NME at 246 closes fragment 1, the new ACE at 315 opens fragment 2, even
    # though 246 > 245 and 315 < 316
    newnum_frag = {n: f for _, n, f in ACE_SITES + NME_SITES}

    all_resnums = sorted(set(body) | set(ace_newnum.values()) | set(nme_newnum.values()))
    out = list(header)
    serial = 1
    prev_frag = None

    for resnum in all_resnums:
        frag = newnum_frag.get(resnum, 1 if resnum <= 245 else 2)
        if prev_frag is not None and frag != prev_frag:
            out.append(f"TER   {serial:>5}\n")
        prev_frag = frag

        matched_ace = next((r for r, n, _ in ACE_SITES if n == resnum), None)
        matched_nme = next((r for r, n, _ in NME_SITES if n == resnum), None)

        if matched_ace is not None:
            for name, x, y, z, occ, temp, elem in ace_by_site[matched_ace]:
                out.append(fmt_atom(serial, name, "ACE", chain, resnum, x, y, z, occ, temp, elem))
                serial += 1
        elif matched_nme is not None:
            for name, x, y, z, occ, temp, elem in nme_by_site[matched_nme]:
                out.append(fmt_atom(serial, name, "NME", chain, resnum, x, y, z, occ, temp, elem))
                serial += 1
        else:
            for line in body[resnum]:
                out.append(line[:6] + f"{serial:>5}" + line[11:])
                serial += 1

    out.append(f"TER   {serial:>5}\n")
    out.append("END\n")

    with open(outfile, "w") as f:
        f.writelines(out)
    print(f"Wrote {outfile}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
