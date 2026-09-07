# receptor_prep

Author: Anna Korda

Fixes the real blocker on `gmx pdb2gmx` for the Uni-GBSA receptor prep --
verified by direct inspection, not assumed. What pdb4amber originally reported
as "missing heavy atoms" and a "structural gap" for c1/c5 turned out to be
mostly false positives: the only genuine gap is one ~70-residue disordered
loop (a truncated ICL3, ~70A from Asp116 -- nowhere near the ligand pocket).
Every one of the 9 "missing atom" residues pdb4amber flagged has its full,
real atom set present when checked directly against the raw file; the actual
problem was that Maestro capped both ends of the truncation with ACE/NME
groups **merged into the flanking residue** under a CHARMM-style naming
scheme GROMACS/AMBER don't recognize -- the same class of issue as the true
N/C-termini, just at 2 more locations.

## c1.pdb / c5.pdb -> fix_caps.py

Both are the same construct: real termini at resSeq 24 (N, ACE-capped) and
422 (C, NME-capped), one internal gap at 245->316 (SER capped NME-style /
CYS capped ACE-style on either side). `fix_caps.py` splits all 4 merged caps
into proper standalone ACE/NME residues with AMBER atom names, drops 2 stray
always-present junk atoms per ACE site (literally named "HY", blank segid --
present in every real file seen so far), and inserts a TER between the two
now-genuinely-disconnected fragments.

Handles 2 real ACE-methyl-hydrogen naming variants seen across different
Maestro export runs of the same construct (`HY1/2/3` and `HAY1/2/3`) and 2
NME amide-hydrogen variants (`HNT`/`HT1-3` and `HT`/`HAT1-3`) -- if a future
receptor file uses a third variant, the script prints an explicit WARNING
(atom count != 6 at a cap site) instead of silently writing a broken cap.

```
python fix_caps.py c1.pdb c1_gmxready.pdb
python fix_caps.py c5.pdb c5_gmxready.pdb
```

## ref1.pdb -> fix_gap_ref1.py

A genuinely different prep (real experimental cryo-EM structure, 7E2Y, not a
Maestro-minimized construct) -- confirmed by direct inspection to have **no**
merged cap atoms anywhere. Its one gap (228->324, ~96 residues) already has
plain charged-terminus-style atoms on both flanking LYS residues; the only
real problem is a missing TER record, which is all this script inserts.
`gmx pdb2gmx -ignh` then rebuilds correct termini automatically, no prompt
needed.

```
python fix_gap_ref1.py ref1.pdb ref1_pre_ter.pdb 228
# (ref1_gmxready.pdb here also has the TER at 228 already applied)
```

## Verified

`gmx pdb2gmx -water tip3p -ff amber03 -ignh` runs clean (topology written,
no errors) on all 3 of `c1_gmxready.pdb`, `c5_gmxready.pdb`,
`ref1_gmxready.pdb` as committed here. This only unblocks receptor prep --
the rest of the Uni-GBSA pipeline (ligand topology, solvation, minimization,
actual MM-GBSA calculation) is still open.
