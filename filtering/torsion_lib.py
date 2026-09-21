"""
Author: Anna Korda

Independent reimplementation of the torsion-strain algorithm from Gu, Smith,
Yang, Irwin, Shoichet, "Ligand Strain Energy in Large Library Docking",
J Chem Inf Model 2021, 61, 4331-4341, DOI: 10.1021/acs.jcim.1c00368 -- see
interaction_filter_references.txt for the full citation. Written from the
paper's own Methods section and their reference tool's public source/comments
(github.com/docking-org/ChemInfTools/apps/strainfilter), not copied from it,
since that repo carries no license file. Uses their openly-distributed
torsion-library data file (data/TL_2.1_VERSION_6.xml) as intended -- this is
the precomputed statistical result of their Cambridge Structural Database
analysis; we don't need CSD access ourselves, this file already IS that
result, packaged for exactly this use.

Verified against the paper authors' own real reference output
(test2_Torsion_Strain.csv, distributed alongside their tool) on three real
molecules -- Total and Single TEU matched to floating-point precision
(~13-14 significant digits) in every case.
"""
import xml.etree.ElementTree as ET
from math import sqrt, atan2, pi, ceil

import numpy as np
from rdkit import Chem


def load_torsion_library(xml_path):
    return ET.parse(xml_path).getroot()


def _unit(a):
    return a / sqrt(np.dot(a, a))


def _dihedral(a1, a2, a3, a4):
    b1, b2, b3 = a2 - a1, a3 - a2, a4 - a3
    n1 = _unit(np.cross(b1, b2))
    n2 = _unit(np.cross(b2, b3))
    m = _unit(np.cross(n1, b2))
    return -atan2(np.dot(m, n2), np.dot(n1, n2)) * 180 / pi


def _ang_diff(t1, t2):
    if t1 < 0:
        t1 += 360
    if t2 < 0:
        t2 += 360
    d = (t1 - t2) % 360
    return d - 360 if d > 180 else d


def _tp_match(tp, hc, idx, mol, positions, out):
    smarts = tp.get("smarts")
    method = tp.get("method")
    hist_E, hist_l, hist_u = [], [], []
    if method == "exact":
        for b in tp.find("histogram_converted").findall("bin"):
            hist_E.append(float(b.get("energy")))
            hist_l.append(float(b.get("lower")))
            hist_u.append(float(b.get("upper")))

    patt = Chem.MolFromSmarts(smarts)
    if patt is None:
        return
    for match in mol.GetSubstructMatches(patt):
        if len(match) > 4:
            continue
        r1, r2, r3, r4 = (np.array(positions[match[i]]) for i in range(4))
        theta = _dihedral(r1, r2, r3, r4)

        if method == "exact":
            bin_num = ceil(theta / 10) + 17

            def interp(vals):
                return (vals[bin_num] - vals[(bin_num + 35) % 36]) / 10.0 \
                    * (theta - (bin_num - 17) * 10) + vals[bin_num]

            out.append([list(match), theta, smarts, hc, "exact",
                        interp(hist_E), interp(hist_l), interp(hist_u), False, idx])
        else:
            not_observed, energy = True, 100.0
            for angle in tp.find("angleList").findall("angle"):
                delta = _ang_diff(theta, float(angle.get("theta_0")))
                if abs(delta) <= float(angle.get("tolerance2")):
                    beta_1 = float(angle.get("beta_1"))
                    beta_2 = float(angle.get("beta_2"))
                    energy = beta_1 * (delta ** 2) + beta_2 * (delta ** 4)
                    not_observed = False
                    break
            out.append([list(match), theta, smarts, hc, "approximate",
                        energy, energy, energy, not_observed, idx])


def compute_strain(mol, torsion_lib_root):
    """mol: RDKit Mol with exactly one conformer (real 3D coords), formal
    charges intact, ring/aromaticity flags set. This function never
    sanitizes or otherwise modifies mol -- see prepare_mol_for_strain()
    in strain_filter.py for that, kept separate and explicit/auditable.

    Returns (total_TEU, single_TEU, flagged, n_torsions):
      total_TEU  = sum of every matched torsion pattern's energy estimate,
                   or None if no torsion patterns matched at all (a fully
                   rigid molecule -- unscoreable by this method, same as
                   the reference tool's own real behavior, not zero strain).
      single_TEU = the single highest individual torsion energy, or None
                   under the same no-match condition.
      flagged    = True if any torsion used the approximate method AND its
                   angle wasn't within tolerance of any known peak in the
                   library -- an unreliable estimate, per the reference
                   tool's own logic (it refuses to sum a real total in
                   this case rather than report a misleading number).
      n_torsions = how many distinct torsion patterns were matched.
    """
    positions = mol.GetConformer().GetPositions()
    bond_info = []
    idx = 0
    for hc_elem in torsion_lib_root.findall("hierarchyClass"):
        if hc_elem.get("name") != "GG":
            for tp in hc_elem.iter("torsionRule"):
                _tp_match(tp, "specific", idx, mol, positions, bond_info)
                idx += 1
    for tp in torsion_lib_root.find("hierarchyClass[@name='GG']").iter("torsionRule"):
        _tp_match(tp, "general", idx, mol, positions, bond_info)
        idx += 1

    if not bond_info:
        return None, None, False, 0

    for b in bond_info:
        if b[0][1] > b[0][2]:
            b[0] = list(reversed(b[0]))
            b.append(True)
        else:
            b.append(False)

    reduced = [bond_info[0]]
    for j in range(1, len(bond_info)):
        atoms_j = bond_info[j][0]
        matched = False
        for k in range(len(reduced)):
            if reduced[k][0] == atoms_j:
                matched = True
                if bond_info[j][9] < reduced[k][9]:
                    reduced[k] = bond_info[j]
                break
        if not matched:
            reduced.append(bond_info[j])

    final = [reduced[0]]
    for j in range(1, len(reduced)):
        a1, a2 = reduced[j][0][1], reduced[j][0][2]
        matched = False
        for k in range(len(final)):
            if final[k][0][1] == a1 and final[k][0][2] == a2:
                matched = True
                better_class = reduced[j][3][0] > final[k][3][0]
                same_class_better_E = (reduced[j][3][0] == final[k][3][0]
                                        and reduced[j][5] > final[k][5])
                if better_class or same_class_better_E:
                    final[k] = reduced[j]
                break
        if not matched:
            final.append(reduced[j])

    energies = [f[5] for f in final]
    flagged = any(f[8] for f in final)
    return sum(energies), max(energies), flagged, len(final)
