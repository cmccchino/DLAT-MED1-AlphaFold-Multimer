"""
dlat_med1_interface_analysis.py
================================
Atomic-level interface analysis of DLAT–MED1 predicted complexes
from AlphaFold-Multimer v3 PDB output files.

Project: Computational investigation of the structural basis of
         DLAT–Mediator complex interaction (Russo et al., 2024, Mol Cell)

Methods:
    - Interface residues defined as residues containing at least one atom
      within 5.0 Angstroms of any atom from the opposing chain
    - Hydrogen bonds: N/O donor-acceptor pairs between 2.5 and 3.5 Angstroms
    - Salt bridges: oppositely charged residue pairs (Arg/Lys/His vs Asp/Glu)
      with charged atom groups within 4.0 Angstroms
    - Hydrophobic contacts: nonpolar residue pairs (Ala/Val/Ile/Leu/Met/
      Phe/Trp/Pro/Tyr) within 5.0 Angstroms

Sequence numbering:
    - DLAT: mature sequence numbering (residue 1 = UniProt P10515 residue 87)
    - MED1 TED fragment: PDB residue + 71 = UniProt Q15648 full sequence
    - MED1 LXXLL fragment: PDB residue + 589 = UniProt Q15648 full sequence

Requirements: Python 3.x standard library only (math, os, sys, collections)

Usage:
    python dlat_med1_interface_analysis.py <pdb_file> <dlat_chains> <med1_chain> <med1_offset>

    Example (Run 1 - DLAT monomer + MED1 TED):
        python dlat_med1_interface_analysis.py run1.pdb A B 71

    Example (Run 3 - DLAT trimer + MED1 LXXLL):
        python dlat_med1_interface_analysis.py run3.pdb A,B,C D 589

Author: Computed as part of a self-initiated computational mini-project
        for the Russo Laboratory application, LMU Munich, May 2026
"""

import math
import sys
import os
from collections import defaultdict


# ============================================================
# CONSTANTS — Amino acid classifications
# ============================================================

HYDROPHOBIC_RES = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'PRO', 'TYR'}

POSITIVE_RES    = {'ARG', 'LYS', 'HIS'}
NEGATIVE_RES    = {'ASP', 'GLU'}

POSITIVE_ATOMS  = {'NZ', 'NH1', 'NH2', 'NE', 'ND1', 'NE2'}
NEGATIVE_ATOMS  = {'OD1', 'OD2', 'OE1', 'OE2'}

HBOND_ELEMENTS  = {'N', 'O'}


# ============================================================
# DLAT DOMAIN BOUNDARIES — mature sequence numbering
# ============================================================

DLAT_DOMAINS = [
    (1,   90,  "Lipoyl domain 1"),
    (91,  131, "Linker between lipoyl domains"),
    (132, 208, "Lipoyl domain 2"),
    (209, 269, "Inter-domain linker"),
    (270, 307, "PSBD (Peripheral Subunit-Binding Domain)"),
    (308, 333, "PSBD-AT linker"),
    (334, 561, "Acetyltransferase catalytic domain"),
]

# MED1 TED fragment domain boundaries (fragment numbering, offset +71 for full seq)
MED1_TED_DOMAINS = [
    (1,   356, "InterPro structured domain (full seq 72-427)"),
    (357, 449, "Extended structured region (full seq 428-520)"),
]

# MED1 LXXLL fragment domain boundaries (fragment numbering, offset +589 for full seq)
MED1_LXXLL_DOMAINS = [
    (1,   14,  "Pre-LXXLL1 region"),
    (15,  20,  "LXXLL motif 1 (full seq 604-608: LTSLL)"),
    (21,  55,  "Between LXXLL motifs"),
    (56,  60,  "LXXLL motif 2 (full seq 645-649: LMNLL)"),
    (61,  141, "Post-LXXLL2 / GATA1 interaction region"),
]


# ============================================================
# PDB PARSING
# ============================================================

def parse_pdb(pdb_file):
    """
    Parse a PDB file and return a list of atom dictionaries.

    Parameters
    ----------
    pdb_file : str
        Path to the PDB file.

    Returns
    -------
    list of dict
        Each dict contains: chain, resi (residue number), resn (residue name),
        atom (atom name), x, y, z (coordinates), elem (element).
    """
    atoms = []
    if not os.path.exists(pdb_file):
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    with open(pdb_file, 'r') as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            try:
                atom = {
                    'chain': line[21],
                    'resi':  int(line[22:26].strip()),
                    'resn':  line[17:20].strip(),
                    'atom':  line[12:16].strip(),
                    'x':     float(line[30:38]),
                    'y':     float(line[38:46]),
                    'z':     float(line[46:54]),
                    'elem':  (line[76:78].strip()
                              if len(line) > 76 and line[76:78].strip()
                              else line[12:14].strip())[0].upper()
                }
                atoms.append(atom)
            except (ValueError, IndexError):
                continue  # skip malformed lines

    print(f"  Parsed {len(atoms)} atoms from {os.path.basename(pdb_file)}")
    return atoms


def get_chain_summary(atoms):
    """Return a dict of chain -> (min_resi, max_resi, n_residues)."""
    chains = defaultdict(set)
    for a in atoms:
        chains[a['chain']].add(a['resi'])
    return {c: (min(r), max(r), len(r)) for c, r in chains.items()}


# ============================================================
# DISTANCE CALCULATION
# ============================================================

def euclidean_distance(a1, a2):
    """Calculate Euclidean distance between two atoms."""
    return math.sqrt(
        (a1['x'] - a2['x'])**2 +
        (a1['y'] - a2['y'])**2 +
        (a1['z'] - a2['z'])**2
    )


# ============================================================
# INTERFACE DETECTION
# ============================================================

def find_interface(atoms, chain1, chain2, cutoff=5.0):
    """
    Find all residues at the interface between two chains.

    Parameters
    ----------
    atoms : list of dict
        Parsed atom list from parse_pdb().
    chain1 : str
        Chain identifier for first protein (e.g. 'A').
    chain2 : str
        Chain identifier for second protein (e.g. 'B').
    cutoff : float
        Distance cutoff in Angstroms (default: 5.0A).

    Returns
    -------
    iface1 : list of (resi, resn) tuples — interface residues in chain1
    iface2 : list of (resi, resn) tuples — interface residues in chain2
    contacts : list of (atom1, atom2, distance) tuples — all atomic contacts
    """
    atoms1 = [a for a in atoms if a['chain'] == chain1]
    atoms2 = [a for a in atoms if a['chain'] == chain2]

    iface1, iface2 = set(), set()
    contacts = []

    for a1 in atoms1:
        for a2 in atoms2:
            d = euclidean_distance(a1, a2)
            if d <= cutoff:
                iface1.add((a1['resi'], a1['resn']))
                iface2.add((a2['resi'], a2['resn']))
                contacts.append((a1, a2, d))

    return sorted(iface1), sorted(iface2), contacts


# ============================================================
# DOMAIN MAPPING
# ============================================================

def get_dlat_domain(resi):
    """Return the DLAT domain name for a given mature sequence residue number."""
    for start, end, name in DLAT_DOMAINS:
        if start <= resi <= end:
            return name
    return "Unknown"


def get_med1_domain(resi, offset):
    """
    Return the MED1 domain name for a given fragment residue number.

    Parameters
    ----------
    resi : int
        Residue number in fragment (PDB) numbering.
    offset : int
        71 for TED fragment, 589 for LXXLL fragment.
    """
    domains = MED1_TED_DOMAINS if offset == 71 else MED1_LXXLL_DOMAINS
    for start, end, name in domains:
        if start <= resi <= end:
            return name
    return "Unknown"


# ============================================================
# PISA-EQUIVALENT ANALYSIS
# ============================================================

def analyze_hbonds(contacts, med1_offset):
    """
    Identify hydrogen bonds between chains.
    Criterion: N or O donor/acceptor pairs within 2.5-3.5 Angstroms.

    Returns list of (atom1, atom2, distance) sorted by distance.
    """
    seen = set()
    hbonds = []

    for a1, a2, d in contacts:
        if a1['elem'] in HBOND_ELEMENTS and a2['elem'] in HBOND_ELEMENTS:
            if 2.5 <= d <= 3.5:
                key = (a1['chain'], a1['resi'], a1['atom'],
                       a2['resi'], a2['atom'])
                if key not in seen:
                    seen.add(key)
                    hbonds.append((a1, a2, d))

    return sorted(hbonds, key=lambda x: x[2])


def analyze_salt_bridges(contacts, med1_offset):
    """
    Identify salt bridges between chains.
    Criterion: oppositely charged residue pairs with charged atoms within 4.0A.

    Returns list of (atom1, atom2, distance, type_str) sorted by distance.
    """
    seen = set()
    sbs = []

    for a1, a2, d in contacts:
        if d > 4.0:
            continue

        # DLAT positive — MED1 negative
        if (a1['resn'] in POSITIVE_RES and a2['resn'] in NEGATIVE_RES and
                a1['atom'] in POSITIVE_ATOMS and a2['atom'] in NEGATIVE_ATOMS):
            key = (a1['chain'], a1['resi'], a1['resn'], a2['resi'], a2['resn'])
            if key not in seen:
                seen.add(key)
                sbs.append((a1, a2, d, 'DLAT(+) — MED1(-)'))

        # DLAT negative — MED1 positive
        if (a1['resn'] in NEGATIVE_RES and a2['resn'] in POSITIVE_RES and
                a1['atom'] in NEGATIVE_ATOMS and a2['atom'] in POSITIVE_ATOMS):
            key = (a1['chain'], a1['resi'], a1['resn'], a2['resi'], a2['resn'])
            if key not in seen:
                seen.add(key)
                sbs.append((a1, a2, d, 'DLAT(-) — MED1(+)'))

    return sorted(sbs, key=lambda x: x[2])


def analyze_hydrophobic(contacts, med1_offset):
    """
    Identify hydrophobic contacts between chains.
    Criterion: both residues are nonpolar and within 5.0A.

    Returns list of (atom1, atom2, distance) deduplicated at residue level,
    sorted by distance.
    """
    seen = set()
    hydro = []

    for a1, a2, d in contacts:
        if a1['resn'] in HYDROPHOBIC_RES and a2['resn'] in HYDROPHOBIC_RES:
            key = (a1['chain'], a1['resi'], a1['resn'], a2['resi'], a2['resn'])
            if key not in seen:
                seen.add(key)
                hydro.append((a1, a2, d))

    return sorted(hydro, key=lambda x: x[2])


# ============================================================
# DOMAIN SUMMARY
# ============================================================

def summarize_dlat_domains(iface_residues):
    """Count interface residues per DLAT domain."""
    domain_counts = defaultdict(list)
    for resi, resn in iface_residues:
        dom = get_dlat_domain(resi)
        domain_counts[dom].append((resi, resn))
    return domain_counts


def summarize_med1_domains(iface_residues, offset):
    """Count interface residues per MED1 domain."""
    domain_counts = defaultdict(list)
    for resi, resn in iface_residues:
        dom = get_med1_domain(resi, offset)
        domain_counts[dom].append((resi, resn, resi + offset))
    return domain_counts


# ============================================================
# REPORTING
# ============================================================

def print_separator(char='=', width=70):
    print(char * width)


def print_section(title):
    print_separator()
    print(f"  {title}")
    print_separator()


def report_run(pdb_file, dlat_chains, med1_chain, med1_offset, run_name):
    """
    Run complete interface analysis for one AlphaFold-Multimer prediction.

    Parameters
    ----------
    pdb_file : str
        Path to the unrelaxed rank 1 PDB file from ColabFold.
    dlat_chains : list of str
        Chain identifiers for DLAT (e.g. ['A'] for monomer, ['A','B','C'] for trimer).
    med1_chain : str
        Chain identifier for MED1 fragment (e.g. 'B' or 'D').
    med1_offset : int
        71 for TED fragment, 589 for LXXLL fragment.
    run_name : str
        Descriptive name for this run.
    """
    print(f"\n\n{'='*70}")
    print(f"  {run_name}")
    print(f"{'='*70}")

    # Parse
    atoms = parse_pdb(pdb_file)
    chain_info = get_chain_summary(atoms)
    print("\nChain summary:")
    for chain, (mn, mx, nr) in sorted(chain_info.items()):
        label = "MED1" if chain == med1_chain else f"DLAT chain {chain}"
        print(f"  Chain {chain} [{label}]: {nr} residues, range {mn}–{mx}")

    # Per-chain interface analysis
    all_dlat_iface = {}     # (chain, resi, resn) -> domain
    all_med1_iface = {}     # (resi, resn) -> domain
    all_contacts   = []

    for dc in dlat_chains:
        iface1, iface2, contacts = find_interface(atoms, dc, med1_chain)

        if not iface1:
            print(f"\n  Chain {dc}: no contacts with MED1 chain {med1_chain}")
            continue

        print(f"\n  ── DLAT chain {dc} → MED1 chain {med1_chain} ──")
        print(f"  DLAT interface residues: {len(iface1)}")

        # Domain breakdown — DLAT
        dlat_doms = summarize_dlat_domains(iface1)
        for dom, res in sorted(dlat_doms.items()):
            res_str = " ".join(f"{resn}{resi}" for resi, resn in sorted(res))
            print(f"    {dom} ({len(res)} res): {res_str}")

        # PSBD check — critical finding
        psbd_hits = [r for r in iface1 if 270 <= r[0] <= 307]
        if not psbd_hits:
            print(f"    *** PSBD (270-307): 0 contacts — PSBD completely free ***")

        print(f"  MED1 interface residues: {len(iface2)}")
        med1_doms = summarize_med1_domains(iface2, med1_offset)
        for dom, res in sorted(med1_doms.items()):
            res_str = " ".join(f"{resn}{resi}(full:{full})"
                               for resi, resn, full in sorted(res))
            print(f"    {dom} ({len(res)} res): {res_str}")

        # Accumulate
        for resi, resn in iface1:
            all_dlat_iface[(dc, resi, resn)] = get_dlat_domain(resi)
        for resi, resn in iface2:
            all_med1_iface[(resi, resn)] = get_med1_domain(resi, med1_offset)
        all_contacts.extend(contacts)

    # PISA-equivalent analysis
    print(f"\n  ── PISA-equivalent interface analysis ──")
    hbonds  = analyze_hbonds(all_contacts, med1_offset)
    sbridges = analyze_salt_bridges(all_contacts, med1_offset)
    hydro   = analyze_hydrophobic(all_contacts, med1_offset)

    print(f"  H-bonds (2.5–3.5A N/O pairs):      {len(hbonds)}")
    print(f"  Salt bridges (charged pairs <4.0A): {len(sbridges)}")
    print(f"  Hydrophobic contacts (<5.0A pairs): {len(hydro)}")

    # Interface character
    if len(hydro) > len(hbonds):
        character = "Hydrophobic-dominated"
    elif len(hbonds) > len(hydro):
        character = "Polar/H-bond-dominated"
    else:
        character = "Mixed"
    print(f"  Dominant interface character:       {character}")

    # Print H-bonds
    if hbonds:
        print(f"\n  Hydrogen bonds (sorted by distance):")
        print(f"  {'DLAT':20s} {'Atom':6s} | {'MED1 (frag)':15s} {'MED1 (full)':12s} {'Atom':8s} {'Dist(A)':8s}")
        print(f"  {'-'*20} {'-'*6}-+-{'-'*15} {'-'*12} {'-'*8} {'-'*8}")
        for a1, a2, d in hbonds:
            dlat_label = f"{a1['resn']}{a1['resi']} (ch {a1['chain']})"
            med1_frag  = f"{a2['resn']}{a2['resi']}"
            med1_full  = f"{a2['resi'] + med1_offset}"
            print(f"  {dlat_label:20s} {a1['atom']:6s} | {med1_frag:15s} {med1_full:12s} {a2['atom']:8s} {d:.2f}")

    # Print salt bridges
    if sbridges:
        print(f"\n  Salt bridges:")
        for a1, a2, d, typ in sbridges:
            dlat_label = f"DLAT-{a1['chain']} {a1['resn']}{a1['resi']}"
            med1_label = f"MED1 {a2['resn']}{a2['resi']} (full seq {a2['resi']+med1_offset})"
            print(f"  {dlat_label} — {med1_label}   {d:.2f}A  [{typ}]")
    else:
        print(f"\n  Salt bridges: none")

    # Print hydrophobic contacts
    if hydro:
        print(f"\n  Hydrophobic contacts (sorted by distance):")
        print(f"  {'DLAT':20s} | {'MED1 (frag)':15s} {'MED1 (full)':12s} {'Dist(A)':8s}")
        print(f"  {'-'*20}-+-{'-'*15} {'-'*12} {'-'*8}")
        for a1, a2, d in hydro:
            dlat_label = f"{a1['resn']}{a1['resi']} (ch {a1['chain']})"
            med1_frag  = f"{a2['resn']}{a2['resi']}"
            med1_full  = f"{a2['resi'] + med1_offset}"
            print(f"  {dlat_label:20s} | {med1_frag:15s} {med1_full:12s} {d:.2f}")

    # Summary counts
    print(f"\n  ── Summary ──")
    print(f"  Total DLAT interface residues: {len(all_dlat_iface)}")
    print(f"  Total MED1 interface residues: {len(all_med1_iface)}")

    return {
        'run': run_name,
        'dlat_residues': len(all_dlat_iface),
        'med1_residues': len(all_med1_iface),
        'hbonds': len(hbonds),
        'salt_bridges': len(sbridges),
        'hydrophobic': len(hydro),
        'character': character,
    }


# ============================================================
# COMMAND-LINE INTERFACE
# ============================================================

def parse_args():
    """Parse command-line arguments."""
    if len(sys.argv) < 5:
        print(__doc__)
        print("\nUsage:")
        print("  python dlat_med1_interface_analysis.py <pdb_file> <dlat_chains> <med1_chain> <med1_offset>")
        print("\nExamples:")
        print("  python dlat_med1_interface_analysis.py run1.pdb A B 71")
        print("  python dlat_med1_interface_analysis.py run3.pdb A,B,C D 589")
        sys.exit(1)

    pdb_file    = sys.argv[1]
    dlat_chains = sys.argv[2].split(',')
    med1_chain  = sys.argv[3]
    med1_offset = int(sys.argv[4])
    run_name    = sys.argv[5] if len(sys.argv) > 5 else os.path.basename(pdb_file)

    return pdb_file, dlat_chains, med1_chain, med1_offset, run_name


# ============================================================
# BATCH MODE — run all four experimental configurations
# ============================================================

def run_all_analysis(pdb_files):
    """
    Run analysis on all four runs and print a comparison summary.

    Parameters
    ----------
    pdb_files : dict
        Dictionary mapping run name to (pdb_path, dlat_chains, med1_chain, med1_offset).

    Example:
        pdb_files = {
            'Run1_monomer_TED':   ('run1.pdb',  ['A'],       'B', 71),
            'Run2_monomer_LXXLL': ('run2.pdb',  ['A'],       'B', 589),
            'Run3_trimer_LXXLL':  ('run3.pdb',  ['A','B','C'],'D', 589),
            'Run4_trimer_TED':    ('run4.pdb',  ['A','B','C'],'D', 71),
        }
    """
    results = []
    for run_name, (pdb, dlat_chains, med1_chain, offset) in pdb_files.items():
        result = report_run(pdb, dlat_chains, med1_chain, offset, run_name)
        results.append(result)

    # Comparison summary table
    print(f"\n\n{'='*70}")
    print("  COMPARISON SUMMARY — ALL RUNS")
    print(f"{'='*70}")
    print(f"  {'Run':<30} {'DLAT':>6} {'MED1':>6} {'H-bonds':>8} {'SaltBr':>8} {'Hydrophob':>10}  Character")
    print(f"  {'-'*30} {'-'*6} {'-'*6} {'-'*8} {'-'*8} {'-'*10}  ---------")
    for r in results:
        name = r['run'][:30]
        print(f"  {name:<30} {r['dlat_residues']:>6} {r['med1_residues']:>6} "
              f"{r['hbonds']:>8} {r['salt_bridges']:>8} {r['hydrophobic']:>10}  {r['character']}")

    return results


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    # ── Single file mode (command-line) ──────────────────────
    if len(sys.argv) > 1:
        pdb_file, dlat_chains, med1_chain, med1_offset, run_name = parse_args()
        report_run(pdb_file, dlat_chains, med1_chain, med1_offset, run_name)

    # ── Batch mode (edit paths below for your local files) ───
    else:
        print("Running in batch mode — edit the pdb_files dictionary below")
        print("to point to your ColabFold output PDB files.\n")

        # ── Edit these paths to match your local directory structure ──
        pdb_files = {
            "Run 1 — DLAT monomer + MED1 TED (ipTM 0.148)": (
                "pdb/run1_dlat_monomer_med1_ted_rank001.pdb",
                ['A'], 'B', 71
            ),
            "Run 2 — DLAT monomer + MED1 LXXLL (ipTM 0.344)": (
                "pdb/run2_dlat_monomer_med1_lxxll_rank001.pdb",
                ['A'], 'B', 589
            ),
            "Run 3 — DLAT trimer + MED1 LXXLL (ipTM 0.505)": (
                "pdb/run3_dlat_trimer_med1_lxxll_rank001.pdb",
                ['A', 'B', 'C'], 'D', 589
            ),
            "Run 4 — DLAT trimer + MED1 TED (ipTM 0.412)": (
                "pdb/run4_dlat_trimer_med1_ted_rank001.pdb",
                ['A', 'B', 'C'], 'D', 71
            ),
        }

        # Filter to existing files only
        available = {k: v for k, v in pdb_files.items()
                     if os.path.exists(v[0])}

        if not available:
            print("No PDB files found at the configured paths.")
            print("Please edit the pdb_files dictionary in main() with your file paths.")
            print("Or run in single-file mode:")
            print("  python dlat_med1_interface_analysis.py <pdb_file> <dlat_chains> <med1_chain> <offset>")
        else:
            run_all_analysis(available)
