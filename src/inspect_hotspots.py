#!/usr/bin/env python3
"""
inspect_hotspots.py — analyze hotspot residues in a target PDB to decide
how to crop it for binder design.

Usage:
    python inspect_hotspots.py 5XWP_clean.pdb A 411,421,473,995

Reports, for the given chain:
  * total residue count and numbering range
  * CA-CA distances between all hotspot pairs (are they clustered?)
  * a suggested crop window that contains all hotspots plus padding
"""
import sys
from itertools import combinations

try:
    from Bio.PDB import PDBParser
except ImportError:
    sys.exit("Biopython not found. Run inside the BindCraft env.")

def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: python inspect_hotspots.py <pdb> <chain> <r1,r2,...>")
    pdb, chain_id, hotspot_str = sys.argv[1], sys.argv[2], sys.argv[3]
    hotspots = [int(x) for x in hotspot_str.split(",")]

    structure = PDBParser(QUIET=True).get_structure("t", pdb)
    model = structure[0]
    if chain_id not in model:
        sys.exit(f"Chain {chain_id} not in {pdb}. Chains: {[c.id for c in model]}")
    chain = model[chain_id]

    # collect standard residues with a CA
    resnums = []
    ca = {}
    for res in chain:
        if res.id[0] != " ":      # skip hetero/water
            continue
        if "CA" in res:
            n = res.id[1]
            resnums.append(n)
            ca[n] = res["CA"].coord

    resnums.sort()
    print(f"Chain {chain_id}: {len(resnums)} residues, "
          f"numbering {resnums[0]}..{resnums[-1]}")

    # check hotspots exist
    for h in hotspots:
        tag = "OK" if h in ca else "*** NOT FOUND (check numbering!) ***"
        print(f"  hotspot {h}: {tag}")

    present = [h for h in hotspots if h in ca]
    if len(present) < 2:
        sys.exit("Need >=2 valid hotspots to compute geometry.")

    print("\nCA-CA distances between hotspots (Angstrom):")
    for a, b in combinations(present, 2):
        import numpy as np
        d = float(np.linalg.norm(ca[a] - ca[b]))
        flag = "  <-- far apart (likely different domains)" if d > 40 else ""
        print(f"  {a} - {b}: {d:6.1f}{flag}")

    pad = 20
    lo, hi = min(present) - pad, max(present) + pad
    lo = max(lo, resnums[0])
    hi = min(hi, resnums[-1])
    print(f"\nSuggested crop window (hotspot span +/- {pad}): {lo}..{hi}")
    print(f"  -> {hi - lo + 1} residues kept (was {len(resnums)})")
    print("\nNOTE: if hotspots are far apart, a single contiguous crop may be")
    print("      too large or may not make biological sense. In that case,")
    print("      decide which hotspots are the true interface and crop to those.")

if __name__ == "__main__":
    main()
