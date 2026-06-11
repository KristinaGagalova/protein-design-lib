#!/usr/bin/env python3
"""
clean_pdb.py -- clean a PDB for binder design / structure prediction.

Default behaviour:
  * keep all standard protein ATOM records
  * remove waters (HOH) and other HETATM records
  * keep TER records (chain terminators) where chains actually end
  * drop hydrogens (optional), altlocs (keep 'A'/' '), and insertion handling

Examples:
  # strip all HETATM + waters, keep everything else:
  python clean_pdb.py in.pdb 5XWP_clean.pdb

  # keep only chain A:
  python clean_pdb.py in.pdb 5XWP_clean.pdb --chains A

  # strip waters but KEEP other heteroatoms (e.g. a ligand/metal you need):
  python clean_pdb.py in.pdb out.pdb --keep-hetatm

  # also keep specific HETATM residue names even in default strip mode:
  python clean_pdb.py in.pdb out.pdb --keep-resn ZN,HEM
"""
import argparse
import sys

WATER = {"HOH", "WAT", "DOD", "TIP", "TIP3", "TIP4", "SOL"}

# Standard 20 amino acids (three-letter). RFdiffusion only understands these.
STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

# Common modified residues -> standard equivalent. These appear as HETATM or
# non-standard ATOM names in crystal structures; converting (not deleting)
# avoids leaving chain gaps that break RFdiffusion.
MODIFIED_TO_STANDARD = {
    "MSE": "MET",   # selenomethionine -> methionine (very common)
    "SEC": "CYS",   # selenocysteine
    "PYL": "LYS",   # pyrrolysine
    "HYP": "PRO",   # hydroxyproline
    "PTR": "TYR",   # phosphotyrosine
    "SEP": "SER",   # phosphoserine
    "TPO": "THR",   # phosphothreonine
    "MLY": "LYS",   # methyllysine
    "CSO": "CYS",   # oxidized cysteine
}

def main():
    ap = argparse.ArgumentParser(description="Clean a PDB file.")
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--chains", default=None,
                    help="Comma-separated chains to KEEP (e.g. A or A,B). "
                         "Default: keep all chains.")
    ap.add_argument("--keep-hetatm", action="store_true",
                    help="Keep all HETATM records (only waters removed).")
    ap.add_argument("--keep-resn", default="",
                    help="Comma-separated HETATM resnames to keep even when "
                         "stripping HETATM (e.g. ZN,HEM,NAD).")
    ap.add_argument("--keep-water", action="store_true",
                    help="Keep water molecules (removed by default).")
    ap.add_argument("--keep-hydrogens", action="store_true",
                    help="Keep hydrogen atoms (removed by default).")
    ap.add_argument("--keep-altloc", action="store_true",
                    help="Keep all altloc conformers (default keeps only "
                         "the first: altloc ' ' or 'A').")
    ap.add_argument("--rfdiffusion", action="store_true",
                    help="STRICT preset for RFdiffusion/RFD3: keep ONLY standard "
                         "protein residues (backbone). Removes ALL heteroatoms, "
                         "waters, hydrogens, and extra altlocs; converts common "
                         "modified residues (MSE->MET etc.) to standard; drops "
                         "any residue not in the standard 20. Overrides keep-* "
                         "flags for heteroatoms.")
    args = ap.parse_args()

    keep_chains = set(args.chains.split(",")) if args.chains else None
    keep_resn = set(r.strip().upper() for r in args.keep_resn.split(",") if r.strip())

    # In rfdiffusion mode, force strict heteroatom removal regardless of other flags.
    if args.rfdiffusion:
        keep_resn = set()
        args.keep_hetatm = False
        args.keep_water = False
        args.keep_hydrogens = False
        # altlocs already collapsed to first by default

    kept_atom = kept_het = 0
    dropped_water = dropped_het = dropped_h = dropped_alt = dropped_chain = 0
    converted = dropped_nonstd = 0

    out_lines = []
    with open(args.infile) as fh:
        for line in fh:
            rec = line[:6].strip()

            if rec == "TER":
                out_lines.append(line)
                continue

            if rec not in ("ATOM", "HETATM"):
                out_lines.append(line)
                continue

            # --- ATOM/HETATM coordinate line ---
            altloc = line[16]
            resn = line[17:20].strip().upper()
            chain = line[21]
            element = line[76:78].strip()
            atom_name = line[12:16].strip()

            # chain filter
            if keep_chains is not None and chain not in keep_chains:
                dropped_chain += 1
                continue

            # altloc filter
            if not args.keep_altloc and altloc not in (" ", "A"):
                dropped_alt += 1
                continue

            # hydrogen filter
            is_h = element == "H" or (not element and atom_name.startswith(("H", "1H", "2H", "3H")))
            if not args.keep_hydrogens and is_h:
                dropped_h += 1
                continue

            # -------- RFdiffusion strict mode --------
            if args.rfdiffusion:
                target_resn = resn
                # convert known modified residues to their standard equivalent
                if resn in MODIFIED_TO_STANDARD:
                    target_resn = MODIFIED_TO_STANDARD[resn]
                    # rewrite resname field (cols 18-20) and turn HETATM into ATOM
                    line = "ATOM  " + line[6:17] + f"{target_resn:>3}" + line[20:]
                    converted += 1
                # after conversion, must be a standard amino acid; else drop
                if target_resn not in STANDARD_AA:
                    dropped_nonstd += 1
                    if rec == "HETATM":
                        dropped_het += 1
                    continue
                kept_atom += 1
                out_lines.append(line)
                continue

            # -------- normal (non-rfdiffusion) mode --------
            if rec == "HETATM":
                if resn in WATER and not args.keep_water:
                    dropped_water += 1
                    continue
                if not args.keep_hetatm and resn not in keep_resn:
                    dropped_het += 1
                    continue
                kept_het += 1
            else:
                kept_atom += 1

            out_lines.append(line)

    with open(args.outfile, "w") as fh:
        fh.writelines(out_lines)

    print(f"Wrote {args.outfile}", file=sys.stderr)
    if args.rfdiffusion:
        print(f"  RFdiffusion strict mode", file=sys.stderr)
        print(f"  kept   : {kept_atom} standard-protein atoms", file=sys.stderr)
        print(f"  converted modified residues: {converted} (e.g. MSE->MET)", file=sys.stderr)
        print(f"  dropped: {dropped_water} water, {dropped_het} HETATM, "
              f"{dropped_nonstd} non-standard-residue atoms, "
              f"{dropped_h} hydrogens, {dropped_alt} altlocs, "
              f"{dropped_chain} unwanted-chain atoms", file=sys.stderr)
    else:
        print(f"  kept   : {kept_atom} ATOM, {kept_het} HETATM", file=sys.stderr)
        print(f"  dropped: {dropped_water} water, {dropped_het} other HETATM, "
              f"{dropped_h} hydrogens, {dropped_alt} altlocs, "
              f"{dropped_chain} atoms in unwanted chains", file=sys.stderr)

if __name__ == "__main__":
    main()
