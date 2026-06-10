#!/usr/bin/env python3
"""
crop_target.py — crop a target PDB to a residue window for binder design.

Keeps a contiguous residue range from one chain, drops everything else
(other chains, hetero/water), and writes the result. By default it KEEPS
the original residue numbering so your hotspot IDs (e.g. A411) stay valid.

Usage:
    python crop_target.py IN.pdb OUT.pdb A 391 1015
    #                       in     out   chain lo hi

    # optional: renumber the kept residues to start at 1
    python crop_target.py IN.pdb OUT.pdb A 391 1015 --renumber

If you renumber, remember to update the hotspot residue IDs in your
settings JSON to the new numbers.
"""
import sys
from Bio.PDB import PDBParser, PDBIO, Select

def main():
    args = sys.argv[1:]
    renumber = "--renumber" in args
    args = [a for a in args if a != "--renumber"]
    if len(args) != 5:
        sys.exit("Usage: python crop_target.py IN.pdb OUT.pdb CHAIN LO HI [--renumber]")
    inp, outp, chain_id = args[0], args[1], args[2]
    lo, hi = int(args[3]), int(args[4])

    structure = PDBParser(QUIET=True).get_structure("t", inp)
    model = structure[0]
    if chain_id not in model:
        sys.exit(f"Chain {chain_id} not found. Present: {[c.id for c in model]}")

    class Crop(Select):
        def accept_chain(self, chain):
            return chain.id == chain_id
        def accept_residue(self, res):
            return res.id[0] == " " and lo <= res.id[1] <= hi
        def accept_atom(self, atom):
            return not atom.is_disordered() or atom.get_altloc() in ("", "A")

    if renumber:
        chain = model[chain_id]
        kept = [r for r in chain if r.id[0] == " " and lo <= r.id[1] <= hi]
        # renumber sequentially from 1; offset shown for hotspot remap
        offset = lo - 1
        for r in kept:
            het, num, ins = r.id
            r.id = (het, num - offset, ins)
        print(f"Renumbered: old {lo}..{hi} -> new 1..{hi-lo+1}  (subtract {offset} from hotspots)")

    io = PDBIO()
    io.set_structure(structure)
    io.save(outp, Crop() if not renumber else None)
    # when renumbering we already mutated ids; still restrict chain/range:
    if renumber:
        io.save(outp, Crop())

    print(f"Wrote {outp}: chain {chain_id}, residues {lo}..{hi} "
          f"({hi-lo+1} residues)")

if __name__ == "__main__":
    main()
