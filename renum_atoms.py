#!/usr/bin/env python
# coding: UTF-8
# Copyright (C) 2011, Carlos Ríos V. <crosvera@gmail.com>
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import PDB

def renum_atoms(structure, output, model_id=0, start=1):
    """
    renum_atoms will renumber the atoms starting with start...
    """

    model = structure[model_id]
    atoms = list(model.get_atoms())
    for a in atoms:
        a.set_serial_number(start)
        start += 1

    w = PDB.PDBIO()
    w.set_structure(structure)
    w.save(output, conserve_atoms_number=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('pdb_in', 
        help='PDB file to be renumbered')
    parser.add_argument('pdb_out', 
        help='PDB file where the renumbered atoms will be stored.')
    parser.add_argument('--start', dest='start',
        type=int, default=1, help='Starting number for renumber.')
    parser.add_argument('--model', dest='model',
        type=int, default=0, help='Model (id) where the renumber will be applied.')

    opts =  parser.parse_args()

    structure = PDB.PDBParser().get_structure(opts.pdb_in, opts.pdb_in)
    renum_atoms(structure, opts.pdb_out, opts.model, opts.start)
