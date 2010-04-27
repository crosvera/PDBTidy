#!/usr/bin/env python
# coding: UTF-8
# Copyright (C) 2010, Carlos Ríos V. <crosvera@gmail.com>
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.  

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.PDB.Polypeptide import to_one_letter_code
from Bio import PDB

def get_sequence(structure, format, chain=None):
    """ Return a sequence from a structure with the given format.

    Arguments:
    o structure - Structure object, the protein structurewhere the sequence
                  will be extracted.
    o format - string, the format which the sequence will be returned.
    o chain - string, optional, the chain that will be extracted, if None
              will extract the sequences from all the chains.
    """

    records = []
    for model in structure:
        for c in model:
            if chain == c.id or chain == None:
                s = "".join( to_one_letter_code.get(res.resname,"X") \
                          for res in c if "CA" in res.child_dict )
                record = SeqRecord( Seq(s, IUPAC.protein),
                        id = structure.id+":"+c.id,
                        name = structure.header['head'],
                        description = structure.header['name'] )

                records.append( record )

    seqs = ""
    for r in records:
        seqs += r.format(format)
        seqs += "\n"

    return seqs



if __name__ == "__main__":
    import sys, os

    if len(sys.argv) > 4 or len(sys.argv) < 3:
        print "Invalid Arguments!!!!"
        sys.exit(-1)


    # argv[1] -> PDB file
    # argv[2] -> Format
    # argv[3] -> Chain

    # Example:
    #   % python get_sequence.py 13GS.pdb fasta
    # Or:
    #   % python get_sequence.py 13GS.pdb fasta B

    parser = PDB.PDBParser()

    n = sys.argv[1].split(os.sep)[-1].strip(".pdb")
    structure = parser.get_structure(n, sys.argv[1])

    if len(sys.argv) == 3:
        print get_sequence(structure, sys.argv[2])
    else:
        print get_sequence(structure, sys.argv[2], sys.argv[3])
