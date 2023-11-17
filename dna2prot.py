#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author  : Emmanuel Clostres
Mail    : emmanuel.clostres@univ-rennes.fr
Python  : v3.8
Date    : 18/11/2023

Description :
Translate nucleotide sequences.
Frame 1 to 3 are translated from the forward strand, and frame 4 to 6 are translated from the
reverse strand (using the reverse complement of the forward strand).

LICENSE :
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see
<https://www.gnu.org/licenses/>.
"""

import argparse
from Bio import SeqIO

VERSION = "1.0"

GENCODE = {
      'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
      'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
      'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
      'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
      'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
      'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
      'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
      'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
      'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
      'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
      'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
      'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
      'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
      'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
      'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
      'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}


def parse_arguments():
    """
    Parse the system arguments.
    """
    parser = argparse.ArgumentParser(description=f'dna2prot v{VERSION} translates a nucleotide sequence into a protein '
                                                 f'in all frames')
    parser.add_argument('-i', '--input', required=True, help='path to input nucleotide fasta file')
    parser.add_argument('-o', '--output', required=True, help='path to output amino acid fasta file')
    parser.add_argument('-3', '--forward', default=False, action='store_true', help='translate forward strand only')
    args = parser.parse_args()
    return args


def translate(sequence: str) -> str:
    """
    Translate a DNA sequence using the genetic code.
    """
    prot = []
    for n in range(len(sequence) // 3):
        prot.append(GENCODE.get(sequence[3*n:3*n+3], 'X'))
    return ''.join(prot)


def rev_comp(sequence: str) -> str:
    """
    This function returns a reverse complement of a DNA sequence.
    """
    # Complement strand
    sequence = sequence.lower()
    sequence = sequence.replace("a", "T").replace("c", "G").replace("g", "C").replace("t", "A")
    # Reverse strand
    sequence = sequence[::-1]
    return sequence.upper()


def main():
    # Import parameters
    args = parse_arguments()
    inpath = args.input
    outpath = args.output
    forward = args.forward

    print(f"inpath  : {inpath}")
    print(f"outpath : {outpath}")
    print(f"forward : {forward}")

    # Processing
    infile = SeqIO.parse(open(inpath), 'fasta')
    outfile = open(outpath, 'wt')

    for fasta in infile:
        header, seq = fasta.id, str(fasta.seq)
        seq = seq.upper()
        # print(f"header   : {header}")
        # print(f"sequence : {seq}")

        # Forward strand
        for n in range(3):
            outfile.write(f">{header}_frame{n+1}\n{translate(seq[n:])}\n")

        # Reverse strand
        if not forward:
            revcomp_seq = rev_comp(seq)
            for n in range(3):
                outfile.write(f">{header}_frame{n+4}\n{translate(revcomp_seq[n:])}\n")

    outfile.close()


if __name__ == '__main__':
    main()

# example = "ATCTGTTGAAGGATTCAGCCATGGCCGTCCTACGG"
# example_with_N = "ATCNGTTGAAGGATTCAGCCATGGCCGTCCTACGG"
# for seq in [example, example_with_N]:
#     print(f"\nseq :")
#     print(f"phase 1 (+) : {translate(seq[0:])}")
#     print(f"phase 2 (+) : {translate(seq[1:])}")
#     print(f"phase 3 (+) : {translate(seq[2:])}")
#     # print(f"phase 4 (-) : {translate(seq[::-1][0:])}")
#     # print(f"phase 5 (-) : {translate(seq[::-1][1:])}")
#     # print(f"phase 6 (-) : {translate(seq[::-1][2:])}")
