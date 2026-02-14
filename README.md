# dna2prot

<pre>
Author  : Emmanuel Clostres
Mail    : emmanuel.clostres@univ-rennes.fr
Python  : v3.8+
Version : 1.0
</pre>


## Description

dna2prot is a lightweight Python tool designed to translate nucleotide sequences across all 3 or 6 possible reading frames.

Unlike tools like EMBOSS `getorf` that selectively extract Open Reading Frames (ORFs) bounded by START and STOP codons, this script performs a raw, uninterrupted translation of the entire input sequence from end to end. Whenever it encounters a STOP codon, it simply inserts an asterisk (*) and continues translating the remaining nucleotides without breaking the sequence. This brute-force approach guarantees exactly 3 or 6 continuous protein strings per input, making it ideal for researchers who need to visualize the complete translational landscape of a sequence without any algorithmic filtering or size thresholds.


## Installation

```bash
git clone https://github.com/ntzv-git/dna2prot.git
```


### Requirements

- python3.8


### Package dependencies

- argparse
- Bio


## Usage

```
usage: dna2prot.py [-h] -i INPUT -o OUTPUT [-3]

dna2prot v1.0 translates a nucleotide sequence into a protein in all frames

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to input nucleotide fasta file
  -o OUTPUT, --output OUTPUT
                        path to output amino acid fasta file
  -3, --forward         translate forward strand only
```
