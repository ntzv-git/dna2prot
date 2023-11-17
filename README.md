# dna2prot

<pre>
Author  : Emmanuel Clostres
Mail    : emmanuel.clostres@univ-rennes.fr
Python  : v3.8+
Version : 1.0
</pre>


## Description

Translate nucleotide sequences. Frame 1 to 3 are translated from the forward strand, and frame 4 to 6 are translated from the reverse strand (using the reverse complement of the forward strand).


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
