## Description
This project involves 2 main purposes:
1. To get the DNA sequence of MC1R gene by querying mygene and ensembl databases, find the longest open reading frame (ORF), and then convert it to amino acid sequence.
2. Find out what other species have genes that are homologous to MC1R.

## Input
There are no specific input files. We are querying data through API.

## Output
1. `mc1r.fasta` -- a fasta file that contains full DNA sequence of MC1R gene, DNA sequence of the longest ORF, and amino acid sequence of the longest ORF.
2. `mc1r_homology_list.txt` -- a txt file with the full list of unique species that have genes homologous to MC1R.

## Usage
This script can be called with the following command: `python3 assignment4_2.7.py`.