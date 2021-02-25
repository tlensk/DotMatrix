# DotMatrix
### Tatiana Lenskaia and Daniel Boley

This repository contains code and examples for the dot-matrix display implemented using dictionary-based methods.

## Dot.py
The dot-matrix display script in Python 3 to visualize similarity between sequences.
The script visualizes strings of length _m_ shared between sequences.

To run the script, one can run the following command in the directory that contains input files for sequences:

`python3 Dot.py m file1.fasta file2.fasta`

to compare two individual sequences

or

`python3 Dot.py m multifastafile`

for pairwise comparison between the sequences in multifasta file



Color scheme:
* blue - exact match of length _m_ on a forward strand
* black - the end of continious match on a forward strand
* red - exact match of length _m_ on a reverse strand
* purple - the end of continious match on a reverse strand
* green - SNP on either forward or reverse strand


## Example 1
This folder contains figures generated by Dot.py for a bacterium-varus pair example.


## Example 2
This folder contains pictures generated by Dot.py for Coronavirus example.

