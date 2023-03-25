"""Tester"""
# import time
# from typing import List

# import subprocess

import os

import CrisprACT3Finder

PATH = os.path.dirname(os.path.realpath(__file__))

input_path = "/groups/itay_mayrose/josefbrook/projects/Gene_Activation_Library/in"
output_path = "/groups/itay_mayrose/josefbrook/projects/Gene_Activation_Library/test"
in_fasta_file = "/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
in_gff_file = "/Arabidopsis_thaliana.TAIR10.54.gff3"
gene_ids_file = "/groups/itay_mayrose/josefbrook/projects/Gene_Activation_Library/in/gene_ids_test.txt"
in_pams = ('AGG', 'GGG', 'CGG', 'TGG')
results = (CrisprACT3Finder.sgrnas_for_genes(output_path, input_path, gene_ids_file, input_path+in_fasta_file, input_path+in_gff_file,
                                             "ucrispr", in_pams, 300, 0, 45, 60, 1, 3))
for cand in results:
    print(str(cand))
