"""Tester"""
import CrisprACT3Finder

input_path = "/groups/itay_mayrose/josefbrook/projects/Gene_Activation_Library/in"
output_path = "/groups/itay_mayrose/josefbrook/projects/Gene_Activation_Library/out"
in_fasta_file = "/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
in_gff_file = "/Arabidopsis_thaliana.TAIR10.54.gff3"
in_genes_list = ["AT1G01010", "AT1G01110", "AT1G01020"]
in_pams = ['AGG', 'GGG', 'CGG', 'TGG']
results = (CrisprACT3Finder.sgrnas_for_genes(output_path, input_path, in_genes_list, input_path+in_fasta_file, input_path+in_gff_file,
                                             "ucrispr", in_pams, 200, 0, 45, 60))
for cand in results:
    print(str(cand))
