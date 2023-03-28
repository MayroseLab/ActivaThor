"""Find sgRNAs for CRISPR-Act3"""
import argparse
import subprocess
import time
import warnings
from itertools import groupby

from typing import Tuple
import gffpandas.gffpandas as gffpd
from operator import attrgetter

from DeepHF.LoadDeepHF import load_deephf
from GetCasSites import get_sites
import ScoringFunctions
from FindOffTargets import *

warnings.filterwarnings('ignore')


def choose_scoring_function(input_scoring_function: str):
    """
    This function translates the chosen input for the scoring function and returns its pointer in the algorithm files
    and whether the function takes PAMs into its calculation.

    :param input_scoring_function: chosen scoring function by the user
    :return: scoring function pointer to use in the algorithm
    :rtype: function, boolean
    """
    if input_scoring_function == "ucrispr":
        return ScoringFunctions.ucrispr
    elif input_scoring_function == "deephf":
        return load_deephf()
    else:
        print("Did not specify valid scoring function")


def gene_ids_to_list(file_path: str) -> List[str]:
    """
    Read the genes IDs file and save the IDs into a list.

    :param file_path: full path to the gene IDs txt file
    :return: a list of gene IDs
    """
    # opening the file in read mode
    genes_file = open(f"{file_path}", "r")
    # reading the file
    data = genes_file.read().upper()
    # replacing end splitting the text when newline ('\n') is seen.
    gene_ids_into_list = data.split("\n")
    genes_file.close()
    return gene_ids_into_list


def add_promoter_region_indices(genes_filt_df, upstream: int, downstream: int):
    """

    :param genes_filt_df:
    :param upstream:
    :param downstream:
    """
    # create promoter site start indices
    genes_filt_df['site_start'] = genes_filt_df.apply(
        lambda x: x['start'] - upstream - 1 if x['strand'] == '+' else x['end'] + downstream,
        axis=1)
    # create promoter site end indices
    genes_filt_df['site_end'] = genes_filt_df.apply(
        lambda x: x['start'] + downstream - 1 if x['strand'] == '+' else x['end'] + upstream,
        axis=1)


def filter_gff_file(gff_file: str, genes_ids: List[str], upstream: int, downstream: int, out_path: str):
    """
    Parse the input GFF file with gffpandas, filter only the annotations for genes and add columns of promoter site
    start and site end.

    :param gff_file: path to input GFF format file of the genes annotations
    :param genes_ids: list of input gene IDs
    :param upstream: number of base-pairs upstream to the TSS to search potential sgRNA target sites
    :param downstream: number of base-pairs downstream to the TSS to search potential sgRNA target sites
    :param out_path:
    :return: gffpandas data frame
    """
    # open annotations file with gff pandas
    annotations = gffpd.read_gff3(gff_file)
    # filter relevant genes. Inform the user if an input gene ID was not found in the GFF file
    filtered_annotations = annotations.filter_feature_of_type(['gene'])
    # genes_filt_df = filtered_annotations.get_feature_by_attribute('gene_id', genes_ids).df
    genes_filt_gff = filtered_annotations.get_feature_by_attribute('gene_id', genes_ids)
    filt_annotations_attr = filtered_annotations.attributes_to_columns()
    for gene_id in genes_ids:
        if gene_id not in list(filt_annotations_attr['gene_id']):
            print(f"No matching gene ID={gene_id} was found in the GFF file")
    # add promoter site start and site end indices to the data frame
    genes_filt_df = genes_filt_gff.attributes_to_columns()
    add_promoter_region_indices(genes_filt_df, upstream, downstream)
    genes_filt_df['start'] = genes_filt_df.apply(lambda x: x['site_start'], axis=1)
    genes_filt_df['end'] = genes_filt_df.apply(lambda x: x['site_end'], axis=1)
    # change the type of the added sites to "promoter"
    genes_filt_df.loc[genes_filt_df["type"] == "gene", "type"] = "promoter"
    # append the promoter annotations to the GFF
    new_anno_df = filt_annotations_attr.append(genes_filt_df, ignore_index=True)
    new_anno_df['attributes'] = new_anno_df.apply(lambda x: x['gene_id'], axis=1)
    new_anno_df = new_anno_df.drop(new_anno_df.iloc[:, 9:], axis=1)
    annotations.df = new_anno_df
    gff_with_proms_path = "/gff_with_proms.gff3"
    annotations.to_gff3(out_path + gff_with_proms_path)
    # annotations.df = genes_filt_df
    # promoters_df = annotations.attributes_to_columns()
    return genes_filt_df, gff_with_proms_path


def get_upstream_sites(out_path: str, fasta_file: str, filtered_gene_df) -> List[str]:
    """Find the sequences upstream to the gene's TSS

    :param out_path: the directory path to which the algorithm will store the results
    :param fasta_file: path to input FASTA format file of the genome
    :param filtered_gene_df:
    :return:
    """
    # create BED format file
    filtered_gene_df.to_csv(out_path + '/genes.bed', sep='\t',
                            columns=['seq_id', 'start', 'end', 'gene_id', 'score', 'strand'],
                            header=False, index=False)
    bed_file = out_path + "/genes.bed"  # TODO validate PATHs
    # run bedtools
    seq = subprocess.run(['bedtools', 'getfasta', '-fi', fasta_file, '-bed', bed_file, '-nameOnly', '-s'],
                         stdout=subprocess.PIPE)
    sites_list = seq.stdout.decode().split()
    return sites_list


def get_targets_from_site(upstream_site: str, pams: Tuple) -> Tuple[List[str], List[str], List[str], List[str]]:
    """
    Find potential CRISPR target sites in a given DNA sequence, and return a list of the target sequences.

    :param upstream_site: DNA sequence of a gene TSS upstream site
    :param pams: list of PAMs by which potential sgRNA target sites will be found
    :return: a list of targets sequences that CRISPR can target
    """
    found_fwd_targets_first, found_fwd_targets_second, found_rev_targets_first, found_rev_targets_second = get_sites(upstream_site, pams)
    return found_fwd_targets_first, found_fwd_targets_second, found_rev_targets_first, found_rev_targets_second


def scores_batch(sgrnas: List[str], scoring_function, out_path: str) -> List[float]:
    """
    Calculate on target scores for sgRNAs in 'sgrnas' using the given scoring function.

    :param out_path: the directory path to which the algorithm will store the results
    :param sgrnas: list of candidate sgRNAs
    :param scoring_function: chosen on-target scoring function
    :return:
    """
    scores = scoring_function(sgrnas, out_path)
    return scores


def create_candidates(promoter_range_rank: int, targets_lst: List[str], gene_id: str,
                      chromosome: str, act_site_start: int, act_site_end: int, gene_strand: str, targets_strand: str,
                      min_gc_cont: int, max_gc_cont: int) -> List[ActivationCandidate]:
    """

    :param promoter_range_rank: defines the range where the target was found in the promoter. 1 for 0-200, 2 for 201-300.
    :param targets_lst: list of target sequences found in a specific gene promoter

    :param gene_id: gene ID from the gff annotations file
    :param chromosome: number from the gff annotations file
    :param act_site_start: promoter region start location in the chromosome
    :param act_site_end: promoter region end location in the chromosome
    :param gene_strand: gene strand ( + or - ) from the gff annotations file
    :param targets_strand: strand relatively to the gene strand, on which the targets were found
    :param min_gc_cont: percentage of minimum GC content by which to filter sgRNA candidates
    :param max_gc_cont: percentage of maximum GC content by which to filter sgRNA candidates
    :return: list of sgRNAs as ActivationCandidate objects of the specific gene
    """

    candidates_list = []
    for i in range(len(targets_lst)):
        target_strand = "+"
        if (gene_strand == "-" and targets_strand == "+") or (gene_strand == "+" and targets_strand == "-"):
            target_strand = "-"
        candidate = ActivationCandidate(targets_lst[i][:20], targets_lst[i][20:], gene_id, chromosome, act_site_start,
                                        act_site_end, target_strand, promoter_range_rank, 0, 0, {}, 0)
        candidate.calc_GC_content(min_gc_cont, max_gc_cont)
        candidate.calc_seq_repetitions()

        candidates_list.append(candidate)
    return candidates_list


def return_candidates(out_path: str, sites_lst: List[str], genes_filt_df, scoring_function, pams: Tuple, min_gc_cont: int,
                      max_gc_cont: int, top_num: int = 3) -> List[ActivationCandidate]:
    """

    :param out_path: the directory path to which the algorithm will store the results
    :param sites_lst: list of short genes annotations and promoter sites
    :param genes_filt_df:
    :param scoring_function: chosen on-target scoring function
    :param pams: list of PAMs by which potential sgRNA target sites will be found
    :param min_gc_cont: percentage of minimum GC content by which to filter sgRNA candidates
    :param max_gc_cont: percentage of maximum GC content by which to filter sgRNA candidates
    :param top_num: number of top scored candidates per gene
    :return: list of sgRNAs as ActivationCandidate objects
    """
    candidates_list = []
    for i in range(0, len(sites_lst), 2):
        gene_candidates = []
        # create gene annotations
        gene_id = sites_lst[i][1:-3]  # TODO find with regex or any other tool of validation
        chromosome = genes_filt_df[genes_filt_df['gene_id'] == gene_id]['seq_id'].values[0]
        act_site_start = int(genes_filt_df[genes_filt_df['gene_id'] == gene_id]['start'].values[0])
        act_site_end = int(genes_filt_df[genes_filt_df['gene_id'] == gene_id]['end'].values[0])
        strand = genes_filt_df[genes_filt_df['gene_id'] == gene_id]['strand'].values[0]
        # extract target sequences from each TSS upstream site
        targets_fwd_1, targets_fwd_2, targets_rev_1, targets_rev_2 = get_targets_from_site(sites_lst[i + 1].upper(), pams)
        # loop through targets in forward strand
        gene_candidates += create_candidates(1, targets_fwd_1, gene_id, chromosome, act_site_start, act_site_end,
                                             strand, "+", min_gc_cont, max_gc_cont)
        gene_candidates += create_candidates(2, targets_fwd_2, gene_id, chromosome, act_site_start, act_site_end,
                                             strand, "+", min_gc_cont, max_gc_cont)
        gene_candidates += create_candidates(1, targets_rev_1, gene_id, chromosome, act_site_start, act_site_end,
                                             strand, "-", min_gc_cont, max_gc_cont)
        gene_candidates += create_candidates(2, targets_rev_2, gene_id, chromosome, act_site_start, act_site_end,
                                             strand, "-", min_gc_cont, max_gc_cont)
        candidates_list += gene_candidates
    # calculate on-target scores for the candidate sgRNAs
    t0 = time.perf_counter()
    scores = scores_batch([candidate.seq+candidate.pam for candidate in candidates_list], scoring_function, out_path)
    t1 = time.perf_counter()
    print(f"scoring function for {len(candidates_list)} candidates ran in {t1 - t0} seconds")
    for i in range(len(candidates_list)):
        candidates_list[i].on_score = scores[i]
    # sort candidates list for each gene ID by GC content, nucleotide repetitions and on-target score
    gene_id_attr_sort = sorted(candidates_list, key=attrgetter('gene'))
    # Group the objects by 'gene_id'
    grouped_candidates = {gene_id: list(objects) for gene_id, objects in groupby(gene_id_attr_sort, key=attrgetter('gene'))}
    # Create a dictionary with gene_id as key and value as a list of top 3 ranked objects
    candidates_list_sorted = []
    for gene_id, objects in grouped_candidates.items():
        # candidates_list_sorted += sorted(objects, key=lambda c: (c.gene, c.promoter_range_rank,
        #                             c.gc_content_category, c.nuc_rep_score, -c.on_score))[:top_num]
        candidates_list_sorted += sorted(objects, key=lambda c: (c.gene, c.promoter_range_rank, -c.on_score))[:top_num]
    return candidates_list_sorted


def results_to_dataframe(candidates_list: List[ActivationCandidate], out_path: str, with_off_targets: int = 1) -> pd.DataFrame:
    """

    :param candidates_list:
    :param out_path: the directory path to which the algorithm will store the results
    :param with_off_targets: choose whether to search for off-targets in the algorithm run
    :return:
    """
    df = pd.DataFrame([candidate.to_dict() for candidate in candidates_list]).reset_index()
    df['Rank'] = df.groupby('gene')['index'].transform('rank')
    if with_off_targets == 1:
        header = ['Rank', 'seq', "pam", "gene", "chromosome", "act_site_start", "act_site_end", "strand",
              "promoter_range_rank", "gc_content", "gc_content_category", 'nucleotide_repetitions', 'nuc_rep_score', "on_score", "off1", "off2"]
    else:
        header = ['Rank', 'seq', "pam", "gene", "chromosome", "act_site_start", "act_site_end", "strand",
                  "promoter_range_rank", "gc_content", "gc_content_category", 'nucleotide_repetitions', 'nuc_rep_score', "on_score"]
    df.to_csv(out_path + "/results.csv", sep=',', columns=header, index=False)
    return df


def sgrnas_for_genes(out_path: str, in_path: str, genes_ids_file: str, fasta_file: str, gff_file: str,
                     scoring_function: str = 'ucrispr',
                     pams: Tuple = ('AGG', 'GGG', 'CGG', 'TGG'), bp_upstream: int = 300, bp_downstream: int = 0, min_gc_cont: int = 45,
                     max_gc_cont: int = 60, with_off_targets: int = 1, top_num: int = 3) -> List[ActivationCandidate]:
    """
    Given a list of genes ensembl IDs, a FASTA file of full or partial genome sequence, a GFF file of annotations
    of the genome and a chosen scoring function - search and return sgRNA to target the wanted genes, with their sequence
    and the efficiency score as a list of ActivationCandidate objects.

    :param out_path: the directory path to which the algorithm will store the results
    :param in_path: the directory path from which the input files will be found
    :param genes_ids_file: path to input txt format file of the gene IDs
    :param fasta_file: path to input FASTA format file of the genome
    :param gff_file: path to input GFF format file of the genes annotations
    :param scoring_function: chosen on-target scoring function
    :param pams: list of PAMs by which potential sgRNA target sites will be found
    :param bp_upstream: number of base-pairs upstream to the TSS to search potential sgRNA target sites
    :param bp_downstream: number of base-pairs downstream to the TSS to search potential sgRNA target sites
    :param min_gc_cont: percentage of minimum GC content by which to filter sgRNA candidates
    :param max_gc_cont: percentage of maximum GC content by which to filter sgRNA candidates
    :param with_off_targets: choose whether to search for off-targets in the algorithm run
    :param top_num: number of top scored candidates per gene
    :return: list of potential sgRNA as ActivationCandidate objects to target the input genes
    """
    # choosing the scoring function
    scoring_function = choose_scoring_function(scoring_function)
    # parse gene IDs text file and save into a list of gene IDs
    genes_ids = gene_ids_to_list(genes_ids_file)
    # handling the GFF annotations file - filter only the annotations for genes and add columns of promoter site start and site end.
    promoters_df, gff_with_proms_path = filter_gff_file(gff_file, genes_ids, bp_upstream, bp_downstream, out_path)
    # handling the genome FASTA file - extract promoter sites using bedtools getfasta
    sites = get_upstream_sites(out_path, fasta_file, promoters_df)
    # loop through the genes in the target sites list and create sgRNA candidates
    candidates_list = return_candidates(out_path, sites, promoters_df, scoring_function, pams, min_gc_cont, max_gc_cont, top_num)
    # get off-targets for the sgRNA candidates
    if with_off_targets == 1:
        get_off_targets(candidates_list, in_path, out_path, gff_with_proms_path)
    # results to DataFrame: add Rank column ranking the candidates of each gene, and save the results to CSV
    df = results_to_dataframe(candidates_list, out_path, with_off_targets)
    print(df)
    return candidates_list


def parse_arguments(parser_obj: argparse.ArgumentParser):
    """
    using a pars_obj object this function parses command line strings into python objects. the chosen parameters for the
    algorithm run are taken from this function.

    :param parser_obj: object for parsing command line strings into python objects
    :return: parsed parameters for the algorithm run
    :rtype: argparse.Namespace
    """
    parser_obj.add_argument('out_path', type=str, metavar='<out_path>',
                            help='The output path to the directory in which the output files will be written')
    parser_obj.add_argument('in_path', type=str, metavar='<out_path>',
                            help='The directory path from which the input files will be found')
    parser_obj.add_argument('genes', type=str, default='', help='')
    parser_obj.add_argument('fasta_file', type=str, metavar='<fasta_file>', help='The path to the input FASTA file')
    parser_obj.add_argument('gff_file', type=str, metavar='<gff_file>', help='The path to the input GFF file')
    parser_obj.add_argument('--scoring_function', '-n', type=str, default='ucrispr',
                            help='the on scoring function of the targets. Optional scoring systems are: deephf, ucrispr. '
                                 'Additional scoring function may be added by the user or by request.')
    parser_obj.add_argument('--pams', type=List, default=['AGG', 'GGG', 'CGG', 'TGG'],
                            help='list of PAMs by which potential sgRNA target sites will be found')
    parser_obj.add_argument('--bp_upstream', '-u', type=int, default=300,
                            help='number of base-pairs upstream to the TSS to search potential sgRNA target sites')
    parser_obj.add_argument('--bp_downstream', '-d', type=int, default=0,
                            help='number of base-pairs downstream to the TSS to search potential sgRNA target sites')
    parser_obj.add_argument('--min_gc_cont', '-i', type=int, default=45,
                            help='percentage of minimum GC content by which to filter sgRNA candidates')
    parser_obj.add_argument('--max_gc_cont', '-x', type=int, default=60,
                            help='percentage of maximum GC content by which to filter sgRNA candidates')
    parser_obj.add_argument('--with_off_targets', '-o', type=int, default=1,
                            help='choose whether to search for off-targets in the algorithm run')
    parser_obj.add_argument('--top_num', '-t', type=int, default=1,
                            help='number of top scored candidates per gene')
    arguments = parser_obj.parse_args()
    return arguments


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parse_arguments(parser)
    sgrnas_for_genes(out_path=args.out_path,
                     in_path=args.in_path,
                     genes_ids_file=args.genes,
                     fasta_file=args.fasta_file,
                     gff_file=args.gff_file,
                     scoring_function=args.scoring_function,
                     pams=args.pams,
                     bp_upstream=args.bp_upstream,
                     bp_downstream=args.bp_downstream,
                     min_gc_cont=args.min_gc_cont,
                     max_gc_cont=args.max_gc_cont,
                     with_off_targets=args.with_off_targets,
                     top_num=args.top_num
                     )

# crispritz chromosomes directory should contain only fasta format files
