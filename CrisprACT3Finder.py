"""Find sgRNAs for CRISPR-Act3"""
import argparse
import subprocess
import warnings

import gffpandas.gffpandas as gffpd

from DeepHF.LoadDeepHF import load_deephf
from GetCasSites import get_sites
import ScoringFunctions
from ActivationCandidate import ActivationCandidate

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


def filter_gff_file(gff_file: str, genes_ids: list[str], upstream: int, downstream: int):
    """

    :param gff_file: path to input GFF format file of the genes annotations
    :param genes_ids:
    :param upstream: number of base-pairs upstream to the TSS to search potential sgRNA target sites
    :param downstream: number of base-pairs downstream to the TSS to search potential sgRNA target sites
    :return:
    """
    # open annotations file with gff pandas
    annotation = gffpd.read_gff3(gff_file)
    # filter the dataframe by feature 'gene' - leave only annotations for genes
    filter_genes_df = annotation.filter_feature_of_type(['gene'])
    # save the attributes as columns
    attr_to_columns = filter_genes_df.attributes_to_columns()
    # filter relevant genes
    genes_filt_df = attr_to_columns[attr_to_columns['gene_id'].isin(genes_ids)]  # TODO check validity
    # create promoter site start indices
    genes_filt_df['site_start'] = genes_filt_df.apply(
        lambda x: x['start'] - upstream - 1 if x['strand'] == '+' else x['end'] + downstream,
        axis=1)
    # create promoter site end indices
    genes_filt_df['site_end'] = genes_filt_df.apply(
        lambda x: x['start'] + downstream - 1 if x['strand'] == '+' else x['end'] + upstream,
        axis=1)
    return genes_filt_df


def get_upstream_sites(out_path: str, fasta_file: str, filtered_gene_df) -> list[str]:
    """Find the sequences upstream to the gene's TSS

    :param out_path: the directory path to which the algorithm will store the results
    :param fasta_file: path to input FASTA format file of the genome
    :param filtered_gene_df:
    :return:
    """
    # create BED format file
    filtered_gene_df.to_csv(out_path + '/genes.bed', sep='\t',
                            columns=['seq_id', 'site_start', 'site_end', 'gene_id', 'score', 'strand'],
                            header=False, index=False)
    bed_file = out_path + "/genes.bed"  # TODO validate PATHs
    # run bedtools
    seq = subprocess.run(['bedtools', 'getfasta', '-fi', fasta_file, '-bed', bed_file, '-nameOnly', '-s'],
                         stdout=subprocess.PIPE)
    sites_list = seq.stdout.decode().split()
    return sites_list


def get_targets_from_site(upstream_site: str, pams: list[str]) -> tuple[list[str], list[str]]:
    """
    Find potential CRISPR target sites in a given DNA sequence, and return a list of the target sequences.

    :param upstream_site: DNA sequence of a gene TSS upstream site
    :param pams: list of PAMs by which potential sgRNA target sites will be found
    :return: a list of targets sequences that CRISPR can target
    """
    found_fwd_targets, found_rev_targets = get_sites(upstream_site, pams)
    return found_fwd_targets, found_rev_targets


def scores_batch(sgrnas: list[str], scoring_function) -> list[float]:
    """
    Calculate on target scores for sgRNAs in 'sgrnas' using the given scoring function

    :param sgrnas: list of candidate sgRNAs
    :param scoring_function: chosen on-target scoring function
    :return:
    """
    scores = scoring_function(sgrnas)
    return scores


def create_candidates(targets_lst: list[str], scores_lst: list[float], gene_id, chromosome, act_site_start, act_site_end,
                      gene_strand, targets_strand, min_gc_cont, max_gc_cont) -> list[ActivationCandidate]:
    """

    :param targets_lst: list of target sequences found in a specific gene promoter
    :param scores_lst: on-target function scores for the targets in the list of targets
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
    j = 0
    for target in targets_lst:
        target_strand = "+"
        if (gene_strand == "-" and targets_strand == "+") or (gene_strand == "+" and targets_strand == "-"):
            target_strand = "-"
        candidate = ActivationCandidate(target[:20], target[20:], gene_id, chromosome, act_site_start, act_site_end,
                                        target_strand, 0, 0, {}, scores_lst[j])
        candidate.calc_GC_content(min_gc_cont, max_gc_cont)
        candidate.calc_seq_repetitions()

        j += 1
        candidates_list.append(candidate)
    return candidates_list


def return_candidates(sites_lst: list[str], genes_filt_df, scoring_function, pams: list[str], min_gc_cont: int,
                      max_gc_cont: int) -> list[ActivationCandidate]:
    """

    :param sites_lst: list of short genes annotations and promoter sites
    :param genes_filt_df:
    :param scoring_function: chosen on-target scoring function
    :param pams: list of PAMs by which potential sgRNA target sites will be found
    :param min_gc_cont: percentage of minimum GC content by which to filter sgRNA candidates
    :param max_gc_cont: percentage of maximum GC content by which to filter sgRNA candidates
    :return: list of sgRNAs as ActivationCandidate objects
    """
    candidates_list = []
    for i in range(0, len(sites_lst), 2):
        # create gene annotations
        gene_id = sites_lst[i][1:-3]  # TODO find with regex or any other tool of validation
        chromosome = genes_filt_df[genes_filt_df['gene_id'] == gene_id]['seq_id'].values[0]
        act_site_start = int(genes_filt_df[genes_filt_df['gene_id'] == gene_id]['site_start'].values[0])
        act_site_end = int(genes_filt_df[genes_filt_df['gene_id'] == gene_id]['site_end'].values[0])
        strand = genes_filt_df[genes_filt_df['gene_id'] == gene_id]['strand'].values[0]
        # extract target sequences from each TSS upstream site
        targets_fwd, targets_rev = get_targets_from_site(sites_lst[i + 1].upper(), pams)
        # calculate on-target scores for the candidate sgRNAs
        scores_fwd = scores_batch(targets_fwd, scoring_function)
        scores_rev = scores_batch(targets_rev, scoring_function)
        # loop through targets in forward strand
        gene_candidates_list = create_candidates(targets_fwd, scores_fwd, gene_id, chromosome, act_site_start, act_site_end,
                                             strand, "+", min_gc_cont, max_gc_cont)
        gene_candidates_list += create_candidates(targets_rev, scores_rev, gene_id, chromosome, act_site_start, act_site_end,
                                             strand, "-", min_gc_cont, max_gc_cont)
        candidates_list += sorted(gene_candidates_list, key=lambda c: (c.gc_content_category, -sum(c.nucleotide_repetitions.values()), -c.on_score))[:10]
    return candidates_list


def sgrnas_for_genes(out_path: str, genes_ids: list[str], fasta_file: str, gff_file: str, scoring_function: str,
                     pams: list, bp_upstream: int, bp_downstream: int, min_gc_cont: int = 45,
                     max_gc_cont: int = 60) -> list[ActivationCandidate]:
    """
    Given a list of genes ensembl IDs, a FASTA file of full or partial genome sequence, a GFF file of annotations
    of the genome and a chosen scoring function - search and return sgRNA to target the wanted genes, with their sequence
    and the efficiency score as a list of ActivationCandidate objects.

    :param out_path: the directory path to which the algorithm will store the results
    :param genes_ids: path to input txt format file of the gene IDs
    :param fasta_file: path to input FASTA format file of the genome
    :param gff_file: path to input GFF format file of the genes annotations
    :param scoring_function: chosen on-target scoring function
    :param pams: list of PAMs by which potential sgRNA target sites will be found
    :param bp_upstream: number of base-pairs upstream to the TSS to search potential sgRNA target sites
    :param bp_downstream: number of base-pairs downstream to the TSS to search potential sgRNA target sites
    :param min_gc_cont: percentage of minimum GC content by which to filter sgRNA candidates
    :param max_gc_cont: percentage of maximum GC content by which to filter sgRNA candidates
    :return: list of potential sgRNA as ActivationCandidate objects to target the input genes
    """
    # choosing the scoring function
    scoring_function = choose_scoring_function(scoring_function)
    # handling the GFF annotations file
    genes_filt_df = filter_gff_file(gff_file, genes_ids, bp_upstream, bp_downstream)
    # handling the genome FASTA file
    sites = get_upstream_sites(out_path, fasta_file, genes_filt_df)
    # loop through the genes in the target sites list and create sgRNA candidates
    candidates_list = return_candidates(sites, genes_filt_df, scoring_function, pams, min_gc_cont, max_gc_cont)
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
    parser_obj.add_argument('genes_ids', type=str, default='', help='')
    parser_obj.add_argument('fasta_file', type=str, metavar='<fasta_file>', help='The path to the input FASTA file')
    parser_obj.add_argument('gff_file', type=str, metavar='<gff_file>', help='The path to the input GFF file')
    parser_obj.add_argument('--scoring_function', '-n', type=str, default='ucrispr',
                            help='the on scoring function of the targets. Optional scoring systems are: deephf, ucrispr. '
                                 'Additional scoring function may be added by the user or by request.')
    parser_obj.add_argument('--pams', type=list, default=['AGG', 'GGG', 'CGG', 'TGG'],
                            help='list of PAMs by which potential sgRNA target sites will be found')
    parser_obj.add_argument('--bp_upstream', '-u', type=int, default=200,
                            help='number of base-pairs upstream to the TSS to search potential sgRNA target sites')
    parser_obj.add_argument('--bp_downstream', '-d', type=int, default=0,
                            help='number of base-pairs downstream to the TSS to search potential sgRNA target sites')
    parser_obj.add_argument('--min_gc_cont', '-i', type=int, default=45,
                            help='percentage of minimum GC content by which to filter sgRNA candidates')
    parser_obj.add_argument('--max_gc_cont', '-x', type=int, default=60,
                            help='percentage of maximum GC content by which to filter sgRNA candidates')
    arguments = parser_obj.parse_args()
    return arguments


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parse_arguments(parser)
    sgrnas_for_genes(out_path=args.out_path,
                     genes_ids=args.genes,
                     fasta_file=args.fasta_file,
                     gff_file=args.gff_file,
                     scoring_function=args.scoring_function,
                     pams=args.pams,
                     bp_upstream=args.bp_upstream,
                     bp_downstream=args.bp_downstream,
                     min_gc_cont=args.min_gc_cont,
                     max_gc_cont=args.max_gc_cont
                     )
