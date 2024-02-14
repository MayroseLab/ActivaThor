from typing import List
from pybedtools import BedTool


def create_off_targets_bed(list_of_candidates, lower_intersect_limit: int = 10,
                           upper_intersect_limit: int = 20) -> BedTool:
    """
    :param list_of_candidates: A list of CandidateWithOffTargets objects.
    :param lower_intersect_limit: The lower limit for the intersection between the off-target and a genomic region.
    :param upper_intersect_limit: The upper limit for the intersection between the off-target and a genomic region.
    :return: A BedTool object containing all off-target sequences
    """
    off_targets_bed_string = ""
    for candidate in list_of_candidates:
        for i, off_target in enumerate(candidate.off_targets_list):
            off_targets_bed_string += f"{off_target.chromosome} {off_target.start_position + lower_intersect_limit} {off_target.start_position + upper_intersect_limit} {candidate.seq} {i}\n"
    return BedTool(off_targets_bed_string, from_string=True)


def run_bedtools(list_of_candidates: List, genes_with_proms_csv_path: str,
                 lower_intersect_limit: int = 10, upper_intersect_limit: int = 20) -> BedTool:
    """
    This function assumes that the genome used for the CRISPys input and the annotation file are from
    the same database.
    :param list_of_candidates: A list of ActivationCandidate objects
    :param genes_with_proms_csv_path: A path to a csv file containing structural annotation of a given genome.
    IMPORTANT: make sure that the chromosome naming used in the annotation file and the genome file are the same.
    :param lower_intersect_limit: The lower limit for the intersection between the off-target and a genomic region.
    :param upper_intersect_limit: The upper limit for the intersection between the off-target and a genomic region.
    :return: A BedTool object containing the intersection between the off-targets and the genomic regions
    in the gff file
    """
    off_targets_bed = create_off_targets_bed(list_of_candidates, lower_intersect_limit, upper_intersect_limit)
    annotation_bed = BedTool(genes_with_proms_csv_path)
    intersect_bed = off_targets_bed.intersect(annotation_bed, wb=True, stream=True).saveas()
    return intersect_bed


def add_genomic_regions_to_off_targets(intersect_bed, sequence_to_candidate_dict):
    """
    :param intersect_bed: A BedTool object containing the intersection between the off-targets and the genomic regions
    in the gff file
    :param sequence_to_candidate_dict: sequence -> a CandidateWithOffTargets object with the proper sequence
    :return:
    """
    for interval in intersect_bed:
        sgRNA_candidate = sequence_to_candidate_dict[interval.name]
        # get the off target by the index created in the bed object
        off_target = sgRNA_candidate.off_targets_list[int(interval.score)]
        # check if off-target is a candidate
        genomic_region = (interval.fields[7], interval.fields[13])
        # add the genomic region to offtarget.genomic_region
        off_target.genomic_region = genomic_region
    return
