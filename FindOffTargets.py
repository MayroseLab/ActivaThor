"""Module for finding off-targets"""
import os
import sys

import pandas as pd
from typing import Dict, List

from ActivationCandidate import ActivationCandidate, OffTarget
from OffTargetsAnnotations import add_genomic_regions_to_off_targets, run_bedtools
from ScoringFunctions import moff


def create_crispritz_input_file(list_of_candidates: List, crispritz_path: str) -> str:
    """
    :param list_of_candidates: A list of CandidateWithOffTargets objects
    :param crispritz_path: A path to the crispys result folder where a folder for crispritz will be created
    :return: A path to the input for xxx and will write the input file to xxx
    """
    crispritz_infile = os.path.join(crispritz_path, 'crispritz_infile.txt')
    out = ""
    for candidate in list_of_candidates:  # go over each candidate and get the guide sequence
        out += f"{candidate.seq}NNN\n"
    with open(crispritz_infile, 'w') as f:
        f.write(out)
    return crispritz_infile


def create_pam_file(pam_file_path: str) -> str:
    """
    create the pam file for crispritz (for NGG pam)

    :param pam_file_path: path to folder location
    :return: pam file path
    """
    pam_file = f"{pam_file_path}/pamNGG.txt"
    if not os.path.exists(pam_file_path):
        os.makedirs(pam_file_path)
        with open(pam_file, "w") as f:
            f.write("NNNNNNNNNNNNNNNNNNNNNGG 3")
        return pam_file
    else:
        return pam_file


def run_crispritz(list_of_candidates: List, crispritz_script_path: str, output_path: str,
                  genome_by_chr_path: str,
                  pam_file_path: str, max_number_of_mismatches: int = 4, threads: int = 1) -> pd.DataFrame:
    """
    This function runs crispritz and returns its output
    nucleotide codes (e.g. 'NRG' matches both NGG and NAG).

    :param crispritz_script_path: path to crispritz python script path.
    :param list_of_candidates: A list of CandidateWithOffTargets objects
    :param output_path: A path containing the output of CRISPys
    :param genome_by_chr_path: path to the folder where the files of each chromosome fasta file.
    :param pam_file_path: a path to folder where pam file (created if not exist)
    :param max_number_of_mismatches: The maximum number of mismatches allowed between a sgRNA
     and a potential off-target.
    :param threads: number of threads to use
    :return: The output of crispritz as pd.DataFrame, where each row is a potential offtarget.
    """
    crispritz_path = f"{output_path}/crispritz"
    os.makedirs(crispritz_path, exist_ok=True)
    # create crispritz input file
    crispritz_infile = create_crispritz_input_file(list_of_candidates, crispritz_path)
    pam_file = create_pam_file(pam_file_path)
    # run crispritz
    os.system(
        f"python {crispritz_script_path} search {genome_by_chr_path}/ {pam_file} {crispritz_infile} {crispritz_path}/crispritz -mm {max_number_of_mismatches} -r -th {threads}")
    # get results
    crispritz_results = pd.read_csv(f"{crispritz_path}/crispritz.targets.txt", sep="\t")
    print(f"Number of off-targets: {crispritz_results.shape[0]}")
    return crispritz_results


# create dictionary of sequence:candidate
def create_sequence_to_candidate_dict(list_of_candidates):
    """
    This function takes a list of sgRNA candidates, and returns a sequence to candidate dictionary.

    :param list_of_candidates: a list of sgRNA candidates
    :return: a dictionary: sequence -> an ActivationCandidate object where candidate.seq = sequence.
    """
    assert all([isinstance(element, ActivationCandidate) for element in list_of_candidates])
    sequence_to_candidate_dict = {}
    for i, candidate in enumerate(list_of_candidates):
        sequence_to_candidate_dict[candidate.seq] = list_of_candidates[i]
    return sequence_to_candidate_dict


# add to each "Candidate" its off-targets
def add_crispritz_off_targets(crispritz_results, sequence_to_candidate_dict: Dict[str, ActivationCandidate]):
    """
    This function adds all found off-targets to each CandidateWithOffTargets using the crispritz results.

    :param crispritz_results: The output of crispritz as a pd datatable, where each row is a potential offtarget.
    :param sequence_to_candidate_dict: sequence -> a CandidateWithOffTargets object with the proper sequence
    """
    # apply the 'get_off_target' function on each row in the crispritz table results
    crispritz_results.apply(get_off_target, args=(sequence_to_candidate_dict,), axis=1)
    return


# This function will be used in an apply command' it reads crispritz results and create an offtarget object out of each one
def get_off_target(x, sequence_to_candidate_dict: Dict[str, ActivationCandidate]):
    """
    A function to use with apply on crispritz result table
    it takes a row of crispritz results and a dictionary of sequence:candidate, and make an OffTarget
    object from the crispritz results and add it to the candidate offtargets list.

    :param x: a row in crispritz results table
    :param sequence_to_candidate_dict: a dictionary of sequence:candidate
    """
    candidate = sequence_to_candidate_dict[x['crRNA'][:20]]
    off_target = OffTarget(x['DNA'].upper(), candidate.chromosome, int(x['Position']), x['Direction'], int(x['Mismatches']))
    legit_letters = True
    for char in off_target.seq:
        if char not in {"A", "C", "T", "G"}:
            legit_letters = False
            break
    if legit_letters:
        candidate.off_targets_list.append(off_target)
    return


def calculate_scores(candidates_list: List[ActivationCandidate]):
    """
    Calculate the off-target scores for each off-target of each candidate and store the scores in each off-target's
    score attribute.

    :param candidates_list: a list of sgRNA candidates
    """
    batch_off_targets_list = []
    batch_candidates_list = []
    for candidate in candidates_list:
        batch_candidates_list += [candidate.seq + candidate.pam for _ in range(len(candidate.off_targets_list))]
        batch_off_targets_list += [off_target.seq for off_target in candidate.off_targets_list]
    scores = moff(batch_off_targets_list, batch_candidates_list)
    for i in range(len(candidates_list)):
        for j in range(len(candidates_list[i].off_targets_list)):
            candidates_list[i].off_targets_list[j].score = round(scores[i+j], 4)
        # sort the off-targets in the off_targets_list by score from highest to lowest
        candidates_list[i].sort_off_targets()


def get_off_targets(candidates_list: List[ActivationCandidate], in_path: str, out_path: str, gff_with_proms_path: str):
    """
    Find the off-targets for each candidate using crispritz, store them in the candidate's off_targets_list attribute
    as a list of OffTarget objects. Then calculates the off-target scores for each off-target of each candidate and
    store the scores in each off-target's score attribute.

    :param candidates_list: a list of sgRNA candidates
    :param in_path: path to directory in which the chromosome FASTA files are stored
    :param out_path: output path for the algorithm results
    :param gff_with_proms_path:
    """
    # run crispritz
    off_targets_pd = run_crispritz(candidates_list, sys.exec_prefix + "/bin/crispritz.py", out_path,
                                   in_path + "/chromosomes", out_path + "/pams")  # TODO fix the paths
    # create a dictionary of sequence -> candidate
    sequence_to_candidate_dict = create_sequence_to_candidate_dict(candidates_list)
    # add the found off-targets of each candidate to the candidate's off_targets_list
    add_crispritz_off_targets(off_targets_pd, sequence_to_candidate_dict)
    #
    intersect_bed = run_bedtools(candidates_list, out_path+gff_with_proms_path)
    # add the genomic regions to the OffTargets attributes
    add_genomic_regions_to_off_targets(intersect_bed, sequence_to_candidate_dict)
    # calculate the off-target scores for each off_target of each candidate
    calculate_scores(candidates_list)
