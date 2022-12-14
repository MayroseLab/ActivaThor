"""Module for finding off-targets"""
import os

import pandas as pd
from typing import Dict


def create_crispritz_input_file(list_of_candidates: list, crispritz_path: str) -> str:
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
    Args:
        pam_file_path: path to folder location

    Returns:
        pam file path
    """
    pam_file = f"{pam_file_path}/pamNGG.txt"
    if not os.path.exists(pam_file_path):
        os.makedirs(pam_file_path)
        with open(pam_file, "w") as f:
            f.write("NNNNNNNNNNNNNNNNNNNNNGG 3")
        return pam_file
    else:
        return pam_file


def run_crispritz(list_of_candidates: list, crispritz_script_path: str, crispys_output_path: str,
                  genome_by_chr_path: str,
                  pam_file_path, max_number_of_mismatches: int = 4, threads: int = 1) -> pd.DataFrame:
    """
    This function runs crispritz and returns its output
    nucleotide codes (e.g. 'NRG' matches both NGG and NAG).
    :param crispritz_script_path: path to crispritz python script path.
    :param list_of_candidates: A list of CandidateWithOffTargets objects
    :param crispys_output_path: A path containing the output of CRISPys
    :param genome_by_chr_path: path to the folder where the files of each chromosome fasta file.
    :param pam_file_path: a path to folder where pam file (created if not exist)
    :param max_number_of_mismatches: The maximum number of mismatches allowed between a sgRNA
     and a potential off-target.
    :param threads: number of threads to use
    :return: The output of crispritz as pd.DataFrame, where each row is a potential offtarget.
    """
    crispritz_path = f"{crispys_output_path}/crispritz"
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
    :return: a dictionary: sequence -> a CandidateWithOffTargets object where candidate.seq = sequence.
    """
    assert all([isinstance(element, Candidate) for element in list_of_candidates])
    sequence_to_candidate_dict = {}
    for i, candidate in enumerate(list_of_candidates):
        sequence_to_candidate_dict[candidate.seq] = list_of_candidates[i]
    return sequence_to_candidate_dict


# add to each "Candidate" its off-targets
def add_crispritz_off_targets(crispritz_results, sequence_to_candidate_dict: Dict):
    """
    This function adds all found off-targets to each CandidateWithOffTargets using the crispritz results
    :param crispritz_results: The output of crispritz as a pd datatable, where each row is a potential offtarget.
    :param sequence_to_candidate_dict: sequence -> a CandidateWithOffTargets object with the proper sequence
    :return: None
    """
    # apply the 'get_off_target' function on each row in the crispritz table results
    crispritz_results.apply(get_off_target, args=(sequence_to_candidate_dict,), axis=1)
    return


# This function will be used in an apply command' it reads crispritz results and create an offtarget object out of each one
def get_off_target(x, sequence_to_candidate_dict):
    """
    A function to use with apply on crispritz result table
    it takes a row of crispritz results and a dictionary of sequence:candidate, and make an OffTarget
    object from the crispritz results and add it to the candidate offtargets list
    Args:
        x: a row in crispritz results table
        sequence_to_candidate_dict: a dictionary of sequence:candidate
    Returns:
        none
    """
    candidate = sequence_to_candidate_dict[x['crRNA'][:20]]
    off_target = OffTarget(x['DNA'].upper(), x['Chromosome'], int(x['Position']), x['Direction'], int(x['Mismatches']))
    candidate.off_targets_list.append(off_target)
    return
