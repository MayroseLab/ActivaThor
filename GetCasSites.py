"""Extracting target sequences from upstream sites"""
import regex
from typing import List, Tuple


def give_complementary(seq: str) -> str:
    """Given a DNA sequence (5' to 3') this function returns its antisense sequence (also 5' to 3'). This is used to
    find possible cut sites for CRISPR in the antisense strand of a given DNA sequence.
    :param seq: input DNA sequence
    :return: antisense sequence for the input
    """
    complementary_seq_list = []
    for i in range(len(seq)):
        if seq[len(seq) - 1 - i] == 'A':
            complementary_seq_list.append('T')
        elif seq[len(seq) - 1 - i] == 'T':
            complementary_seq_list.append('A')
        elif seq[len(seq) - 1 - i] == 'C':
            complementary_seq_list.append('G')
        elif seq[len(seq) - 1 - i] == 'G':
            complementary_seq_list.append('C')
        elif seq[len(seq) - 1 - i] == 'N':
            complementary_seq_list.append('N')
    return ''.join(complementary_seq_list)


def get_sites(upstream_site: str, pams: Tuple) -> Tuple[List[Tuple[str, int]], List[Tuple[str, int]], List[Tuple[str, int]], List[Tuple[str, int]]]:
    """
    This function is used to find CRISPR target site sequences from an input DNA sequence. Using regex this
    function searches for all the patterns of 23 letters long strings with all the PAM sequences in 'pams' in their end,
    in the sense and the antisense strands of the input DNA sequence. The function then returns a list of all the found
    potential targets.

    :param upstream_site: DNA sequence of a gene TSS upstream site
    :param pams: type of PAM
    :return: a list of targets sequences that CRISPR can target
    """
    target_len = 20
    found_fwd_targets_first = []
    found_fwd_targets_second = []
    found_rev_targets_first = []
    found_rev_targets_second = []
    # loop over different PAM's
    for i in range(len(pams)):
        target_and_pam = "." * target_len + pams[i]
        compiled = regex.compile(target_and_pam)
        found_sense_targets_first = regex.finditer(compiled, upstream_site[78:])  # 1-222 bps from TSS
        found_sense_targets_second = regex.finditer(compiled, upstream_site[:100])  # 201-300 bps from TSS
        found_antisense_targets_first = regex.finditer(compiled, give_complementary(upstream_site[:222]))  # 1-222 bps from TSS
        found_antisense_targets_second = regex.finditer(compiled, give_complementary(upstream_site[200:]))  # 201-300 bps from TSS

        found_fwd_targets_first += [(seq.group(0), 222-seq.start()) for seq in found_sense_targets_first if 'N' not in seq.group(0)]
        found_fwd_targets_second += [(seq.group(0), 300-seq.start()) for seq in found_sense_targets_second if 'N' not in seq.group(0)]
        found_rev_targets_first += [(seq.group(0), seq.start()) for seq in found_antisense_targets_first if 'N' not in seq.group(0)]
        found_rev_targets_second += [(seq.group(0), 200+seq.start()) for seq in found_antisense_targets_second if 'N' not in seq.group(0)]
    return found_fwd_targets_first, found_fwd_targets_second, found_rev_targets_first, found_rev_targets_second
