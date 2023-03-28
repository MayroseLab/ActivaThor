"""ActivationCandidate class"""
from typing import Dict


class ActivationCandidate:
    """A class representing a sgRNA candidate to target CRISPR activation sites in a gene promoter."""

    def __init__(self, seq: str, pam: str, gene: str, chromosome: str, act_site_start: int, act_site_end: int,
                 strand: str, promoter_range_rank: int, gc_content: float, gc_content_category: int, nucleotide_repetitions: Dict, nuc_rep_score: int, on_score: float = 0):
        self.seq = seq
        self.pam = pam
        self.gene = gene
        self.chromosome = chromosome
        self.act_site_start = act_site_start
        self.act_site_end = act_site_end
        self.strand = strand
        self.promoter_range_rank = promoter_range_rank
        self.gc_content = gc_content
        self.gc_content_category = gc_content_category
        if len(nucleotide_repetitions) == 0:
            self.nucleotide_repetitions = {'5*A': 0, '5*C': 0, '4*G': 0, '4*T': 0}
        self.nuc_rep_score = nuc_rep_score
        self.on_score = on_score
        self.off_targets_list = []

    def __repr__(self):
        return f"seq = {self.seq}, gene = {self.gene}, on-score = {self.on_score}"

    def calc_GC_content(self, min_gc_cont: int, max_gc_cont: int):
        """
        Calculate the percentage of G and C nucleotides in the sgRNA candidate sequence, and calculate the GC content
        category.

        :param min_gc_cont: percentage of minimum GC content by which to filter sgRNA candidates
        :param max_gc_cont: percentage of maximum GC content by which to filter sgRNA candidates
        """
        gc_num = self.seq.count("G") + self.seq.count("C")
        self.gc_content = round(gc_num / 20, 2)
        gc_cont_category = 1
        for i in range(10):
            if (min_gc_cont - 5 * i)/100 <= self.gc_content <= (max_gc_cont + 5 * i)/100:
                self.gc_content_category = gc_cont_category
                break
            gc_cont_category += 1

    def calc_seq_repetitions(self):
        """
        Fill the nucleotide repetitions dictionary
        """
        self.nucleotide_repetitions['5*A'] = self.seq.count('AAAAA')
        self.nucleotide_repetitions['5*C'] = self.seq.count('CCCCC')
        self.nucleotide_repetitions['4*G'] = self.seq.count('GGGG')
        self.nucleotide_repetitions['4*T'] = self.seq.count('TTTT')
        self.nuc_rep_score = sum(self.nucleotide_repetitions.values())

    def to_dict(self):
        """Create a dictionary of the ActivationCandidate object"""
        self_dict = self.__dict__
        if len(self.off_targets_list) >= 2:
            self_dict["off1"] = f"{self.off_targets_list[0].score}:{self.off_targets_list[0].number_of_mismatches}:{self.off_targets_list[0].seq}:{self.off_targets_list[0].start_position}:{self.off_targets_list[0].genomic_region}"
            self_dict["off2"] = f"{self.off_targets_list[1].score}:{self.off_targets_list[1].number_of_mismatches}:{self.off_targets_list[1].seq}:{self.off_targets_list[1].start_position}:{self.off_targets_list[1].genomic_region}"
        self_dict["nucleotide_repetitions"] = list(self.nucleotide_repetitions.items())
        return self_dict

    def sort_off_targets(self):
        """
        Sort the off_targets_list by the off-targets scores from highest to lowest
        """
        self.off_targets_list.sort(key=lambda off: -off.score)


class OffTarget:
    """
    This class contains off-target information for an sgRNA candidate, it is intended to be part of a
    ActivationCandidate object as an item of the off_target_list
    """
    def __init__(self, seq: str, chromosome: str, start_position: int, strand: str, number_of_mismatches: int):
        self.seq = seq
        self.chromosome = chromosome
        self.start_position = start_position
        self.strand = strand
        self.number_of_mismatches = number_of_mismatches
        self.genomic_region = ""
        self.score = -1

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return f"[{self.seq}, {str(self.chromosome)}, {self.start_position}, {self.strand}, {self.number_of_mismatches}, {round(float(self.score), 4)}, {self.genomic_region}]"
