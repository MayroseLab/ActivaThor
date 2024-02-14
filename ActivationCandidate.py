"""ActivationCandidate class"""
from typing import Dict


class ActivationCandidate:
    """A class representing a sgRNA candidate to target CRISPR activation sites in a gene promoter."""

    def __init__(self, seq: str, pam: str, gene: str, gene_strand: str, chromosome: str, act_site_start: int, act_site_end: int, dist_from_TSS: int,
                 pam_strand: str, promoter_range_rank: int, gc_content: float, gc_content_category: int, nucleotide_repetitions: Dict, nuc_rep_score: int, on_score: float = 0):
        self.seq = seq
        self.pam = pam
        self.gene = gene
        self.gene_strand = gene_strand
        self.chromosome = chromosome
        self.act_site_start = act_site_start
        self.act_site_end = act_site_end
        self.dist_from_TSS = dist_from_TSS
        self.pam_strand = pam_strand
        self.promoter_range_rank = promoter_range_rank
        self.gc_content = gc_content
        self.gc_content_category = gc_content_category
        if len(nucleotide_repetitions) == 0:
            self.nucleotide_repetitions = {'5*A': 0, '5*C': 0, '4*G': 0, '4*T': 0}
        self.nuc_rep_score = nuc_rep_score
        self.on_score = on_score
        self.off_targets_list = []

    def __repr__(self):
        return f"seq = {self.seq}, pam = {self.pam}, gene = {self.gene}, dist_from_TSS = {self.dist_from_TSS}"

    def __lt__(self, other):
        if self.gc_content_category > other.gc_content_category:
            return True
        elif self.gc_content_category == other.gc_content_category:
            if self.on_score < other.on_score:
                return True
            else:
                return False
        else:
            return False

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
        self.nucleotide_repetitions['4*G'] = self.seq.count('GGGG')  # Most critical
        self.nucleotide_repetitions['4*T'] = self.seq.count('TTTT')
        self.nuc_rep_score = sum(self.nucleotide_repetitions.values())
        # self.nuc_rep_score = self.seq.count('GGGG')  # count only 4*G repetitions
        # print('nuc_rep_score CALCULATED BY "GGGG" COUNT ONLY!!')

    def to_dict(self):
        """Create a dictionary of the ActivationCandidate object"""
        self_dict = self.__dict__
        self_dict["nucleotide_repetitions"] = list(self.nucleotide_repetitions.items())
        if len(self.off_targets_list) >= 2:
            self_dict.update(self.off_targets_list[0].to_dict("1"))
            self_dict.update(self.off_targets_list[1].to_dict("2"))
        elif len(self.off_targets_list) == 1:
            self_dict.update(self.off_targets_list[0].to_dict("1"))
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
        self.score = -1
        self.number_of_mismatches = number_of_mismatches
        self.seq = seq
        self.start_position = start_position
        self.genomic_region = ""
        self.chromosome = chromosome
        self.strand = strand

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return f"[{self.seq}, {str(self.chromosome)}, {self.start_position}, {self.strand}, {self.number_of_mismatches}, {round(float(self.score), 4)}, {self.genomic_region}]"

    def to_dict(self, num: str):
        """Create a dictionary of the ActivationCandidate object"""

        self_dict = {f"off{num} score": self.score, f"off{num} mms": self.number_of_mismatches, f"off{num} seq": self.seq,
                     f"off{num} position": self.start_position, f"off{num} chr or gene": self.genomic_region[0],
                     f"off{num} chr num or gene ID": self.genomic_region[1]}
        return self_dict
