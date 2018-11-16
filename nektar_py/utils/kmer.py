import math
import sys
import numpy as np
import jellyfish
import common as uc


class Kmer(object):

    def __init__(self, seq, score=0):
        self._seq = seq.rstrip('\n')
        self._score = float(score)
        self._len = len(seq)

        self._count = None
        self._vertex = None

    @property
    def seq(self):
        return self._seq

    @property
    def score(self):
        return self._score

    @property
    def len(self):
        return self._len

    @property
    def entropy(self):
        return uc.get_entropy(self.seq)

    @property
    def switch(self):
        return float(uc.get_switch(self.seq)) / self.len

    @property
    def rev_seq(self):
        return self.get_reverse_complement()

    @property
    def count(self):
        return self._count

    @property
    def vertex(self):
        return self._vertex

    @vertex.setter
    def vertex(self, value):
        self._vertex = int(value)

    @vertex.deleter
    def vertex(self):
        self._vertex = None

    def init_count(self, num_sample):
        self._count = np.zeros(num_sample)

    def add_count(self, count, num_sample):
        if len(self._count) > num_sample:
            self._count[num_sample] = count
        else:
            sys.stderr.write("In Kmer the table count need to be init with enough number of samples.\n")
            sys.exit()

    def del_count(self):
        del self._count

    def get_count(self, num_sample):
        return self._count[num_sample]

    def set_count(self, count):
        self._count = count

    def get_count_jf(self, jf):
        res_k = jellyfish.MerDNA(self._seq)
        res_k.canonicalize()
        return jf[res_k]

    def get_str_count(self):
        str_res = ""
        it_c = iter(self._count)
        c = it_c.next()
        str_res = str(int(c))
        for c in it_c:
            str_res += "\t" + str(int(c))

        return str_res

    def get_reverse_complement(self):
        return uc.get_reverse_complement(self._seq)

    def switch_filter(self, min_switch):
        cpt_switch = self.switch

        return (cpt_switch >= min_switch)

    def set_score_entropy(self):
        self._score = self.entropy

    def entropy_filter(self, min_entropy):
        return (min_entropy < self.entropy)

    def print_kmer(self):
        return self._seq + "\t" + str(self._score) + "\n"

    def print_info(self):
        return self._seq + "\t" + str(self._score) + "\t" + str(self.entropy) +\
            "\t" + str(self.switch) + "\n"

    def get_pearson(self, kmer):
        res = np.corrcoef(self._count, kmer.count)[0, 1]
        return res

    def get_coef_expression(self, kmer):
        sum_k1 = np.sum(self._count)
        sum_k2 = np.sum(kmer.count)

        res = min(sum_k1, sum_k2) / max(sum_k1, sum_k2)
        return res

    def can_assembly(self, next_seq, dict_kmers, pearson, expression):
        next_seq_rev = uc.get_reverse_complement(next_seq)

        if next_seq in dict_kmers:
            if pearson == 0 or pearson <= self.get_pearson(dict_kmers[next_seq]):
                if expression == 0 or expression <= self.get_coef_expression(dict_kmers[next_seq]):
                    return dict_kmers[next_seq]
        if next_seq_rev in dict_kmers:
            if pearson == 0 or pearson <= self.get_pearson(dict_kmers[next_seq_rev]):
                if expression == 0 or expression <= self.get_coef_expression(dict_kmers[next_seq_rev]):
                    return dict_kmers[next_seq_rev]

        return None

    def get_assembly(self, dict_kmers, min_len=0, pearson=0, expression=0, depth=0, all_seed=True):
        sys.setrecursionlimit(1000000)
        words = dict()
        # kmer_rev = self.get_reverse_complement()
        # del dict_kmers[kmer_rev]

        # dict_visit save kmer already use to avoid endless loop
        dict_visit = dict()
        prefixs = self._get_prefix(self._seq, dict_kmers, dict_visit, pearson, expression,
                                   depth, depth, all_seed)
        del dict_visit
        dict_visit = dict()
        suffixs = self._get_suffix(self._seq, dict_kmers, dict_visit, pearson, expression,
                                   depth, depth, all_seed)

        if len(prefixs) > 0 and len(suffixs) > 0:
            for prefix, p_value in prefixs.items():
                for suffix, s_value in suffixs.items():
                    word_seq = prefix + self._seq + suffix
                    if len(word_seq) >= min_len:
                        word_score = abs(self._score) + p_value + s_value
                        words[word_seq] = self._get_assembly_score(word_score, word_seq)
        else:
            if len(suffixs) > 0:
                for suffix, s_value in suffixs.items():
                    word_seq = self._seq + suffix
                    if len(word_seq) >= min_len:
                        word_score = abs(self._score) + s_value
                        words[word_seq] = self._get_assembly_score(word_score, word_seq)
            if len(prefixs) > 0:
                for prefix, p_value in prefixs.items():
                    word_seq = prefix + self._seq
                    if len(word_seq) >= min_len:
                        word_score = abs(self._score) + p_value
                        words[word_seq] = self._get_assembly_score(word_score, word_seq)
        return words

    def _next_in_dict(self, next_seq, dict_kmers, dict_visit, pearson, expression, all_seed):
        next_seq_rev = uc.get_reverse_complement(next_seq)
        if next_seq in dict_visit or next_seq_rev in dict_visit:
            # print("kmer already used")
            return None

        find = False
        return_kmer = None
        if next_seq in dict_kmers:
            if pearson == 0 or pearson <= self.get_pearson(dict_kmers[next_seq]):
                if expression == 0 or expression <= self.get_coef_expression(dict_kmers[next_seq]):
                    # print("_next_in_dict next valide")
                    if all_seed:
                        return_kmer = dict_kmers[next_seq]
                    else:
                        return_kmer = dict_kmers.pop(next_seq)
                    dict_visit[next_seq] = None
                    find = True
        if next_seq_rev in dict_kmers:
            if find or pearson == 0 or pearson <= self.get_pearson(dict_kmers[next_seq_rev]):
                if expression == 0 or expression <= self.get_coef_expression(dict_kmers[next_seq_rev]):
                    # print("_next_in_dict next_rev valide")

                    if all_seed:
                        return_kmer = dict_kmers[next_seq_rev]
                    else:
                        return_kmer = dict_kmers.pop(next_seq_rev)
                    dict_visit[next_seq_rev] = None
                    find = True

        return return_kmer

    # I need seq to be sure to extend on same strands even if the kmer is saved on reverse
    def _get_prefix(self, seq, dict_kmers, dict_visit, pearson, expression, depth, max_depth, all_seed):
        # print("_get_prefix: " + self._seq + " depth: " + str(depth) + "_" + str(max_depth) + " visit: " + str(len(dict_visit)))
        words = dict()
        boo_find_kmer = False
        for nuc in "A" "C" "G" "T":
            next_seq = nuc + seq[0:len(seq) - 1]
            next_kmer = self._next_in_dict(next_seq, dict_kmers, dict_visit, pearson,
                                           expression, all_seed)

            if next_kmer is not None:
                boo_find_kmer = True
                prefixs = next_kmer._get_prefix(next_seq, dict_kmers, dict_visit, pearson,
                                                expression, max_depth, max_depth, all_seed)
                if len(prefixs) == 0:
                    words[nuc] = abs(self._score)
                else:
                    for prefix, p_value in prefixs.items():
                        words[prefix + nuc] = abs(self._score) + p_value

        if not boo_find_kmer and depth > 0:
            for nuc in "A" "T" "C" "G":
                next_seq = nuc + seq[0: len(seq) - 1]
                next_kmer = Kmer(next_seq)
                next_kmer.set_count(self._count)
                prefixs = next_kmer._get_prefix(next_seq, dict_kmers, dict_visit, pearson,
                                                expression, depth - 1, max_depth, all_seed)
                if len(prefixs) > 0:
                    for prefix, p_value in prefixs.items():
                        words[prefix + nuc] = p_value + abs(self._score)
                del next_kmer

        return words

    # I need seq to be sure to extend on same strands even if the kmer is saved on reverse
    def _get_suffix(self, seq, dict_kmers, dict_visit, pearson, expression, depth, max_depth, all_seed):
        # print("_get_suffix : " + self._seq + " depth: " + str(depth) + "_" + str(max_depth) + " visit: " + str(len(dict_visit)))
        words = dict()
        boo_find_kmer = False
        for nuc in "A" "C" "G" "T":
            next_seq = seq[1:len(seq)] + nuc
            next_kmer = self._next_in_dict(next_seq, dict_kmers, dict_visit, pearson,
                                           expression, all_seed)

            if next_kmer is not None:
                boo_find_kmer = True
                suffixs = next_kmer._get_suffix(next_seq, dict_kmers, dict_visit, pearson,
                                                expression, max_depth, max_depth, all_seed)
                if len(suffixs) == 0:
                    words[nuc] = abs(self._score)
                else:
                    for suffix, s_value in suffixs.items():
                        words[nuc + suffix] = s_value + abs(self._score)

        if not boo_find_kmer and depth > 0:
            for nuc in "A" "T" "C" "G":
                next_seq = seq[1:len(seq)] + nuc
                next_kmer = Kmer(next_seq)
                next_kmer.set_count(self._count)
                suffixs = next_kmer._get_suffix(next_seq, dict_kmers, dict_visit, pearson,
                                                expression, depth - 1, max_depth, all_seed)
                if len(suffixs) > 0:
                    for suffix, s_value in suffixs.items():
                        words[nuc + suffix] = s_value + abs(self._score)
                del next_kmer

        return words

    def _get_assembly_score(self, word_score, word_seq):
        return word_score / (len(word_seq) - self.len + 1)

