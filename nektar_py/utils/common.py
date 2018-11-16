import time
import csv
import sys
import gc
import os
import math
import collections
import jellyfish

import numpy as np

import utils.sample as us
import utils.kmer as uk

from string import maketrans


# ###########################################################################
# Read functions
def read_histo(name_file):
    sum = 0
    with open(name_file, 'r') as f:
        for line in f:
            p_line = line.split()
            sum += int(p_line[1])

    return sum


def read_project(project_file, samples):
    print(time.strftime('%X') + ": Read project file: " + project_file)
    exist_file = True

    with open(project_file) as cvsfile:
        project_csv = csv.DictReader(cvsfile, delimiter='\t')

        for line in project_csv:
            name_sample = line["SAMPLE"]
            total_kmer = int(line["KMER"])
            group_sample = line.get("GROUP", None)

            files = dict()
            for format in ("JF", "FASTQ1", "FASTQ2"):
                name_file = line.get(format, None)
                if name_file is not None and name_file.lower() != "none":
                    if not os.path.isfile(name_file):
                        sys.stderr.write("In read_project function, file doen't exist: "
                                         + str(name_file) + "\n")

                        exist_file = False

                    files[format] = name_file
                else:
                    files[format] = None

            samples.append(us.Sample(name_sample, group_sample, total_kmer, files))

    if exist_file is False:
        sys.exit()


def read_kmer(file_name):
    print(time.strftime('%X') + ": Read kmer file: " + file_name)
    list_kmer = []

    with open(file_name) as file_r:
        # Check if the first line is not an header
        line = file_r.readline().rstrip('\n').split("\t")
        if len(line[0]) != 0 and line[0].lower() != "kmer":
            l_split = line
            score = l_split[1] if len(l_split) > 1 else "0"
            list_kmer.append(uk.Kmer(l_split[0], score))

        for line in file_r:
            l_split = line.split("\t")
            score = l_split[1] if len(l_split) > 1 else "0"
            list_kmer.append(uk.Kmer(l_split[0], score))

    return list_kmer


def read_fasta(fasta_file):
    dic_seqs = collections.OrderedDict()
    with open(fasta_file, 'r') as r:
        name = ""
        seq = ""
        for line in r:
            if line[0] == ">":
                if name != "":
                    dic_seqs[name] = seq
                    seq = ""
                name = line.split(">")[1].rstrip('\n')
            else:
                seq += line.rstrip('\n')

    return dic_seqs


# ###########################################################################
# Entropy functions
def entropy_ideal(length):
    prob = 1.0 / length
    return -1.0 * length * prob * math.log(prob) / math.log(2.0)


# ###########################################################################
# get functions

# kmer
def get_path_kmer_split_file(path_db, name_sample, num_split, max=False):
    file_name = ""
    if max:
        file_name = name_sample + "_" + num_split + "_max.out"
    else:
        file_name = name_sample + "_" + num_split + ".out"

    dir = os.path.join(path_db, name_sample)
    if not os.path.exists(dir):
        os.makedirs(dir)

    return os.path.join(dir, file_name)


def get_kmer_from_seq(seq, lkmer):
    dict_kmers = [collections.OrderedDict((seq[i:i + lkmer], uk.Kmer(seq[i:i + lkmer]))) for i in range(len(seq) - lkmer)]

    return dict_kmers


def get_dict_kmers_from_seq(seq, lkmer):
    dict_kmers = {seq[i:i + lkmer]: uk.Kmer(seq[i:i + lkmer]) for i in range(len(seq) - lkmer)}

    return dict_kmers


def get_list_kmers_from_seq(seq, lkmer):
    list_kmers = [uk.Kmer(seq[i:i + lkmer]) for i in range(len(seq) - lkmer)]

    return list_kmers


def get_count(kmer, jf):
    res_k = jellyfish.MerDNA(kmer.seq)
    res_k.canonicalize()

    return jf[res_k]


def get_count_project(project, list_kmers):
    cpt_sample = 0
    it_sample = project.it_samples
    sample = it_sample.next()
    print(time.strftime('%X') + " " + sample.name)
    cpt_sample = 0
    # Ticky: loading all the file in cache
    os.system("wc -l " + sample.jf_file)
    jf = jellyfish.QueryMerFile(sample.jf_file)
    cpt_kmer = 0
    for kmer in list_kmers:
        cpt_kmer += 1
        kmer.init_count(project.num_samples)
        kmer.add_count(kmer.get_count_jf(jf), cpt_sample)
    del jf

    for sample in it_sample:
        print(time.strftime('%X') + " " + sample.name)
        cpt_sample += 1
        # Ticky: loading all the file in cache
        os.system("wc -l " + sample.jf_file)
        jf = jellyfish.QueryMerFile(sample.jf_file)
        for kmer in list_kmers:
            kmer.add_count(kmer.get_count_jf(jf), cpt_sample)
        del jf

    gc.collect()
    print(time.strftime('%X') + ": Get count done!")


# Group
def get_pairs_group(ids_group):
    pairs_group = []
    for i in range(len(ids_group)):
        for j in range(i + 1, len(ids_group)):
            pairs_group.append([ids_group[i], ids_group[j]])

    return pairs_group


def get_dir_pair(pair_a, pair_b):
    return str(pair_a) + "_vs_" + str(pair_b)


# Sequence/fasta
def get_seqs_from_fasta(fasta_file):
    seqs = dict()
    with open(fasta_file, 'r') as r:

        name = ""
        seq = ""
        for line in r:
            line = line.rstrip('\n')
            if len(line) > 0 and line[0] == ">":
                if len(name) > 0:
                    seqs[name] = seq
                name = line.split(">")[1]
                seq = ""
            else:
                seq = seq + line

        if len(name) > 0:
            seqs[name] = seq

    return seqs


def get_seq_from_fasta(fasta_file, chr, pos, window):
    seq = ""
    with open(fasta_file, 'r') as r:
        start = max(0, pos - window)
        end = pos + window + 1  # +1 to select the las bp we want (line[x:y+1])

        name = ""
        found_chr = False
        found_pos = False
        cpt_bp = 0
        for line in r:
            line = line.rstrip('\n')
            if line[0] == ">":
                name = line.split(">")[1]
                if name == chr:
                    found_chr = True
                else:
                    found_chr = False
            else:
                if found_chr:
                    cpt_bp += len(line)

                    if cpt_bp >= start:
                        found_pos = True
                        p1 = (len(line) - 1) - (cpt_bp - start)
                        p1 = max(0, p1)

                        p2 = (len(line) - 1) - (cpt_bp - end)
                        p2 = min(p2, len(line))

                        seq = seq + line[p1:p2]

                    if cpt_bp >= end:
                        return seq

            if found_pos and not found_chr:
                return seq

    if seq == "":
        sys.stderr.write("In extract_seq_fasta function, the position (" +
                         str(chr) + " " + str(pos) + ") isn't found.\n")
        sys.exit()

def get_reverse_complement (seq, tb=None):
   '''
   Complements a DNA sequence, returning the reverse complement. Faster than get_reverse_complement_seq().
   '''
   if tb is None:
        tb = maketrans("ACGTRYMKWSBDHVNacgtrymkwsbdhvn",
                       "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn")        
   return seq[-1::-1].translate(tb)

def get_entropy(seq):
    # http://stackoverflow.com/questions/2979174/how-do-i-compute-the-approximate-entropy-of-a-bit-string
    prob = [float(seq.count(c)) / len(seq) for c in dict.fromkeys(list(seq))]
    entropy = - sum([p * math.log(p) / math.log(2.0) for p in prob])

    return entropy


def get_switch(seq):
    prev_char = seq[0]
    cpt_switch = 0.0
    for char in seq:
        if char != prev_char:
            prev_char = char
            cpt_switch += 1

    return cpt_switch


# ###########################################################################
# write functions
def write_seq(project, name_seq, seq, dict_kmer, args, w_tab, num_graph="NA"):
    list_kmer = [seq[i:i + args.LKMER] for i in range(len(seq) - args.LKMER + 1)]

    cpt_line = 1
    for seq_kmer in list_kmer:
        kmer = None
        if seq_kmer in dict_kmer:
            kmer = dict_kmer[seq_kmer]
        else:
            kmer = dict_kmer[get_reverse_complement(seq_kmer)]

        line = str(name_seq) + "\t" + str(num_graph) + "\t" + seq_kmer + "\t" + str(kmer.entropy) + "\t" + str(kmer.switch) + "\t" + kmer.get_str_count()

        cpt_line += 1
        w_tab.write(line + "\n")


def write_kmer(project, dict_seqs, prefix, args, path_dir):
    boo_header = False

    if not os.path.exists(path_dir):
        os.makedirs(path_dir)

    file_out_tab = os.path.join(path_dir, str(prefix) + "_count.tab")
    with open(file_out_tab, 'w') as w_tab:
        for name, seq in list(dict_seqs.items()):
            list_kmer = [seq[i:i + args.LKMER] for i in range(len(seq) - args.LKMER)]

            list_line = [""] * (len(list_kmer) + 1)
            list_line[0] = "ID_ASSEMBLY\tKMER\tENTROPY\tSWITCH"
            cpt_line = 1
            for kmer in list_kmer:
                list_line[cpt_line] = str(name) + "\t" + kmer.seq + "\t" +\
                    str(kmer.entropy) + "\t" + str(kmer.switch)

                cpt_line += 1

            file_out_r = os.path.join(path_dir, str(name) + "_count_ggplot.csv")
            with open(file_out_r, 'w') as w_ggplot:
                w_ggplot.write("POS,KMER,SAMPLE,GROUP,LOG_COUNT,COUNT\n")

                for sample in project.samples:
                    list_line[0] = list_line[0] + "\t" + sample.name
                    jf = jellyfish.QueryMerFile(sample.jf_file)

                    cpt_line = 1
                    for kmer in list_kmer:
                        count = get_count(kmer, jf)
                        log_count = np.log10(count * args.LOG_F + args.LOG_C) /\
                            np.log10(sample.num_kmer + args.LOG_C)

                        w_ggplot.write(str(cpt_line) + "," + kmer.seq + "," +
                                       sample.name + "," + sample.group + "," +
                                       str(log_count) + "," + str(count) + "\n")

                        list_line[cpt_line] = list_line[cpt_line] + "\t" + str(count)
                        cpt_line += 1

            if boo_header:
                del list_line[0]
            else:
                boo_header = True

            for line in list_line:
                w_tab.write(line + "\n")
            del list_line

#            if (args.GGPLOT):
#                file_out_gg = os.path.join(path_dir, str(name) + "_count_ggplot.pdf")
#                data = pd.read_csv(file_out_r, sep=",")
#                gg = ggplot(aes(x="POS", y="LOG_COUNT", color="GROUP", group="GROUP"), data=data)
#                gg = gg + geom_point(alpha=0.5)
#                # gg = gg + stat_smooth()
#                gg = gg + stat_smooth(se=False, size=3)
#                width = len(data["POS"]) * 10 / 3000
#                width = min(23, width)
#                width = max(8, width)
#                ggsave(filename=file_out_gg, plot=gg, width=width)
