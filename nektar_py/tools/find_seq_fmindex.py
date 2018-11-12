import time
import sys
import collections
import os
import pandas as pd

folder_fm = os.path.dirname(os.path.realpath(__file__))
while not folder_fm.endswith('nektar'):
    folder_fm = os.path.dirname(folder_fm)
if folder_fm not in sys.path:
    sys.path.append(folder_fm)
import fmindex.FmIndex as fm


# ###########################################################################
# First functions
def load_fmindex(args):
    fmindex = None
    if args.PROT:
        fmindex = fm.FmIndex("AA8")
    else:
        fmindex = fm.FmIndex("Nuc4")
    fmindex.Load(args.FMINDEX)
    return fmindex


def kmer_check(fmindex, seq, prefix_out):
    print(time.strftime('%X') + ": kmer_check: " + seq)
    # res = fmindex.SearchSeq(seq)
    res = fmindex.SearchSeq(seq, 50)

    if prefix_out == "":
        print(seq + "\t" + str(res) + "\n")
    else:
        file_out = prefix_out + "_by_kmer.out"
        with open(file_out, 'w') as w:
            w.write(seq + "\t" + str(res) + "\n")


def add_to_disct_match(dict_match, res, id_seq):
    res = res.split(",")
    if len(res) > 1:
        res = res[1:len(res)]
        for match in res:
            match = match.split(":")[0]
            if match in dict_match:
                dict_match[match] += "," + id_seq
            else:
                dict_match[match] = id_seq


def file_check(fmindex, seq_file, prefix_out, add_col):
    print(time.strftime('%X') + ": file_check: " + seq_file)
    dict_match = collections.OrderedDict()

    new_seq_file = os.path.splitext(seq_file)[0]
    new_seq_file += "_" + add_col + "." + "tab"
    file_out_match = prefix_out + "_by_match.out"

    if seq_file.lower().endswith(".fa") or seq_file.lower().endswith(".fasta"):
        with open(seq_file, 'r') as f_in, open(new_seq_file, 'w') as f_out:
            line = f_in.readline()
            while line:
                id_seq = line[1:len(line)].rstrip()
                line = f_in.readline()
                seq = ""
                while line and line[0] != ">":
                    seq += line.rstrip()
                    line = f_in.readline()
                res = fmindex.SearchSeq(seq)
                add_to_disct_match(dict_match, res, id_seq)
                f_out.write(id_seq + "\t" + res + "\n")
    else:
        df_seq = pd.read_csv(seq_file, sep="\t")
        df_seq.insert(len(df_seq.columns.values) - 1, add_col, "")

        for index, row in df_seq.iterrows():
            res = fmindex.SearchSeq(row["SEQ"])
            df_seq[add_col][index] = res
            add_to_disct_match(dict_match, res, str(row["ID"]))

        df_seq.to_csv(new_seq_file, sep="\t", index=False)

    with open(file_out_match, 'w') as w:
        for match, ids in dict_match.items():
            num_match = len(ids.split(","))
            w.write(match + "\t" + str(num_match) + "\t" + ids + "\n")


# ###########################################################################
# Main function
def main_find_seq_fmindex(args):
    print(time.strftime('%X') + ": Start find seq in the FMindex: " + args.FMINDEX)

    boo_file = False
    boo_seq = False
    cpt_args = 0

    if args.FILE != "":
        cpt_args += 1
        boo_file = True
        if args.PREFIX_OUT == "":
            sys.stderr.write("Error: -f, need an output file (-o)\n")
            sys.exit()

    if args.SEQ != "":
        cpt_args += 1
        boo_seq = True

    if cpt_args == 0 or cpt_args > 1:
        sys.stderr.write("Error: you need to choose one of them: -f or -s\n")
        sys.exit()

    fmindex = load_fmindex(args)

    if args.PREFIX_OUT != "":
        if not os.path.exists(os.path.dirname(args.PREFIX_OUT)):
            os.makedirs(os.path.dirname(args.PREFIX_OUT))

    if boo_file:
        file_check(fmindex, args.FILE, args.PREFIX_OUT, args.ADD_COL)

    if boo_seq:
        kmer_check(fmindex, args.SEQ, args.PREFIX_OUT)

    print(time.strftime('%X') + ": End find seq in the FMindex")
