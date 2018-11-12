
import argparse

from parser.create_fmindex import *
from parser.find_seq_fmindex import *

import tools.create_fmindex as tc
import tools.find_seq_fmindex as tfs

# ###########################################################################
# Main function
def main():
    print("\n--------------------------------------------------------------")
    print("nektar.py: Tools for kmer created by Jellyfish.")
    print("This program was written by Eric Audemard")
    print("For more information, contact: eric.audemard@umontreal.ca")
    print("----------------------------------------------------------------\n")

    parser = argparse.ArgumentParser(prog='PROG')
    subparsers = parser.add_subparsers(help='sub-command help')

    # create the parser for the "create_fmindex" command
    p_create_fm = subparsers.add_parser('create_fmindex', help='Create fmindex for a fasta file')
    p_create_fm.set_defaults(func=tc.main_create_fmindex)
    get_parser_create_fmindex(p_create_fm)

    # create the parser for the "find_seq_fm" command
    p_find_seq_fm = subparsers.add_parser('find_seq_fm', help='Find if a seq(s) is(are) in an fmindex')
    p_find_seq_fm.set_defaults(func=tfs.main_find_seq_fmindex)
    get_parser_find_seq_fmindex(p_find_seq_fm)

    # recover arguments
    args = parser.parse_args()

    # set the remaining number of thread/process
    if hasattr(args, "THREAD"):
        args.THREAD = args.THREAD - 1

    # execute the command
    args.func(args)

main()
