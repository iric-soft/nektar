from parser.common import *


def get_parser_find_seq_fmindex(parser):
    parser.add_argument("-fm",
                        dest="FMINDEX",
                        help="FMindex file of the ref fasta",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-f",
                        dest="FILE",
                        help="Fasta or tabulated (with a column SEQ) file of sequences/kmers",
                        nargs='?',
                        default="",
                        type=str)
    parser.add_argument("-add",
                        dest="ADD_COL",
                        help="Name of the added column in output tabulated file (Default:FMINDEX)",
                        default="FMINDEX",
                        type=str)

    parser.add_argument("-s",
                        dest="SEQ",
                        help="Sequence/kmer to be checked",
                        nargs='?',
                        default="",
                        type=str)
    parser.add_argument("-o",
                        dest="PREFIX_OUT",
                        help="Output prefix file",
                        nargs='?',
                        default="",
                        type=str)

    parser.add_argument("-prot",
                        dest="PROT",
                        help="switch on AA sequence",
                        action='store_true')

    parser.set_defaults(PROT=False)
