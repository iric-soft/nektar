from parser.common import *


def get_parser_create_fmindex(parser):
    parser.add_argument("-f",
                        dest="FASTA",
                        help="Fasta file",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-o",
                        dest="PREFIX_OUT",
                        help="Prefix of output file")
    parser.add_argument("-b",
                        dest="BUCKET",
                        help="Bucket size (Default: 128)",
                        nargs='?',
                        default="128",
                        type=int)
    parser.add_argument("-m",
                        dest="MARK",
                        help="Mark frequence (Default: 32)",
                        nargs='?',
                        default="32",
                        type=int)

    parser.add_argument("-prot",
                        dest="PROT",
                        help="switch on AA sequence",
                        action='store_true')

    parser.add_argument("-r",
                        dest="RANDOM",
                        help="replace N by a random Nuc (Default: True)",
                        action='store_false')

    parser.add_argument("-rev",
                        dest="REVERSE",
                        help="disable the reverse complement of the fasta file (Default: False)",
                        action='store_false')

    parser.set_defaults(PROT=False)
    parser.set_defaults(RANDOM=True)
    parser.set_defaults(REVERSE=True)
