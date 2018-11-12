import time
import sys
import os

folder_fm = os.path.dirname(os.path.realpath(__file__))
while not folder_fm.endswith('nektar'):
    folder_fm = os.path.dirname(folder_fm)
if folder_fm not in sys.path:
    sys.path.append(folder_fm)
import fmindex.FmIndex as fm


# ###########################################################################
# First functions
def create_fmindex(args):

    fmindex = None
    if args.PROT:
        fmindex = fm.FmIndex(args.FASTA, "AA8", args.RANDOM, args.REVERSE)
    else:
        fmindex = fm.FmIndex(args.FASTA, "Nuc4", args.RANDOM, args.REVERSE)

    fmindex.CreateIndex(args.BUCKET, args.MARK)
    fmindex.Save(args.PREFIX_OUT + "_" + str(args.BUCKET) + "-" + str(args.MARK) + ".fm")
    del fmindex


# ###########################################################################
# Main function
def main_create_fmindex(args):
    print(time.strftime('%X') + ": Start create FM index on:" + args.FASTA)

    create_fmindex(args)

