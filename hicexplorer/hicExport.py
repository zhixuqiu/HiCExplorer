from __future__ import division
from os.path import splitext
import sys
import argparse
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
from scipy import sparse


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Saves a matrix in .npz format using a plain text format.')

    # define the arguments
    parser.add_argument('--inFile', '-in',
                        help='input file(s). Could be one or many files. '
                        'Multiple input files are allowed for hicexplorer or lieberman format. '
                        ' In case of multiple input files, they will be combined. ',
                        nargs='+',
                        required=True)

    parser.add_argument('--inputFormat',
                        help='file format for input file. \n'
                             '(options : npz, lieberman)',
                        default='npz')

    parser.add_argument('--chrNameList',
                        help='list of chromosome names (only if input format is lieberman), eg : 1 2 .',
                        nargs='+',
                        )

    parser.add_argument('--outFileName', '-o',
                        help='File name to save the plain text matrix. In the case of "lieberman" '
                             'output format this should be the path of a folder where the information '
                             'per chromosome is stored.',
                        required=True)

    parser.add_argument('--chromosomeOrder',
                        help='Chromosomes and order in which the '
                             'chromosomes should be plotted. This option '
                             'overrides --region and --region2 ',
                        nargs='+')

    parser.add_argument('--bplimit', '-b',
                        help='When merging many matrices : maximum limit (in base pairs) after '
                              'which the matrix will be truncated. i.e. TADs bigger than this '
                              'size will not be shown. For Matrices with very high resolution, '
                              'truncating the matrix after a limit helps in saving memory '
                              'during processing, without much loss of data. You can use '
                              'bplimit of 2 x size of biggest expected TAD. ',
                        type=int,
                        metavar='INT bp',
                        default=None)

    parser.add_argument('--outputFormat',
                        help='Output format. The possibilities are "dekker",  "ren" and "hicexplorer". '
                             'The dekker format outputs the whole matrix where the '
                             'first column and first row are the bin widths and labels. '
                             'The "ren" format is a list of tuples of the form '
                             'chrom, bin_star, bin_end, values. '
                             'The lieberman format writes separate files for each chromosome,'
                             'with three columns : contact start, contact end, and raw observed score. '
                             'This corresponds to the RawObserved files from lieberman group. ',
                        default='dekker',
                        choices=['dekker', 'ren', 'lieberman', 'hicexplorer'])

    parser.add_argument('--clearMaskedBins',
                        help='if set, masked bins are removed from the matrix. Masked bins '
                             'are those that do not have any values, mainly because they are'
                             'repetitive regions of the genome',
                        action='store_true')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def main():
    args = parse_arguments().parse_args()

    ## create hiC matrix with given input format
    # additional file needed for lieberman format
    if args.inputFormat == 'lieberman':
        if (args.chrNameList is None ):
            exit("Error: --chrNameList is required when the input format is lieberman. ")
        else:
            hic_ma = hm.hiCMatrix(matrixFile=args.inFile, format='lieberman', chrnameList=args.chrNameList)

    elif args.inputFormat == 'npz' and len(args.inFile) > 1: # assume npz_multi format
        if (args.bplimit is None ):
            print """\nNo limit given for maximum depth. Using all available data. \n
                     Set an appropriate resolution limit  the file doesn't save"""
        else:
            print "\nCutting maximum matrix depth to {} for saving".format(args.bplimit)

        hic_ma = hm.hiCMatrix(matrixFile=args.inFile, format='npz_multi', bplimit=args.bplimit)

    else:
        hic_ma = hm.hiCMatrix(matrixFile=args.inFile[0], format='npz')

    if args.chromosomeOrder:
        hic_ma.keepOnlyTheseChr(args.chromosomeOrder)

    if args.clearMaskedBins:
        hic_ma.maskBins(hic_ma.nan_bins)

    sys.stderr.write('saving...\n')

    if args.outputFormat == 'dekker':
        hic_ma.save_dekker(args.outFileName)
    elif args.outputFormat == 'ren':
        hic_ma.save_bing_ren(args.outFileName)
    elif args.outputFormat == 'lieberman':
        hic_ma.save_lieberman(args.outFileName)
    else:
        hic_ma.save(args.outFileName)
