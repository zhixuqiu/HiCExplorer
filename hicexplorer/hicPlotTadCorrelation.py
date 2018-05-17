#!/usr/bin/env python

import os
import sys
import argparse

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     description='Plots TAD correlations with a matrix computed by deeptools.')

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument("--deeptools_matrix",
                                "-m",
                                type=str,
                                help="A matrix obtained from deeptools computeMatrix scale-regions "
                                "with padded regions on both side. The input for the matrix is a TADboundery.bed "
                                "from hicexplorer and an active_histonmark.bw",
                                required=True)

    parserRequired.add_argument("--outFileName", "-o",
                                type=str,
                                help="File name to save the image.",
                                required=True)
    parserRequired.add_argument("--tadRegion", "-ts",
                                type=str,
                                help="Tad region; format is: chr:start-end",
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--dpi',
                           help='Optional parameter: Resolution for the image in case the'
                           'output is a raster graphics image (e.g. png, jpg)',
                           type=int,
                           default=300)
    parserOpt.add_argument('--colorMap',
                           help='Color map to use for the heatmap. Available '
                           'values can be seen here: '
                           'http://matplotlib.org/examples/color/colormaps_reference.html',
                           default='RdYlGn')
    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)
    table = pd.read_csv(args.deeptools_matrix, sep='\t', skiprows=1, header=None)
    values = table.iloc[:, 6:]
    values = np.array(values)
    values[np.isnan(values)] = 0

    matrix_size = values.shape[1]
    correlation_table = np.zeros((matrix_size, matrix_size))
    # log.info('table {}'.format(table[:, :4]))
    for col1 in range(matrix_size):
        for col2 in range(matrix_size):
            if col1 >= col2:
                correlation_table[col1, col2] = spearmanr(values[:, col1], values[:, col2])[0]

            # if col1 < 10 and col2 < 1:
            #     print(values[:, col1])
    final_correlation_table = correlation_table + correlation_table.T - np.diag(correlation_table.diagonal())
    heatmap = plt.figure(figsize=(6.4, 4.8))

    ax = heatmap.add_subplot(111)
    cax = ax.imshow(final_correlation_table, interpolation='nearest', vmin=0, vmax=1, cmap=args.colorMap)
    cbar = plt.colorbar(cax, ticks=[0, 0.2, 0.4, 0.6, 0.8, 1], fraction=0.046, pad=0.04)

    # This is how a correct labeling is done
    # the ticks are set from left (0) to right (len(matrix))
    # Somehow a mapping of the real TAD_start from genome position (chrX:123-456) to where this
    # start pos in matrix is needs to be computed.
    # ax.set_xticks([0, TAD_start, TAD_end, len(matrix)])
    # xticklabels = [None] * 4
    # xticklabels[0] = relabelTicks((int(referencePoint[1]) - region_start) * (-1))
    # xticklabels[1] = referencePoint[0] + ":" + relabelTicks(int(referencePoint[1]))
    # xticklabels[2] = relabelTicks(region_end - int(referencePoint[1]))
    # xticklabels[3] = relabelTicks(region_end - int(referencePoint[1]))

    # plt.xticks([0, 14, 29, 44], ['-' + args.regionLength, 'TAD start', 'TAD end', args.regionLength])
    # plt.yticks([0, 14, 29, 44], ['-' + args.regionLength, 'TAD start', 'TAD end', args.regionLength])
    plt.tight_layout()
    plt.savefig(args.outFileName, dpi=args.dpi)


if __name__ == "__main__":
    main()
