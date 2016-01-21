#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import numpy as np
from matplotlib import use as mplt_use
mplt_use('Agg')
import matplotlib.pyplot as plt

from deeptools.correlation import Correlation
from deeptools._version import __version__


def parse_arguments(args=None):
    basic_args = plot_correlation_args()
    heatmap_parser = heatmap_options()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Tool for the analysis and visualization of sample correlations based on the output of multiBamCoverage or
multiBigwigSummary. Pearson or Spearman methods are available to compute correlation
coefficients. Results can be saved into a heatmap image or as multiple
scatter plots depicting the pairwise correlations.
Optionally, the values can be saved as tables, too.


detailed help:

  plotCorrelation -h

""",
        epilog='example usages:\n'
               'plotCorrelation -in results_file --whatToPlot heatmap --corMethod pearson -o heatmap.png\n\n'
               ' \n\n',
        parents=[basic_args, heatmap_parser])

    return parser


def plot_correlation_args():
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')

    # define the arguments
    required.add_argument('--corData', '-in',
                          metavar='FILE',
                          help='Compressed matrix of values generated by multiBigwigSummary or multiBamCoverage',
                          required=True)

    required.add_argument('--plotFile', '-o',
                          help='File to save the heatmap to. The file extension determines the format, '
                          'so heatmap.pdf will save the heatmap in PDF format. '
                          'The available formats are: .png, '
                          '.eps, .pdf and .svg.',
                          type=argparse.FileType('w'),
                          metavar='FILE',
                          required=True)

    required.add_argument('--corMethod', '-c',
                          help="Correlation method.",
                          choices=['spearman', 'pearson'],
                          required=True)

    required.add_argument('--whatToPlot', '-p',
                          help="Choose between a heatmap or pairwise scatter plots",
                          choices=['heatmap', 'scatterplot'],
                          required=True)

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--skipZeros',
                          help='By setting this option, genomic regions '
                          'that have zero or missing (nan) values in all samples '
                          'are excluded.',
                          action='store_true',
                          required=False)

    optional.add_argument('--labels', '-l',
                          metavar='sample1 sample2',
                          help='User defined labels instead of default labels from '
                          'file names. '
                          'Multiple labels have to be separated by spaces, e.g. '
                          '--labels sample1 sample2 sample3',
                          nargs='+')

    optional.add_argument('--plotTitle', '-T',
                          help='Title of the plot, to be printed on top of '
                          'the generated image. Leave blank for no title.',
                          default='')

    optional.add_argument('--plotFileFormat',
                          metavar='FILETYPE',
                          help='Image format type. If given, this option '
                          'overrides the image format based on the plotFile '
                          'ending. The available options are: png, '
                          'eps, pdf and svg.',
                          choices=['png', 'pdf', 'svg', 'eps'])

    optional.add_argument(
        '--removeOutliers',
        help='If set, bins with very large counts are removed. '
             'Bins with abnormally high reads counts artificially increase '
             'pearson correlation; that\'s why, by default, multiBamCoverage tries '
             'to remove outliers using the median absolute deviation (MAD) '
             'method applying a threshold of 200 to only consider extremely '
             'large deviations from the median. The ENCODE blacklist page '
             '(https://sites.google.com/site/anshulkundaje/projects/blacklists) '
             'contains useful information about regions with unusually high counts.',
        action='store_true')

    optional.add_argument('--version', action='version',
                          version='%(prog)s {}'.format(__version__))

    group = parser.add_argument_group('Output optional options')

    group.add_argument('--outFileCorMatrix',
                       help='Save matrix with pairwise correlation values to a tab-separated file.',
                       metavar='FILE',
                       type=argparse.FileType('w'))

    return parser


def heatmap_options():
    """
    Options for generating the correlation heat map
    """
    parser = argparse.ArgumentParser(add_help=False)
    heatmap = parser.add_argument_group('Heatmap options')

    heatmap.add_argument('--zMin', '-min',
                         default=None,
                         help='Minimum value for the heatmap intensities. '
                              'If not specified, the value is set automatically',
                         type=float)

    heatmap.add_argument('--zMax', '-max',
                         default=None,
                         help='Maximum value for the heatmap intensities.'
                              'If not specified, the value is set automatically',
                         type=float)

    heatmap.add_argument(
        '--colorMap', default='jet',
        metavar='',
        help='Color map to use for the heatmap. Available values can be '
             'seen here: '
             'http://matplotlib.org/examples/color/colormaps_reference.html')

    heatmap.add_argument('--plotNumbers',
                         help='If set, then the correlation number is plotted '
                         'on top of the heatmap. This option is only valid when plotting a heatmap.',
                         action='store_true',
                         required=False)

    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)

    corr = Correlation(args.corData,
                       args.corMethod,
                       labels=args.labels,
                       remove_outliers=args.removeOutliers,
                       log1p=None,
                       skip_zeros=args.skipZeros)

    if args.corMethod == 'pearson':
        # test if there are outliers and write a message recommending the removal
        if len(corr.get_outlier_indices(np.asarray(corr.matrix).flatten())) > 0:
            if args.removeOutliers:
                            sys.stderr.write("\nOutliers were detected in the data. They "
                                             "will be removed to avoid bias "
                                             "in the pearson correlation.\n")

            else:
                sys.stderr.write("\nOutliers were detected in the data. Consider "
                                 "using the --removeOutliers parameter to avoid a bias "
                                 "in the pearson correlation.\n")

    if args.colorMap:
        try:
            plt.get_cmap(args.colorMap)
        except ValueError as error:
            sys.stderr.write(
                "A problem was found. Message: {}\n".format(error))
            exit()

    if args.outFileCorMatrix:
        corr.save_corr_matrix(args.outFileCorMatrix)

    args.plotFile.close()
    if args.whatToPlot == 'scatterplot':
        corr.plot_scatter(args.plotFile.name,
                          plot_title=args.plotTitle,
                          image_format=args.plotFileFormat)
    else:
        corr.plot_correlation(args.plotFile.name,
                              vmax=args.zMax,
                              vmin=args.zMin,
                              colormap=args.colorMap,
                              plot_title=args.plotTitle,
                              image_format=args.plotFileFormat,
                              plot_numbers=args.plotNumbers)
