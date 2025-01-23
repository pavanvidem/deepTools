import deeptools.plotEnrichment
import deeptools.utilities
import os.path
from os import unlink
from matplotlib.testing.compare import compare_images

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"


def test_plotEnrichment():
    """
    Test minimal command line args for plotEnrichment
    """
    BAM = ROOT + "bowtie2_test1.bam"
    BED = ROOT + "multiBamSummary_regions.bed"
    outpng = '/tmp/test_plot.png'
    outcounts = '/tmp/test_counts.txt'

    args = "--bamfiles {0} {0} --BED {1} {1} --minMappingQuality 0 --labels bowtie2_test1.bam bowtie2_test1.bam --regionLabels up down --plotWidth 20 --plotHeight 20 --numPlotsPerRow 4 --alpha 0.9 --plotFileFormat png --plotFile {2} --outRawCounts {3}".format(BAM, BED, outpng, outcounts).split()
    deeptools.plotEnrichment.main(args)

    res = compare_images(ROOT + 'plotEnrichment_output.png', outpng, 50)


    assert res is None, res
    unlink(outpng)

    _foo = open(outcounts, 'r')
    resp = _foo.readlines()
    _foo.close()

    expected = ['file\tfeatureType\tpercent\tfeatureReadCount\ttotalReadCount\n',
                'bowtie2_test1.bam\tup\t100.00\t47\t47\n',
                'bowtie2_test1.bam\tdown\t100.00\t47\t47\n',
                'bowtie2_test1.bam\tup\t100.00\t47\t47\n',
                'bowtie2_test1.bam\tdown\t100.00\t47\t47\n']

    assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
    unlink(outcounts)
