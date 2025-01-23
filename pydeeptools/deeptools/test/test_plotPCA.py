import deeptools.plotPCA
import deeptools.utilities
import os.path
from os import unlink
from matplotlib.testing.compare import compare_images
import numpy as np
import numpy.testing as nt

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"

def test_plotPCA():
    """
    Test minimal command line args for plotPCA
    """
    COR = ROOT + "multiBamSummary_result2b.npz"
    outpng = '/tmp/test_plot.png'
    outdata = '/tmp/test_plot_data.txt'

    args = "--corData {} --plotTitle Test-Plot --plotFile {} --outFileNameData {}".format(COR,outpng, outdata).split()
    deeptools.plotPCA.main(args)

    res = compare_images(ROOT + 'plotPCA_result.png', outpng, 50)
    assert res is None, res
    unlink(outpng)

    resp = np.loadtxt(outdata, skiprows=2, delimiter='\t')

    expected = np.array([[1,-0.7071067811865476,-0.7071067811865475,4.0],
                         [2,-0.7071067811865475,0.7071067811865476,1.2325951644078315e-28]])

    nt.assert_allclose(resp, expected, atol=1e-10)
    unlink(outdata)
