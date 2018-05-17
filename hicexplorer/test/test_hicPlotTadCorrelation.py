from tempfile import NamedTemporaryFile

import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
import os.path
import pytest
from psutil import virtual_memory
mem = virtual_memory()
memory = mem.total / 2**30
import hicexplorer.hicPlotTadCorrelation
tolerance = 60  # default matplotlib pixed difference tolerance
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"

# memory in GB the test computer needs to have to run the test case
LOW_MEMORY = 2
MID_MEMORY = 7
HIGH_MEMORY = 200

REMOVE_OUTPUT = True


@pytest.mark.skipif(LOW_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_hicPlotMatrix_region_region2_log1p_clearMaskedBins_and_bigwig():

    outfile = NamedTemporaryFile(suffix='.png', prefix='test_hicPlotTadCorrelation', delete=False)

    args = "--deeptools_matrix {} -ts {} --outFileName {}".format(ROOT + 'hicPlotTadCorrelation/deeptoolMatrix.mat.gz',
                                                                  'chrX:1-2', outfile.name).split()
    test_image_path = ROOT + "hicPlotTadCorrelation" + '/plot.png'

    hicexplorer.hicPlotTadCorrelation.main(args)
    res = compare_images(test_image_path, outfile.name, tolerance)
    assert res is None, res

    if REMOVE_OUTPUT:
        os.remove(outfile.name)
