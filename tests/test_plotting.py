import os
import matplotlib.pyplot as plt
from meshslice.utils.map_axes import map_axes
from meshslice.utils.plot_map import plot_map
from meshslice import MeshPlotSph, read_mesh


DATADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')


def test_plotting():

    plt.figure()
    map_axes()
    plot_map()


def test_plot_mesh():

    # Testfile
    testfile = os.path.join(DATADIR, "test.vtu")

    # Read Mesh
    M = read_mesh(testfile)

    # Plot
    MeshPlotSph(M, test=True)
