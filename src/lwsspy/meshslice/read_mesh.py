
import pyvista as pv


def read_mesh(filename: str):

    mesh = mesh = pv.read(filename)

    return mesh
