import vtkmodules


def generate_sphere(radius, origin):
    """Return a _vtk.vtkSphere."""
    sphere = vtkmodules.vtkCommonDataModel.vtkSphere()
    sphere.SetRadius(radius)
    sphere.SetCenter(origin)
    return sphere
