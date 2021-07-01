import numpy as np


def sph2cart(r: float or np.ndarray or list,
             theta: float or np.ndarray or list,
             phi: float or np.ndarray or list) \
    -> (float or np.ndarray or list,
        float or np.ndarray or list,
        float or np.ndarray or list):
    """Computes Cartesian coordinates from geographical coordinates.

    Parameters
    ----------
    r : float or numpy.ndarray or list
        Radius
    theta : float or numpy.ndarray or list
        Theta (0, pi)
    phi : float or numpy.ndarray or list
        phi (0, 2*pi)

    Returns
    -------
    float or np.ndarray or list, float or np.ndarray or list, float or np.ndarray or list
        (x, y, z)
    """

    if type(r) is list:
        r = np.array(r)
        theta = np.array(theta)
        phi = np.array(phi)

    # Compute Transformation
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    if type(r) is list:
        x = x.tolist()
        y = y.tolist()
        z = z.tolist()

    return x, y, z
