import numpy as np


def cart2sph(x: float or np.ndarray or list,
             y: float or np.ndarray or list,
             z: float or np.ndarray or list) \
    -> (float or np.ndarray or list,
        float or np.ndarray or list,
        float or np.ndarray or list):
    """Computes Cartesian coordinates from geographical coordinates.

    Parameters
    ----------
    x : float or numpy.ndarray or list
        Radius
    y : float or numpy.ndarray or list
        Latitude (-90, 90)
    z : float or numpy.ndarray or list
        Longitude (-180, 180)

    Returns
    -------
    float or np.ndarray or list, float or np.ndarray or list, float or np.ndarray or list
        (r, theta, phi)
    """

    if type(x) is list:
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)

    # Compute Transformation
    r = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arctan2(y, x)

    # Catch phis below 0
    phi = np.where(phi < 0, phi + 2*np.pi, phi)

    # Catch corner case of latitude computation
    theta = np.where((r == 0), 0, np.arccos(z/r))

    if type(r) is list:
        r = r.tolist()
        theta = theta.tolist()
        phi = phi.tolist()

    return r, theta, phi
