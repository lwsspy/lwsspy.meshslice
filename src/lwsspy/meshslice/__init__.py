import sys

# Global imports when running the modules for testing
if "-m" not in sys.argv:

    # -------- CONSTANTS ------------------------------------------------------
    from .utils.constants import EARTH_RADIUS_M  # noqa
    from .utils.constants import EARTH_RADIUS_KM  # noqa
    from .utils.constants import EARTH_CIRCUM_M  # noqa
    from .utils.constants import EARTH_CIRCUM_KM  # noqa
    from .utils.constants import DEG2KM  # noqa
    from .utils.constants import KM2DEG  # noqa
    from .utils.constants import DEG2M  # noqa
    from .utils.constants import M2DEG  # noqa

    # Import the functions so that they are accesible from the main module
    from .read_mesh import read_mesh  # noqa
    from .meshslice_sph import MeshPlotSph  # noqa
    # from .meshslice_cart import MeshPlotCart
