from cartopy.crs import PlateCarree, Mollweide
import matplotlib.pyplot as plt


def map_axes(proj: str = "moll", central_longitude=0.0) -> plt.Axes:
    """Creates matplotlib axes with map projection taken from cartopy.

    Parameters
    ----------
    proj: str, optional
        shortname for mapprojection
        'moll', 'carr', by default "moll"

    Returns
    -------
    plt.Axes
        Matplotlib axes with projection

    Raises
    ------
    ValueError
        If non supported shortname for axes is given

    Notes
    -----

    :Author:
        Lucas Sawade (lsawade@princeton.edu)

    :Last Modified:
        2021.01.13 20.30

    Examples
    --------

    >>> from lwsspy import map_axes
    >>> map_axes()

    """

    # Check whether name is supported.
    if proj not in ['moll', 'carr']:
        raise ValueError(f"Either 'moll' for mollweide "
                         f"or 'carr' for PlateCarree.\n'{proj}'"
                         f"is not supported.")

    if proj == 'moll':
        projection = Mollweide(central_longitude=central_longitude)
    elif proj == 'carr':
        projection = PlateCarree(central_longitude=central_longitude)

    ax = plt.axes(projection=projection)

    return ax
