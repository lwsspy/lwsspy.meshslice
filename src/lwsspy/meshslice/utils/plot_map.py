import matplotlib.pyplot as plt
import cartopy
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.crs import PlateCarree, Mollweide


def plot_map(fill=True, zorder=None, labelstopright: bool = True,
             labelsbottomleft: bool = True, borders: bool = True,
             rivers: bool = False, lakes: bool = False):
    """Plots map into existing axes.

    Parameters
    ----------
    fill : bool, optional
        fills the continents in light gray, by default True
    zorder : int, optional
        zorder of the map, by default -10
    projection : cartopy.crs.projection, optional
        projection to be used for the map.
    labelstopright : bool, optional
        flag to turn on or off the ticks
    labelsbottomleft : bool, optional
        flag to turn on or off the ticks
    borders : bool
        plot borders. Default True
    rivers : bool
        plot rivers. Default False
    lakes : bool 
        plot lakes. Default True

    Returns
    -------
    matplotlib.pyplot.Axes
        Axes in which the map was plotted

    Notes
    -----

    :Author:
        Lucas Sawade (lsawade@princeton.edu)

    :Last Modified:
        2020.01.07 16.30


    """

    ax = plt.gca()

    # Put lables all around
    if isinstance(ax.projection, cartopy.crs.PlateCarree):

        # Set xticks Should be automated, but I just don't know how rn
        ax.set_xticks([-180, -135, -90, -45, 0, 45,
                       90, 135, 180], crs=ax.projection)
        ax.set_yticks([-90, -45, 0,  45, 90], crs=ax.projection)

        # Set label formatter
        ax.xaxis.set_major_formatter(cartopy.mpl.ticker.LongitudeFormatter())
        ax.yaxis.set_major_formatter(cartopy.mpl.ticker.LatitudeFormatter())

        ax.tick_params(
            labeltop=labelstopright, labelright=labelstopright,
            labelbottom=labelsbottomleft, labelleft=labelsbottomleft)

    # Set gridlines
    # gl = ax.gridlines(draw_labels=False, linewidth=1, color='lightgray',
    #
    #     alpha=0.5, linestyle='-', zorder=-1.5)

    # Add land
    if fill:
        ax.add_feature(cartopy.feature.LAND, zorder=zorder, edgecolor='black',
                       linewidth=0.5, facecolor=(0.9, 0.9, 0.9))
    else:
        ax.add_feature(cartopy.feature.LAND, zorder=zorder, edgecolor='black',
                       linewidth=0.5, facecolor=(0, 0, 0, 0))

    if borders:
        ax.add_feature(cartopy.feature.BORDERS,
                       zorder=None if zorder is None else zorder + 1,
                       facecolor='none', edgecolor=(0.5, 0.5, 0.5),
                       linewidth=0.25)

    if rivers:
        ax.add_feature(cartopy.feature.RIVERS, zorder=zorder,
                       edgecolor=(0.3, 0.3, 0.7),)
        #    edgecolor=(0.5, 0.5, 0.7) )

    if lakes:
        ax.add_feature(cartopy.feature.LAKES,
                       zorder=None if zorder is None else zorder + 1,
                       edgecolor='black', linewidth=0.5,
                       facecolor=(1.0, 1.0, 1.0))
    return ax
