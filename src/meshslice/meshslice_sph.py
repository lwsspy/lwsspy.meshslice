
# External
import os
from typing import Optional, List
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation, LinearTriInterpolator
from copy import deepcopy
import pyvista as pv
from pyvista import _vtk
from pyvista.utilities.helpers import generate_plane


# Internal imports
from .utils.rotation_matrix import rotation_matrix
from .utils.visibility_callback import SetVisibilityCallback
from .utils.vtksphere import generate_sphere
from .utils.updaterc import updaterc
from .utils.geo2cart import geo2cart
from .utils.cart2geo import cart2geo
from .utils.pol2cart import pol2cart
from .utils.cart2pol import cart2pol
from .utils.map_axes import map_axes
from .utils.plot_map import plot_map
from .utils.SphericalNN import SphericalNN
from . import EARTH_RADIUS_KM
from . import DEG2KM


class MeshPlotSph():

    def __init__(self, mesh,
                 lat: Optional[float] = None,
                 lon: Optional[float] = None,
                 rotangle: float = 0.0,
                 depth: Optional[float] = None,
                 rotlat: Optional[float] = None,
                 rotlon: Optional[float] = None,
                 get_mean: bool = True,
                 opacity: bool = False,
                 meshname: str = 'RF',
                 cmapname: str = "seismic",
                 cmapsym: bool = True,
                 clim: Optional[List[float]] = None,
                 debug: bool = True,
                 outputdir: str = '.',
                 use_ipyvtk: bool = False,
                 figures_only: bool = False,
                 test: bool = False,
                 fmt: str = "pdf"):

        # Init state
        self.not_init_state = False
        self.debug = debug
        self.fmt = fmt

        # Copy to not modify original mesh when rotating
        self.mesh = mesh.copy(deep=True)

        # Function to get activae scalar
        self.meshname = meshname
        self.mesh.set_active_scalars(self.meshname)

        # To rotate data set
        if rotlon is not None:
            pv.core.common.axis_rotation(
                self.mesh.points, rotlon, inplace=True, deg=True, axis='z')
        if rotlat is not None:
            pv.core.common.axis_rotation(
                self.mesh.points, rotlat, inplace=True, deg=True, axis='y')

        # Get bounds
        self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax = \
            self.mesh.bounds

        # Get min/max radius, depth
        self.r = np.sqrt(np.sum(deepcopy(self.mesh.points)**2, axis=1))
        self.rmin = np.min(self.r)
        self.rmax = np.max(self.r)
        self.dmin = self.rmax - self.rmax
        self.dmax = self.rmax - self.rmin

        # Get the depth
        if depth is not None:
            self.depth = depth
        else:
            self.depth = (self.dmin+self.dmax)/2

        # Get colorbounds
        self.minM, self.maxM = self.mesh.get_data_range()

        # Get initial position depending on inputt position or point mean
        if lat is not None and lon is not None:
            # Get location
            self.latitude, self.longitude = lat, lon

            # Convert to vector
            self.initpos = geo2cart(1.0, self.latitude, self.longitude)

        else:
            # Get mean vector
            self.initpos = np.mean(self.mesh.points, axis=0)

            # Get convert to geo location.
            _, self.latitude, self.longitude = cart2geo(*self.initpos)

        # Set center of the slice to be the init
        self.center = self.initpos

        # Only set so that it doesn't have to be redefined
        self.normalx = np.array([1, 0, 0])
        self.normaly = np.array([0, 1, 0])
        self.normalz = np.array([0, 0, 1])

        # Starting parameters
        self.rotmat = np.eye(3)
        self.rotmat_lat = np.eye(3)
        self.rotmat_lon = np.eye(3)
        self.rotmat_slice = np.eye(3)
        self.rotangle = rotangle
        self.init_rotangle = rotangle
        self.radius = (self.rmin+self.rmax)/2

        # Empty slice parameters to be populated
        self.Yslice = dict(
            slc=None,
            actor2D=None,
            actor3D=None,
            normal=[0, 1, 0],
            origin=[0, 0, 0])
        self.Zslice = dict(
            slc=None,
            actor2D=None,
            actor3D=None,
            normal=[0, 1, 0],
            origin=[0, 0, 0])
        self.Rslice = dict(
            slc=None,
            actor2D=None,
            actor3D=None,
            radius=self.radius,
            center=[0, 0, 0])

        self.cmapname = cmapname
        self.cmap = plt.cm.get_cmap(self.cmapname, 256)

        # Set up colorbar
        if clim is not None:
            self.clim = clim
            if -self.clim[0] == self.clim[1]:
                self.cmapsym = True
            else:
                self.cmapsym = False
        else:
            self.clim = clim
            self.cmapsym = cmapsym

        # Find scalar bar max and min
        self.cabsmax = np.max(np.abs([self.minM, self.maxM]))
        self.minM, self.maxM = [-self.cabsmax, self.cabsmax]

        # Find colorbar bounds depending on symmetry
        if self.cmapsym:

            # If clim is defined get it
            if self.clim is not None:
                self.initclim = deepcopy(self.clim)
            else:
                self.clim = [-0.5*self.cabsmax, 0.5*self.cabsmax]
                self.initclim = [-0.5*self.cabsmax, 0.5*self.cabsmax]

            # Set slider ranges
            self.cminrange = [self.minM, 0.0]
            self.cmaxrange = [0, self.maxM]

        else:
            # If clim is defined get it
            if self.clim is not None:
                self.initclim = deepcopy(self.clim)
            else:
                self.clim = [self.minM, self.maxM]
                self.initclim = [self.minM, self.maxM]

            # Set slider ranges
            self.cminrange = [self.minM, self.maxM]
            self.cmaxrange = [self.minM, self.maxM]

        self.scalar_bar_args = dict(
            width=0.4, height=0.05, position_x=0.3, position_y=0.05)

        # Check for illumination and plot it if possible.
        if self.mesh['illumination'] is not None and opacity is True:

            raise ValueError("Not working/not correcly implemented")

            self.init_illum = 50
            self.opacity = 'opacity'
            self.mesh['opacity'] = np.where(
                self.mesh['illumination'] >= self.init_illum, 1.0,
                self.mesh['illumination']/self.init_illum)
            self.mesh.set_active_scalars(self.meshname)

            # Opacity
            self.mesh['opacity'] = np.ones_like(self.mesh['illumination'])
            self.opacity_slider = self.p.add_slider_widget(
                self.opacity_function, value=self.init_illum,
                title='Illumination', rng=[0.0, 500.0],
                pointa=(.525, .75), pointb=(.975, .75))

        else:
            self.opacity: float = 1.0

        # Initiate Plotter
        pv.rcParams['multi_rendering_splitting_position'] = 0.7
        self.p = pv.Plotter(shape='1|4', window_size=(1600, 900))

        # Start first subplot
        self.p.subplot(0)

        # Plot background volume
        self.volume = self.p.add_mesh(
            self.mesh,
            opacity=0.1,
            stitle=self.meshname,
            clim=self.clim,
            cmap=self.cmap,
            show_scalar_bar=False)

        # Plot colorbar
        self.p.add_scalar_bar(
            title=self.meshname, **self.scalar_bar_args)

        # Add orientation widget
        self.p.add_orientation_widget(self.mesh)

        # Get slices
        self.get_slices()

        # Rotation
        self.p.subplot(1)
        self.rotate_slider = self.p.add_slider_widget(
            self.rotate_slice,
            value=self.rotangle,
            title='Rot',
            rng=[0, 90],
            pointa=(.025, .1),
            pointb=(.975, .1))  # , event_type='always')  # This makes the slider an activae listener

        # Longitude value
        self.p.subplot(1)
        self.rotlon_slider = self.p.add_slider_widget(
            self.rotate_lon,
            value=self.longitude,
            title='Lon',
            rng=[-180, 180],
            pointa=(.025, .25),
            pointb=(.975, .25))

        # Longitude value
        self.p.subplot(1)
        self.rotlat_slider = self.p.add_slider_widget(
            self.rotate_lat,
            value=self.latitude,
            title='Lat',
            rng=[-90, 90],
            pointa=(.025, .4),
            pointb=(.975, .4))

        # Depth slice
        self.p.subplot(1)
        self.r_slider = self.p.add_slider_widget(
            self.update_r,
            value=self.depth,
            title='Z',
            rng=[self.dmin, self.dmax],
            pointa=(.025, .55),
            pointb=(.975, .55),
            event_type='always')

        # Add slice actors
        self.add_slice_actors()

        # Colorbar/Scalar Widget
        self.p.subplot(1)

        # Update slider rcParams
        pv.rcParams['slider_style']['modern'].update(
            dict(slider_width=0.01, cap_width=0.01, tube_width=0.001))

        # Update minimum cmap value slider
        self.p.subplot(1)
        self.cmin_slider = self.p.add_slider_widget(
            self.cmin_callback,
            value=self.clim[0],
            title='cmin',
            rng=self.cminrange,
            pointa=(.025, .9),
            pointb=(.475, .9),
            event_type='always')

        # Update maximum cmap value slider
        self.p.subplot(1)
        self.cmax_slider = self.p.add_slider_widget(
            self.cmax_callback,
            value=self.clim[1],
            title='cmax',
            rng=self.cmaxrange,
            pointa=(.525, .9),
            pointb=(.975, .9),
            event_type='always')

        # Plot bounds
        self.p.subplot(0)
        bounds = self.p.show_bounds(bounds=mesh.bounds, location='back')

        # Define button and text positions
        font_size = 14
        x_offset = 80
        y_offset = 50
        size = 40
        checkbox_pos = [30, 30]
        text_pos = [x_offset, 30]

        # Volume visibiliy
        self.p.subplot(2)
        self.p.add_text("Volume", position=text_pos, font_size=font_size)
        self.p.add_checkbox_button_widget(
            SetVisibilityCallback(self.volume),
            value=True,
            size=size,
            position=checkbox_pos,
            color_on='white',
            color_off='grey',
            background_color='grey')
        checkbox_pos[1] += y_offset
        text_pos[1] += y_offset

        # Turn on/off bounds
        self.p.subplot(2)
        self.p.add_text("Bounds", position=text_pos, font_size=font_size)
        self.p.add_checkbox_button_widget(
            SetVisibilityCallback(bounds),
            value=True,
            size=size,
            position=checkbox_pos,
            color_on='white',
            color_off='grey',
            background_color='grey')
        for i in range(2):
            checkbox_pos[1] += y_offset
            text_pos[1] += y_offset

        # Colormap symmetric on/off button.
        self.p.subplot(2)
        self.p.add_text("Symmetric cmap", position=text_pos,
                        font_size=font_size)
        self.p.add_checkbox_button_widget(
            self.csym_callback,
            value=self.cmapsym,
            size=size,
            position=checkbox_pos,
            color_on='white',
            color_off='grey',
            background_color='grey')

        for i in range(2):
            checkbox_pos[1] += y_offset
            text_pos[1] += y_offset

        # Export button for Y Slice
        self.p.subplot(2)
        self.p.add_text("Export Y-Slice", position=text_pos,
                        font_size=font_size)
        self.p.add_checkbox_button_widget(
            self.plot_y_slice,
            value=True,
            size=size,
            position=checkbox_pos,
            color_on='white',
            color_off='white',
            background_color='grey')
        checkbox_pos[1] += y_offset
        text_pos[1] += y_offset

        # Export button for Y Slice
        self.p.subplot(2)
        self.p.add_text("Export Z-Slice", position=text_pos,
                        font_size=font_size)
        self.p.add_checkbox_button_widget(
            self.plot_z_slice,
            value=True,
            size=size,
            position=checkbox_pos,
            color_on='white',
            color_off='white',
            background_color='grey')
        checkbox_pos[1] += y_offset
        text_pos[1] += y_offset

        # Export button for depth slice
        self.p.subplot(2)
        self.p.add_text("Export Depth Slice", position=text_pos,
                        font_size=font_size)
        self.p.add_checkbox_button_widget(
            self.depth_slice_to_map,
            value=True,
            size=size,
            position=checkbox_pos,
            color_on='white',
            color_off='white',
            background_color='grey')

        # Init state
        self.outputdir = outputdir

        # Finalize Window creation
        self.not_init_state = True

        # Necessary for first update of all plots
        self.Rslice['slc'].shallow_copy(self.Rslice['alg'].GetOutput())
        self.Yslice['slc'].shallow_copy(self.Yslice['alg'].GetOutput())
        self.Zslice['slc'].shallow_copy(self.Zslice['alg'].GetOutput())

        # Set view (Updates subplots)
        self.p.subplot(3)
        self.p.view_vector(-self.Yslice['normal'], viewup=self.normalup)
        self.p.subplot(4)
        self.p.view_vector(-self.Zslice['normal'], viewup=self.normalup)

        # Do general update
        self.p.update()

        if test:
            self.p.close()
            return

        # Finally show the damn thing
        self.figures_only = figures_only
        if self.figures_only:
            self.depth_slice_to_map(1.0)
            self.plot_y_slice(1.0)
            self.plot_z_slice(1.0)
            self.p.close()
            plt.close('all')
        else:
            self.p.show(return_viewer=True, use_ipyvtk=use_ipyvtk)

    def get_slices(self):

        if not hasattr(self, "plane_sliced_meshes"):
            self.p.plane_sliced_meshes = []

        self.Yslice['alg'] = _vtk.vtkCutter()  # Construct the cutter object
        # Use the grid as the data we desire to cut
        self.Yslice['alg'].SetInputDataObject(self.mesh)
        self.Yslice['alg'].GenerateTrianglesOff()
        self.Yslice['slc'] = pv.wrap(self.Yslice['alg'].GetOutput())
        self.p.plane_sliced_meshes.append(self.Yslice['slc'])

        self.Zslice['alg'] = _vtk.vtkCutter()  # Construct the cutter object
        # Use the grid as the data we desire to cut
        self.Zslice['alg'].SetInputDataObject(self.mesh)
        self.Zslice['alg'].GenerateTrianglesOff()
        self.Zslice['slc'] = pv.wrap(self.Zslice['alg'].GetOutput())
        self.p.plane_sliced_meshes.append(self.Zslice['slc'])

        self.Rslice['alg'] = _vtk.vtkCutter()  # Construct the cutter object
        # Use the grid as the data we desire to cut
        self.Rslice['alg'].SetInputDataObject(self.mesh)
        self.Rslice['alg'].GenerateTrianglesOff()
        self.Rslice['slc'] = pv.wrap(self.Rslice['alg'].GetOutput())
        self.p.plane_sliced_meshes.append(self.Rslice['slc'])

    def add_slice_actors(self):

        self.p.subplot(0)
        self.Yslice['actor3D'] = self.p.add_mesh(self.Yslice['slc'],
                                                 stitle=self.meshname,
                                                 clim=self.clim,
                                                 cmap=self.cmap,
                                                 opacity=self.opacity)
        self.p.subplot(0)
        self.Zslice['actor3D'] = self.p.add_mesh(self.Zslice['slc'],
                                                 stitle=self.meshname,
                                                 clim=self.clim,
                                                 cmap=self.cmap,
                                                 opacity=self.opacity)

        self.p.subplot(0)
        self.Rslice['actor3D'] = self.p.add_mesh(self.Rslice['slc'],
                                                 stitle=self.meshname,
                                                 clim=self.clim,
                                                 cmap=self.cmap,
                                                 opacity=self.opacity)

        self.p.subplot(3)
        self.p.add_text("Y", position='upper_left', font_size=16)
        self.Yslice['actor2D'] = self.p.add_mesh(self.Yslice['slc'],
                                                 stitle=self.meshname,
                                                 clim=self.clim,
                                                 cmap=self.cmap,
                                                 opacity=self.opacity)
        self.p.subplot(4)
        self.p.add_text("Z", position='upper_left', font_size=16)
        self.Zslice['actor2D'] = self.p.add_mesh(self.Zslice['slc'],
                                                 stitle=self.meshname,
                                                 clim=self.clim,
                                                 cmap=self.cmap,
                                                 opacity=self.opacity)

    def r_slice_update(self):

        # create the plane for clipping
        sphere = generate_sphere(self.Rslice['radius'], self.Rslice['center'])

        # the cutter to use the plane we made
        self.Rslice['alg'].SetCutFunction(sphere)
        self.Rslice['alg'].Update()  # Perform the Cut
        self.Rslice['slc'].shallow_copy(self.Rslice['alg'].GetOutput())

    def y_slice_update(self):
        # Rotation matrix of the normal
        if self.debug:
            print("  Getting Normal up")
        self.normalup = self.rotmat @ self.normalx

        if self.debug:
            print("  Getting Slice normal")
        self.Yslice['normal'] = self.rotmat @ self.normaly

        # create the plane for clipping
        plane = generate_plane(self.Yslice['normal'], self.Yslice['origin'])

        # Update cut function
        self.Yslice['alg'].SetCutFunction(plane)
        self.Yslice['alg'].Update()  # Perform the Cut
        self.Yslice['slc'].shallow_copy(self.Yslice['alg'].GetOutput())

        # Update view in 2D plot
        self.p.subplot(3)
        self.p.view_vector(-self.Yslice['normal'], viewup=self.normalup)

    def z_slice_update(self):
        # Rotation matrix of the normal
        if self.debug:
            print("  Getting Normal up")
        self.normalup = self.rotmat @ self.normalx

        if self.debug:
            print("  Getting Slice normal")
        self.Zslice['normal'] = self.rotmat @ self.normalz

        # create the plane for clipping
        plane = generate_plane(self.Zslice['normal'], self.Zslice['origin'])

        # Update cut function
        self.Zslice['alg'].SetCutFunction(plane)
        self.Zslice['alg'].Update()  # Perform the Cut
        self.Zslice['slc'].shallow_copy(self.Zslice['alg'].GetOutput())

        # Update view in 2D plot
        self.p.subplot(4)
        self.p.view_vector(-self.Zslice['normal'], viewup=self.normalup)

    def opacity_function(self, val):

        raise ValueError("Not working/not correcly implemented")

        if val == 0.0:
            self.mesh['opacity'] = np.ones_like(self.mesh['illumination'])
        else:
            self.mesh['opacity'] = np.where(
                self.mesh['illumination'] >= val, 1.0,
                self.mesh['illumination']/val)
            self.mesh['opacity'] = np.ones_like(self.mesh['illumination'])

        # Update
        self.mesh.set_active_scalars(self.meshname)

    def csym_callback(self, val):

        # Check whether cmap symmetrical or not
        self.cmapsym = val

        # if cmap symmetric find new values and update the colormap
        if self.cmapsym:
            averageclim = np.sum(np.abs(self.clim))/2
            self.clim = [-averageclim, averageclim]

        # Reset colorslider bounds
        self.set_colorslider_bounds()
        self.set_colorslider_values()
        self.p.update_scalar_bar_range(self.clim, name=self.meshname)
        self.p.render()

    def cmin_callback(self, val):

        # Get clims depending on cmap symmetry
        if self.cmapsym:
            self.clim = [val, -val]
        else:
            self.clim = [val, self.clim[1]]

        # Update colormap
        self.p.update_scalar_bar_range(self.clim, name=self.meshname)

        # If not init state, change the value of the opoosit slider
        if self.not_init_state:

            # If symmetric update the slider location on the opposite side
            if self.cmapsym:
                self.cmax_slider.GetRepresentation(
                ).SetValue(-val)

            # Don't do anything if val is equal to the maximum of the
            # maximum slider
            else:
                if val > self.clim[1]:
                    pass
                else:
                    self.cmax_slider.GetRepresentation(
                    ).SetMinimumValue(self.clim[0])

    def cmax_callback(self, val):

        # Get clims depending on cmap symmetry
        if self.cmapsym:
            self.clim = [-val, val]
        else:
            self.clim = [self.clim[0], val]

        # Update colormap
        self.p.update_scalar_bar_range(self.clim, name=self.meshname)

        # If not init state, change the value of the opoosit slider
        if self.not_init_state:

            # If symmetric update the slider location on the opposite side
            if self.cmapsym:
                self.cmin_slider.GetRepresentation(
                ).SetValue(-val)

            # If not make maximum of minimum slider the minimum
            else:

                # Don't do anything if val is equal to the minimum of the
                # minimum slider
                if val < self.clim[0]:
                    pass
                else:
                    self.cmin_slider.GetRepresentation(
                    ).SetMaximumValue(self.clim[1])

    def rotate_lat(self, latitude, update_mat_only=False):

        # Get latitude
        self.latitude = latitude

        # Compute the center
        self.center = geo2cart(1.0, self.latitude, self.longitude)

        # Compute new rotation matrix using the latitude
        self.rotmat_lat = rotation_matrix(
            np.array([0, 1, 0]), -self.latitude/180.0*np.pi)

        # Update slice rotation matrix, since center has changed
        self.rotate_slice(self.rotangle, update_mat_only=True)

        # Update
        if update_mat_only is False:
            self.update()

    def rotate_lon(self, longitude, update_mat_only=False):

        # Get longitude
        self.longitude = longitude

        # Compute center
        self.center = geo2cart(1.0, self.latitude, self.longitude)

        # Compute rotation matrix
        self.rotmat_lon = rotation_matrix(
            np.array([0, 0, 1]), self.longitude/180.0*np.pi)

        # Update slice rotation matrix, since center has changed
        self.rotate_slice(self.rotangle, update_mat_only=True)

        # Update the slices
        if update_mat_only is False:
            self.update()

    def rotate_slice(self, angle, update_mat_only=False):

        # Compute rotangle in radians
        if update_mat_only is False:
            self.rotangle = -angle/180.0*np.pi

        # Update slice rotation matrix
        self.rotmat_slice = rotation_matrix(self.center, self.rotangle)

        # Update slices
        if update_mat_only is False:
            self.update()

    def update_r(self, d):

        # Get radius from depth input
        self.depth = d
        self.Rslice['radius'] = self.rmax - d

        # Update tthe radial slice
        self.r_slice_update()

    def update(self):

        # Compute the total rotation matrix
        if self.debug:
            print("Computing total rotation matrix")
        self.rotmat = self.rotmat_slice @ self.rotmat_lon @ self.rotmat_lat

        # Update the Y slice
        if self.debug:
            print("Update Y slice")
        self.y_slice_update()

        # Update the Z slice
        if self.debug:
            print("Update Z slice")
        self.z_slice_update()

    def set_colorslider_values(self):
        # Set cmap sliders:
        self.cmax_slider.GetRepresentation().SetValue(self.clim[0])
        self.cmax_slider.GetRepresentation().SetValue(self.clim[1])

    def set_colorslider_bounds(self):
        # Set cmap sliders:
        if self.cmapsym:
            self.cmin_slider.GetRepresentation().SetMinimumValue(self.minM)
            self.cmin_slider.GetRepresentation().SetMaximumValue(0.0)
            self.cmax_slider.GetRepresentation().SetMinimumValue(0.0)
            self.cmax_slider.GetRepresentation().SetMaximumValue(self.maxM)
        else:
            self.cmin_slider.GetRepresentation().SetMinimumValue(self.minM)
            self.cmin_slider.GetRepresentation().SetMaximumValue(self.maxM)
            self.cmin_slider.GetRepresentation().SetMinimumValue(self.minM)
            self.cmax_slider.GetRepresentation().SetMaximumValue(self.maxM)

    def reset(self, val):

        #  Get starting values
        self.clim = self.initclim
        self.center = self.initpos
        _, self.latitude, self.longitude = cart2geo(*self.initpos)

        # Reset colorbar boundaries
        self.set_colorslider_bounds()
        self.set_colorslider_values()

        # Reset rotationsliders
        self.rotlat_slider.GetRepresentation().SetValue(self.latitude)
        self.rotlon_slider.GetRepresentation().SetValue(self.longitude)
        self.rotate_slider.GetRepresentation().SetValue(self.init_rotangle)

        # Update matrices
        self.rotate_lat(self.latitude, update_mat_only=True)
        self.rotate_lon(self.latitude, update_mat_only=True)
        self.rotate_slice(self.latitude, update_mat_only=True)

        # Update plot
        self.p.update_scalar_bar_range(self.clim, name=self.meshname)
        self.update()
        self.p.render()

    def depth_slice_to_map(self, val):

        # Depth slice to geo
        slctemp = self.Rslice["slc"].copy(deep=True)

        # Get the geographical points
        _, lat, lon = cart2geo(
            slctemp.points[:, 0], slctemp.points[:, 1], slctemp.points[:, 2])

        # Get bounds
        latmin, latmax = np.min(lat), np.max(lat)
        lonmin, lonmax = np.min(lon), np.max(lon)

        # Get data
        data = deepcopy(slctemp[self.meshname])

        # Set up interpolation
        snn = SphericalNN(lat, lon)
        res = 0.25
        llon, llat = np.meshgrid(
            np.arange(lonmin, lonmax+res, res),
            np.arange(latmin, latmax+res, res))

        # Get the interpolated data using weighted spherical nearest neighbour
        # interpolation
        d = snn.interp(data, llat, llon, no_weighting=False,
                       k=10, maximum_distance=res*2.0)

        # Create Map
        fig = plt.figure()
        ax = map_axes(proj='carr')
        plot_map(zorder=1, borders=False, fill=False)
        ax.set_xlim(lonmin, lonmax)
        ax.set_ylim(latmin, latmax)

        # Plot data onto map
        plt.imshow(d[::-1, :], extent=[lonmin, lonmax, latmin,
                                       latmax], cmap=self.cmapname,
                   zorder=-1, vmin=self.clim[0], vmax=self.clim[1])
        cbar = plt.colorbar(orientation='horizontal', aspect=40)
        label = 'Hitcounts' if self.meshname == 'illumination' else "Amplitude"
        cbar.set_label(label)
        plt.title(f"Z = {int(self.depth):d} km")
        if not self.figures_only:
            plt.show(block=False)
        plt.savefig(os.path.join(self.outputdir,
                                 f"depthslice{int(self.depth):d}_"
                                 f"{self.meshname}"
                                 f".{self.fmt}"))

    def plot_y_slice(self, val):

        self.slice_to_array(
            self.Yslice["slc"],
        )

    def plot_z_slice(self, val):

        if self.debug:
            print("Plotting Z slice")
        self.slice_to_array(
            self.Zslice["slc"],
            offset=90.0,
            name="Z"
        )

    def slice_to_array(self, slc, offset: float = 0.0, name: str = "Y"):
        """Converts a PolyData slice to a 2D NumPy array.

        It is crucial to have the true normal and origin of
        the slicing plane

        Parameters
        ----------
        slc : PolyData
            Slice as slice through a mesh.

        """

        # Get slice
        if self.debug:
            print("   Getting Slice")
        slctemp = slc.copy(deep=True)
        slctemp.points = (self.rotmat.T @ slctemp.points.T).T

        # Interpolate slices
        if self.debug:
            print("   Creating polydata")

        # Depending on the slice, column 1 or 2 will contain the points
        if offset == 0.0:
            points = np.vstack(
                (slctemp.points[:, 0],
                 slctemp.points[:, 2],
                 np.zeros_like(slctemp.points[:, 2]))).T
        else:
            points = np.vstack(
                (slctemp.points[:, 0],
                 slctemp.points[:, 1],
                 np.zeros_like(slctemp.points[:, 1]))).T

        pc = pv.PolyData(points)

        if self.debug:
            print("   Triangulate")
        mesh = pc.delaunay_2d(alpha=0.5*1.5*DEG2KM)
        mesh[self.meshname] = slctemp[self.meshname]

        # Get triangles from delauney triangulation to be
        # used for interpolation
        if self.debug:
            print("   Reshape Triangles")
        xy = np.array(mesh.points[:, 0:2])
        r, t = cart2pol(xy[:, 0], xy[:, 1])
        findlimt = t + 4 * np.pi
        mint = np.min(findlimt) - 4 * np.pi
        maxt = np.max(findlimt) - 4 * np.pi
        cells = mesh.faces.reshape(mesh.n_cells, 4)
        triangles = np.array(cells[:, 1:4])

        # Set maximum slice length (makes no sense otherwise)
        if mint < -11.25:
            mint = -11.25
        if maxt > 11.25:
            maxt = 11.25

        # Set up intterpolation values
        dmin = np.min(EARTH_RADIUS_KM - r)
        dmax = np.max(EARTH_RADIUS_KM - r)
        dsamp = np.linspace(dmin, dmax, 1000)
        tsamp = np.linspace(mint, maxt, 1000)
        tt, dd = np.meshgrid(tsamp, dsamp)

        # you can add keyword triangles here if you have the triangle array,
        # size [Ntri,3]
        if self.debug:
            print("   Triangulation 2")
        triObj = Triangulation(t, r, triangles=triangles)

        # linear interpolation
        fz = LinearTriInterpolator(triObj, mesh[self.meshname])
        Z = fz(tt, EARTH_RADIUS_KM - dd)

        # Get aspect of the figure
        aspect = ((maxt - mint)*180/np.pi*DEG2KM) / (dmax - dmin)
        height = 4.0
        width = height * aspect

        if self.debug:
            print("   Plotting")

        # Create figure
        plt.figure(figsize=(width, height))

        plt.imshow(
            Z,
            extent=[mint*180/np.pi*DEG2KM, maxt*180/np.pi*DEG2KM,
                    dmax, dmin],
            cmap=self.cmapname, vmin=self.clim[0], vmax=self.clim[1],
            rasterized=True)

        from_north = np.abs(self.rotangle*180/np.pi) + offset
        plt.title(rf"$N{from_north:6.2f}^\circ \rightarrow$")
        plt.xlabel('Offset [km]')
        plt.ylabel('Depth [km]')

        if not self.figures_only:
            plt.show(block=False)

        plt.savefig(os.path.join(self.outputdir,
                                 f"slice_{name}_"
                                 f"lat{int(self.latitude+90.0)}_"
                                 f"lon{int(self.longitude+180.0)}_"
                                 f"N{int(from_north)}_{self.meshname}"
                                 f".{self.fmt}"))
