# %%
import re
import uio_tools.uio as uio
import uio_tools.eosinter as eosinter
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from typing import List, Union
import glob
import os
from scipy.interpolate import interp1d, PchipInterpolator

# print('\n'.join(plt.style.available))
plt.style.use('standard-scientific')

sigma_sb = 5.67050e-05  # Stefan-Boltzmann constant, [erg cm^-2 s^-1 K^-4]

# Averaging functions


def spatial_mean(arr: np.ndarray, out_axis=0):
  # Average over 2 spatial axes, 'out_axis' is the axis NOT averaged over
  # e.g. if you want to average over axes (1, 2) (conventionally averaging
  # across (x, y)), pass in 'out_axis = 0'
  axis = [0, 1, 2]
  axis.pop(out_axis)
  return np.mean(arr, axis=tuple(axis))


# -------------------------------------------------------------------------
# Classes for storing plot types
# -------------------------------------------------------------------------


# class UIOData():
#   # Default behaviour is to analyse plots in a specific directory since I
#   # organise separate runs in separate directories
#   def __init__(self, model_path: str, eos_path='') -> None:
#     # Properties
#     self.dataset = None
#     self.box = None
#     self.id = model_path.split('/')[-1].split('.')[0]
#     self.box_keys = None
#     self.eos_path = eos_path if eos_path else None
#     self.eos = None

#     self.model_path = model_path
#     self.model = None
#     self.model_num = '000'

#     self.snap_idx = 0
#     self.first_snap_idx = 0  # unnecessary, it's always 0
#     self.final_snap_idx = 0

#     self.min_max_dict = {}  # min-max tuple for each key

#     # Load the specified model
#     if model_path is not None:
#       self.load_model(model_path)

#     if self.eos_path is not None:
#       self.load_eos()

#     # Save initial and final times
#     self.initial_time = self.get_time()
#     self.final_snapshot()
#     self.final_time = self.get_time()
#     self.first_snapshot()

#   def load_eos(self):
#     # Load the supplied equation of state file
#     self.eos = eosinter.EosInter(self.eos_path)

#   def load_model(self, model_path: str):
#     print(f"Loading model at '{model_path}'")
#     model = uio.File(model_path)
#     self.first_snap_idx = 0
#     self.final_snap_idx = len(model.dataset) - 1
#     self.model = model
#     self.model_num = model_path.split(
#         '/')[-1].split('.')[1]  # assumes a certain format
#     self.update_snapshot(0)

#   def print_properties(self):
#     # Simple function to print properties of the instance as well as
#     # available keys for plotting for this file
#     print("=================")
#     print("-----------------")
#     print("UIOPlot instance properties:")
#     print("-----------------")
#     print(f"Model ID:\t{self.id}")
#     print(f"Located at:\t{self.model_path}")
#     print(f"Model number:\t{self.model_num}")
#     print(f"Current snapshot:\t{self.snap_idx} of {self.final_snap_idx}")
#     print("-----------------")
#     print(f"Available keys:")
#     print('\n'.join(self.box_keys))
#     print("-----------------")

#     print("=================")

#   def get_key_name_units(self, key: str):
#     # For a given key, get the quantity's name and unit as stored in UIO
#     name = self.box[key].params['n'].split(" ")[-1]
#     unit = self.box[key].params['u']

#     return name, unit

#   def get_min_max_quantity(self, key):
#     min_quantity, max_quantity = np.min(
#         self[key]), np.max(self[key])

#     return min_quantity, max_quantity

#   def get_dataset_from_snapshot(self, snap_idx: int):
#     return self.model.dataset[snap_idx]

#   # -------------------------------------------------------------------------
#   # Methods for iterating through snapshots
#   # -------------------------------------------------------------------------

#   def first_snapshot(self):
#     self.update_snapshot(0)

#   def final_snapshot(self):
#     self.update_snapshot(self.final_snap_idx)

#   def prev_snapshot(self):
#     # Look for previous model and make it the current model
#     potential_idx = self.snap_idx - 1
#     if potential_idx < self.first_snap_idx:
#       potential_idx = self.first_snap_idx

#     self.update_snapshot(potential_idx)

#   def next_snapshot(self):
#     # Look for next model and make it the current model
#     potential_idx = self.snap_idx + 1
#     if potential_idx >= self.final_snap_idx:
#       potential_idx = self.final_snap_idx

#     self.update_snapshot(potential_idx)

#   def update_snapshot(self, snap_idx: int):
#     # Update 'snap_idx', 'dataset' and 'box' properties
#     self.snap_idx = snap_idx
#     self.dataset = self.model.dataset[self.snap_idx]
#     self.box = self.dataset.box[0]
#     self.box_keys = [key for key in self.box.keys()]
#   # -----------------------------------------------------------------------
#   # Methods that loop over all snapshots to calculate quantities
#   # -----------------------------------------------------------------------

#   def average_quantity_over_snapshots(self, key: str, set_min_max=False):
#     # Average the specified 'key' from 'self.box' across all snapshots
#     # Has kwarg for setting the 'min-max' dict instance for this key so
#     # that we don't have to iterate over all the snapshots to do this again
#     # (useful for average quantity analyses)
#     num_snapshots = self.final_snap_idx + 1
#     snap_idx = self.snap_idx  # store reference before iterating

#     self.first_snapshot()  # load first snapshot
#     avg_quantity = self.box[key].data

#     # Calculate min and max quantities as well
#     if set_min_max:
#       model_min, model_max = np.min(
#           self.box[key].data), np.max(self.box[key].data)

#     # Loop over snapshots
#     for i in range(num_snapshots):
#       data = self.box[key].data
#       avg_quantity += data

#       if set_min_max:
#         min_data, max_data = np.min(data), np.max(data)
#         if min_data < model_min:
#           model_min = min_data
#         if max_data > model_max:
#           model_max = max_data

#       self.next_snapshot()

#     avg_quantity /= num_snapshots

#     self.update_snapshot(snap_idx)  # revert to original snapshot

#     if set_min_max:
#       self.min_max_dict[key] = (model_min, model_max)

#     return avg_quantity

#   def get_quantities_over_snapshots(self, keys: List[str]):
#     # Create a list of quantity 'key' across all snapshots
#     num_snapshots = self.final_snap_idx + 1
#     snap_idx = self.snap_idx  # store reference before iterating
#     self.first_snapshot()

#     output_quantities = []
#     for i in range(num_snapshots):
#       if len(keys) == 1:
#         quantities = self[keys[0]]
#       else:
#         quantities = [self[key] for key in keys]
#       output_quantities.append(quantities)

#       self.next_snapshot()

#     self.update_snapshot(snap_idx)  # revert to original snapshot

#     return np.array(output_quantities)

#   def min_max_quantity_over_snapshots(self, key: str):
#     # Get min & max of specified 'key' from 'self.box' across all snapshots
#     snap_idx = self.snap_idx  # store reference before iterating

#     self.first_snapshot()  # load first snapshot
#     model_min, model_max = self.get_min_max_quantity(key)

#     # Loop over snapshots
#     for i in range(self.final_snap_idx + 1):
#       data = self[key]
#       min_data, max_data = np.min(data), np.max(data)
#       if min_data < model_min:
#         model_min = min_data
#       if max_data > model_max:
#         model_max = max_data

#       self.next_snapshot()

#     self.update_snapshot(snap_idx)  # revert to original snapshot

#     # Set dict and return? Is this confusing because of the side effect?
#     self.min_max_dict[key] = (model_min, model_max)
#     return model_min, model_max

#   # -----------------------------------------------------------------------
#   # Data manipulation functions
#   # -----------------------------------------------------------------------
#   # Change this function to a dict that gets derived quantities, ignore the
#   # standard 'else' case because that's obvious
#   # Just use the dict to call the correct functions, in the functions, get
#   # the quantities necessary instead of passing them in

#   def __getitem__(self, key: str) -> np.ndarray:
#     # !!!
#     # REFACTOR ALERT: Change all the if/elif/else blocks into a big dict
#     # that holds the different quantities to get (how do I do derived?)
#     # !!!
#     # Get quantity from current model's current snapshot's box
#     data = None
#     opts = ['+', '-', '*', '/', '^']

#     # Check if 'key' has any operators in it
#     # Need to find which operators are here!
#     if any(opt in key for opt in opts):
#       # Split key by operators
#       pattern = r'|'.join([f"\{opt}" for opt in opts])
#       opt_keys = [o_key.strip() for o_key in re.split(pattern, key)]

#       # Replac expressions in 'key'
#       for opt_key in opt_keys:
#         key = key.replace(opt_key, f"self['{opt_key}']")

#       # Evaluate key
#       data = eval(key)
#     else:
#       # Derived quantities
#       if key.lower() in ['x', 'xc1', 'y', 'xc2', 'z', 'xc3', 'xc_all']:
#         x, y, z = self.get_x_vectors(self.box)

#         data_dict = {
#             'x': x,
#             'xc1': x,
#             'y': y,
#             'xc2': y,
#             'z': z,
#             'xc3': z,
#             'xc_all': [x, y, z]
#         }

#         data = data_dict[key]

#       # Here, the first two letters specify the grid and the third letter
#       # (if available) specifies which grid to return
#       elif key.lower() in ['xy_grid', 'xz_grid', 'yz_grid',  # both grids
#                            # single grids
#                            'xyx_grid', 'xyy_grid',
#                            'xzx_grid', 'xzz_grid',
#                            'yzy_grid',  'yzz_grid']:
#         x, y, z = self.get_x_vectors(self.box)

#         data_dict = {
#             'xy_grid': np.meshgrid(x, y),
#             'xz_grid': np.meshgrid(x, z),
#             'yz_grid': np.meshgrid(y, z),
#             'xyx_grid': np.meshgrid(x, y)[0],
#             'xyy_grid': np.meshgrid(x, y)[1],
#             'xzx_grid': np.meshgrid(x, z)[0],
#             'xzz_grid': np.meshgrid(x, z)[1],
#             'yzy_grid': np.meshgrid(y, z)[0],
#             'yzz_grid': np.meshgrid(y, z)[1],
#         }

#         data = data_dict[key]

#       elif key.lower() in ['vx', 'vy', 'vz', 'v', 'velocities']:
#         v1, v2, v3 = self.get_velocity_vectors(self.box)

#         data_dict = {
#             'vx': v1,
#             'vy': v2,
#             'vz': v3,
#             'v': [v1, v2, v3],
#             'velocities': [v1, v2, v3]
#         }

#         data = data_dict[key]

#       elif key.lower() in ['t', 'time', 'modeltime']:
#         data = self.get_time()

#       elif key.lower() == 'kinetic energy':
#         # Calculate kinetic energy with density and velocities
#         density = self.box['rho'].data
#         v1, v2, v3 = self.get_velocity_vectors(self.box)

#         # Full 3D kinetic energy
#         data = calculate_kinetic_energy(density, v1, v2, v3)

#       elif key.lower() == 'momentum':
#         # Calculate momentum with density and velocities
#         density = self.box['rho'].data
#         v1, v2, v3 = self.get_velocity_vectors(self.box)

#         # Full 3D momentum
#         data = calculate_momentum(density, v1, v2, v3)

#       elif key.lower() in ['time', 'modeltime']:
#         data = self.get_time()

#       # EOS quantities
#       elif key.lower() in ['temperature', 'pressure', 'entropy']:
#         rho, ei = self.box['rho'].data, self.box['ei'].data
#         data = self.eos.STP(rho, ei, quantity=key)

#       # Box quantities
#       else:
#         data = self.box[key].data

#     return data.squeeze()

#   # Package below into a separate UIO script (or keep this entire thing
#   # as 'uio_utilities')
#   # -----------------------------------------------------------------------
#   # Convenience functions for box data
#   # -----------------------------------------------------------------------
#   def get_x_vectors(self, box=None, convert_km=True) -> Union[
#           np.ndarray, np.array, np.array]:
#     if not box:
#       box = self.box
#     # Get 'x', 'y', 'z' from box and squeeze empty dimensions
#     x = box['xc1'].data.squeeze()
#     y = box['xc2'].data.squeeze()
#     z = box['xc3'].data.squeeze()

#     if convert_km:
#       x /= 1e5
#       y /= 1e5
#       z /= 1e5

#     return x, y, z

#   def get_velocity_vectors(self, box=None, convert_km=True) -> Union[
#           np.ndarray, np.ndarray, np.ndarray]:
#     if not box:
#       box = self.box
#     # Note 3D array indexing is reversed because of Python vs IDL
#     v1 = box['v1'].data
#     v2 = box['v2'].data
#     v3 = box['v3'].data

#     if convert_km:
#       v1 /= 1e5
#       v2 /= 1e5
#       v3 /= 1e5

#     return v1, v2, v3

#   def get_time(self, dataset=None):
#     if not dataset:
#       dataset = self.dataset
#     return dataset['modeltime'].data

#   def get_time_difference(self, snap_idx1: int, snap_idx2: int):
#     # Calculate the time difference between two snapshots
#     time1 = self.get_dataset_from_snapshot(snap_idx1)
#     time2 = self.get_dataset_from_snapshot(snap_idx2)

#     return time2 - time1


# class FullData(UIOData):
#   def __init__(self, model_path: str, eos_file='') -> None:
#     super().__init__(model_path, eos_file)
#     self.x = self['xc1']
#     self.y = self['xc2']
#     self.z = self['xc3']

#     # Standard x-z grid
#     self.set_grid(*initialise_grid(self.x, self.y, self.z, 'xz'))

#   def add_tau(self, mean_data, gravity: float):
#     # 'mean_data: MeanData' should correspond to the same model snapshot
#     # as the full file; needs to be passed in separately since it contains
#     # quantities used to calculate optical depth
#     # 'gravity' defined on linear scale
#     mean_data.add_qlmean_quantities(gravity)
#     self.tau = mean_data.tau

#   def quantity_at_tau_val(self, quantity: str, tau_val: float):
#     # 'tau_val' on linear scale
#     return interp1d(self.tau, self[quantity])(tau_val)

#   def z_zero_point(self, kind='tau', tau_val=1):
#     # set the 'z' zero point
#     # 'kind' is one of 'tau' (set based on optical depth,
#     # requires a 'tau_val')
#     # or 'bottom' (sets the bottom of the grid to be the zero point)
#     kinds = {
#         'tau': lambda x: self.quantity_at_tau_val('xc3', x),
#         'bottom': lambda x: None
#     }

#     zero_point = kinds[kind](tau_val)
#     return zero_point

#   # -------------------------------------------------------------------------
#   # Setters
#   # -------------------------------------------------------------------------

#   def set_grid(self, X_grid: np.ndarray, Y_grid: np.ndarray):
#     self.X_grid = X_grid
#     self.Y_grid = Y_grid

#   def set_z_zero_point(self, zero_point=None):
#     # Set the 'z' zero point
#     if not zero_point:
#       # Use tau=1 as standard
#       zero_point = self.z_zero_point()

#     self.z = self['z'] - zero_point
#   # -----------------------------------------------------------------------
#   # Plotting methods
#   # -----------------------------------------------------------------------

#   def plot_heatmap(self, ax, plot_values, log_quantity=False, title=None,
#                    plot_type='image', add_cbar=True, cbar=None,
#                    cmap='jet', cbar_label=None, cbar_label_pos='left',
#                    origin='lower',
#                    vlimits=None):
#     # Normalise data to range [0, 1] before applying colours
#     if log_quantity:
#       plot_values = np.log10(plot_values)

#     if vlimits:  # should be a tuple of (vmin, vmax)
#       vmin, vmax = vlimits
#     else:  # determine min,max from plot values
#       vmin, vmax = plot_values.min(), plot_values.max()

#     norm = Normalize(vmin=vmin, vmax=vmax)

#     # Plot a heatmap using X_, Y_grids and a Z datacube 'plot_values'
#     if plot_type == 'mesh':
#       im = ax.pcolormesh(self.X_grid, self.Y_grid,
#                          plot_values, cmap=cmap, norm=norm)

#     elif plot_type == 'contour':
#       im = ax.contourf(self.X_grid, self.Y_grid,
#                        plot_values, cmap=cmap, norm=norm)

#     elif plot_type == 'image':
#       x_limits = [self.X_grid[0][0], self.X_grid[0][-1]]
#       y_limits = [self.Y_grid[0][0], self.Y_grid[-1][0]]
#       extent = (x_limits[0], x_limits[1], y_limits[0], y_limits[1])

#       im = ax.imshow(plot_values, interpolation='bilinear', origin=origin,
#                      cmap=cmap, norm=norm, extent=extent)

#     else:
#       print(f"Warning: Plot type {plot_type} is not valid. Valid choices are \
#         'mesh', 'contour' and 'image'.")
#       print("Defaulting to 'mesh'.")
#       im = ax.pcolormesh(self.X_grid, self.Y_grid, plot_values, cmap='jet')

#     # Set number of ticks
#     ax.xaxis.set_major_locator(MaxNLocator(5))
#     ax.yaxis.set_major_locator(MaxNLocator(5))

#     if add_cbar:
#       # colorbar same height as heatmap
#       divider = make_axes_locatable(ax)
#       cax = divider.append_axes("right", size="5%", pad=0.05)

#       # Create new colorbar
#       if cbar is None:
#         cbar = ax.figure.colorbar(im, ax=ax, cax=cax)

#       # Set default colorbar ticks
#       ticks = np.linspace(vmin, vmax, num=5)
#       cbar.set_ticks(ticks)
#       cbar.ax.set_yticklabels([f"{t:.2f}" for t in ticks])

#       if cbar_label:
#         cbar_label_positions = {
#             'top': cbar.ax.set_title,
#             'right': cbar.set_label,
#             'bottom': cbar.ax.set_xlabel,
#         }

#         if cbar_label_pos in cbar_label_positions.keys():
#           cbar_label_positions[cbar_label_pos](cbar_label)
#         else:
#           print(f"Error: {cbar_label_pos} is not a valid choice."
#                 "'top', 'right' and 'bottom' are valid choices. Using"
#                 "'right'")
#           cbar_label_positions['right'](cbar_label)

#       if title:
#         ax.set_title(title)

#       return cbar

#   def plot_quantity_heatmap(self, ax: plt.Axes, key: str,
#                             plane='xz', set_z_zero_to_tau_1=False,
#                             title=None, cmap='jet', log_quantity=False,
#                             xlabel=None, ylabel=None,
#                             auto_label_axes=False,
#                             origin='lower',
#                             add_cbar=True, cbar_label=None,
#                             cbar_label_pos='right',
#                             average_snaps=False):
#     # Plot quantity in a certain plane as a heatmap. Axis labels are
#     # determined by 'plane'

#     # Set 'min-max' bounds for a key if they have not yet been set
#     set_min_max = False if key in self.min_max_dict else True

#     # Get data from box and set 'min-max' bounds for the key
#     # Average over all snapshots
#     if average_snaps:
#       data = self.average_quantity_over_snapshots(
#           key, set_min_max=set_min_max)
#     else:
#       data = self[key]
#       if set_min_max:
#         self.min_max_quantity_over_snapshots(key)

#     # Average over axis not in plane
#     if data is not None:
#       if len(data.shape) == 3:
#         # For now, assume IDL indexing since 'data' has not been converted to
#         # Python indexing
#         print(f"Averaging data with key {key} in {plane} plane.")
#         avg_data = average_data_to_plane(data, plane, is_idl_idx=True)

#       else:
#         print(
#             f"Data with key {key} has shape {data.shape}, skipping plane average.")
#         avg_data = data

#     else:
#       print(f"Error: Could not compute data for key {key}.")
#       return None

#     # Plot heatmap of data
#     # Make grid points for plotting plane
#     if set_z_zero_to_tau_1:
#       print("Setting z zero point")
#       self.set_z_zero_point()

#     self.set_grid(*initialise_grid(self.x, self.y, self.z, plane))

#     cbar = self.plot_heatmap(ax, avg_data, log_quantity=log_quantity,
#                              plot_type='image', origin=origin,
#                              title=title, cmap=cmap, add_cbar=add_cbar,
#                              cbar_label=cbar_label,
#                              cbar_label_pos=cbar_label_pos)

#     # Set labels
#     if auto_label_axes:
#       xlabel, ylabel = [f"{item} [km]" for item in list(plane)]

#     if xlabel is not None:
#       ax.set_xlabel(xlabel)

#     if ylabel is not None:
#       ax.set_ylabel(ylabel)

#     return cbar


# class MeanData(UIOData):
#   def __init__(self, model_path: str, eos_file='', gravity=None) -> None:
#     # gravity='infer' : infer from model name
#     super().__init__(model_path, eos_file)
#     self.z = self.get_z()

#     if gravity:  # defined on linear scale
#       self.add_qlmean_quantities(gravity)

#   def add_qlmean_quantities(self, gravity: float):
#     # 'gravity' defined on linear scale
#     self.set_box(2)
#     z = self.get_x_vectors(convert_km=False)[2].squeeze()
#     xcm = self['kapparho_xmean']
#     tau0 = xcm[-1] / self['rho_xmean'][-1] * self['p_xmean'][-1] / gravity

#     # # Reverse 'xcm' to be monotonically decreasing to get correct integral
#     interp = PchipInterpolator(z, xcm[::-1])
#     # Integrate kappa-rho over z to get tau
#     # self.tau = cumtrapz(self['kapparho_xmean'], z)
#     self.frad = self['ferb_xmean']
#     self.tau = tau0 + interp.antiderivative()(z)

#   def quantity_at_tau_val(self, quantity: str, tau_val: float):
#     # 'tau_val' on linear scale
#     return interp1d(self.tau, self[quantity])(tau_val)

#   def z_zero_point(self, kind='tau', tau_val=1):
#     # set the 'z' zero point
#     # 'kind' is one of 'tau' (set based on optical depth,
#     # requires a 'tau_val')
#     # or 'bottom' (sets the bottom of the grid to be the zero point)
#     kinds = {
#         'tau': lambda x: self.quantity_at_tau_val('xc3', x),
#         'bottom': lambda x: None
#     }

#     zero_point = kinds[kind](tau_val)
#     return zero_point

#   def set_z_zero_point(self, zero_point=None):
#     # Set the 'z' zero point
#     if not zero_point:
#       # Use tau=1 as standard
#       zero_point = self.z_zero_point()

#     self.z -= zero_point

#   def first_snapshot(self, box_idx=2):
#     self.update_snapshot(0, box_idx)

#   def final_snapshot(self, box_idx=2):
#     self.update_snapshot(self.final_snap_idx, box_idx)

#   def prev_snapshot(self, box_idx=2):
#     # Look for previous model and make it the current model
#     potential_idx = self.snap_idx - 1
#     if potential_idx < self.first_snap_idx:
#       potential_idx = self.first_snap_idx

#     self.update_snapshot(potential_idx, box_idx)

#   def next_snapshot(self, box_idx=2):
#     # Look for next model and make it the current model
#     potential_idx = self.snap_idx + 1
#     if potential_idx >= self.final_snap_idx:
#       potential_idx = self.final_snap_idx

#     self.update_snapshot(potential_idx, box_idx)

#   def update_snapshot(self, snap_idx: int, box_idx=2):
#     # Update 'snap_idx', 'dataset' and 'box' properties
#     self.snap_idx = snap_idx
#     self.dataset = self.model.dataset[self.snap_idx]
#     self.box = self.dataset.box[box_idx]  # '2' is z3, standard
#     self.box_keys = [key for key in self.box.keys()]

#   def set_box(self, box_idx: int):
#     self.box = self.dataset.box[box_idx]

#   def __getitem__(self, key: str) -> np.ndarray:
#     # Get quantity from current model's current snapshot's boxes
#     data = None
#     opts = ['+', '-', '*', '/', '^']

#     if any(opt in key for opt in opts):
#       # Split key by operators
#       pattern = r'|'.join([f"\{opt}" for opt in opts])
#       opt_keys = [o_key.strip() for o_key in re.split(pattern, key)]

#       # Replace expressions in 'key'
#       # Warning: if the 'opt_key' is present more than once, it breaks!
#       for opt_key in opt_keys:
#         key = key.replace(opt_key, f"self['{opt_key}']")

#       # Evaluate key
#       data = eval(key)

#     else:
#       # Derived quantities
#       if key.lower() == 'kinetic energy':
#         # Calculate kinetic energy with density and velocities
#         density = self.box['rho_xmean'].data
#         v1, v2, v3 = self.get_velocity_vectors(self.box)

#         # Full 3D kinetic energy
#         data = calculate_kinetic_energy(density, v1, v2, v3)

#       elif key.lower() == 'momentum':
#         # Calculate momentum with density and velocities
#         density = self.box['rho_xmean'].data
#         v1, v2, v3 = self.get_velocity_vectors(self.box)

#         # Full 3D momentum
#         data = calculate_momentum(density, v1, v2, v3)

#       # # EOS quantities  # Don't use for mean files!
#       # eos_quantities = ['temperature', 'pressure', 'entropy']
#       # eos_quantities.extend(
#       #     [f"{quantity}_xmean" for quantity in eos_quantities])

#       # if key.lower() in eos_quantities:
#       #   rho, ei = self.box['rho_xmean'].data, self.box['ei_xmean'].data
#       #   data = self.eos.STP(rho, ei, quantity=key)

#       # Box quantities
#       else:
#         data = self.box[key].data

#     return data.squeeze()

#   def min_max_quantity_over_snapshots(self, key: str):
#     # Get min & max of specified 'key' from 'self.box' across all snapshots
#     snap_idx = self.snap_idx  # store reference before iterating

#     self.first_snapshot()  # load first snapshot
#     model_min, model_max = self.get_min_max_quantity(key)

#     # Loop over snapshots
#     for i in range(self.final_snap_idx + 1):
#       data = self[key]
#       min_data, max_data = np.min(data), np.max(data)
#       if min_data < model_min:
#         model_min = min_data
#       if max_data > model_max:
#         model_max = max_data

#       self.next_snapshot()

#     self.update_snapshot(snap_idx)  # revert to original snapshot

#     # Set dict and return? Is this confusing because of the side effect?
#     self.min_max_dict[key] = (model_min, model_max)
#     return model_min, model_max

#   def get_min_max_quantity(self, key):
#     min_quantity, max_quantity = np.min(
#         self[key]), np.max(self[key])

#     return min_quantity, max_quantity

#   def get_z(self, unit='km'):
#     # Get the 'z' direction from the box, assuming we are in box z3 (idx 2)
#     z = self.box['xc3'].data.squeeze()

#     if unit == 'km':
#       z /= 1e5  # Change from 'cm' to 'km'

#     elif unit == 'm':
#       z /= 1e3  # Change from 'cm' to 'km'

#     return z
#   # -----------------------------------------------------------------------
#   # Plotting methods
#   # -----------------------------------------------------------------------

#   def plot_keys(self, ax: plt.Axes, key1: str, key2: str,
#                 title=None, average_snaps=False,
#                 label=None, colour=None, ls='-',
#                 xlabel=None, ylabel=None,
#                 log_key_1=False, log_key_2=False,
#                 xscale=None, yscale=None, as_abu=False):
#     # Plot 'key1' vs 'key2' in the z-direction (no implementation for x,y)
#     self.set_box(2)  # set box to 'z3' (z-axis)

#     data_list = []
#     log_keys = [log_key_1, log_key_2]
#     # Set 'min-max' bounds for a key if they have not yet been set
#     for key, log_key in zip([key1, key2], log_keys):
#       set_min_max = False if key in self.min_max_dict else True

#       if key.lower() == 'tau':
#         data = np.log10(self.tau)

#       elif key.lower() == 'z':
#         data = self.z

#       elif not key.endswith('_xmean'):
#         xkey = key + '_xmean'

#       else:
#         xkey = key

#       # Get data from box and set 'min-max' bounds for the key
#       # Average over all snapshots
#       if not key.lower() == 'tau' and not key.lower() == 'z':
#         if average_snaps:
#           data = self.average_quantity_over_snapshots(
#               xkey, set_min_max=set_min_max)
#         else:
#           data = self[key]

#         if set_min_max:
#           if key in ['temperature', 'pressure', 'entropy']:
#             self.min_max_quantity_over_snapshots(key)
#           else:
#             self.min_max_quantity_over_snapshots(xkey)

#         # Convert QUC to abundances if specified
#         # Assumes QUC001 is the hydrogen number density
#         if key.startswith('quc') and as_abu:
#           data = number_density_to_abundance(data, self['quc001_xmean'])

#       if log_key:
#         data = np.log10(data)

#       data_list.append(data)

#     # Plot key1 data vs key2 data
#     if colour:
#       ax.plot(data_list[0], data_list[1], c=colour, label=label, ls=ls)
#     else:
#       ax.plot(data_list[0], data_list[1], label=label, ls=ls)

#     if title:
#       ax.set_title(title)

#     if xlabel:
#       ax.set_xlabel(xlabel)
#     if ylabel:
#       ax.set_ylabel(ylabel)

#     if xscale:
#       ax.set_xscale(xscale)
#     if yscale:
#       ax.set_yscale(yscale)

#   def z_plot(self, ax: plt.Axes, key: str,
#              title=None, average_snaps=False,
#              label=None, colour=None, ls='-',
#              xlabel=None, ylabel=None,
#              xscale='linear', yscale='linear',
#              as_abu=False):
#     # Convenience method for plotting keys against 'z'; just light wrapper
#     # around 'self.plot_keys()'
#     if not xlabel:
#       xlabel = 'z [km]'
#     self.plot_keys(ax, key1='z', key2=key, title=title,
#                    average_snaps=average_snaps,
#                    label=label, colour=colour, ls=ls,
#                    xlabel=xlabel, ylabel=ylabel,
#                    xscale=xscale, yscale=yscale, as_abu=as_abu)

#   def tau_plot(self, ax: plt.Axes, key: str,
#                title=None, average_snaps=False,
#                label=None, colour=None, ls='-',
#                xlabel=None, ylabel=None,
#                xscale='linear', yscale='linear',
#                log_quantity=False,
#                as_abu=False):
#     # Convenience method for plotting keys against 'log_tau'; just light
#     # wrapper around 'self.plot_keys()'
#     if not ylabel:
#       ylabel = r'$\log{tau}$'

#     self.plot_keys(ax, key1=key, key2='tau', title=title,
#                    average_snaps=average_snaps,
#                    log_key_1=log_quantity, log_key_2=False,
#                    label=label, colour=colour, ls=ls,
#                    xlabel=xlabel, ylabel=ylabel,
#                    xscale=xscale, yscale=yscale, as_abu=as_abu)


# class UIOLoader():
#   # List holding filepaths for UIOPlot() instances
#   # Contains helper methods to load next, previous models and holds state
#   # of plotting variables to be the same across all models
#   def __init__(self, model_directory: str, file_type='full', eos_file='') -> None:
#     self.model_dir = model_directory
#     self.model_files = None
#     self.file_type = file_type  # either 'full' or 'mean'
#     self.id = ''
#     self.idx = 0  # current idx in 'model_files'
#     self.current_model = None
#     self.current_model_path = None
#     self.eos_file = eos_file if eos_file else None

#     self.load_model_files()

#   def load_model(self):
#     # Load the current model
#     self.current_model_path = f"{self.model_files[self.idx]}"
#     if self.file_type == 'full':
#       self.current_model = FullData(self.current_model_path, self.eos_file)
#     elif self.file_type == 'mean':
#       self.current_model = MeanData(self.current_model_path, self.eos_file)
#     else:
#       print("Invalid file type. Valid types: 'full', 'mean'")

#   def load_first_model(self):
#     self.idx = 0
#     self.load_model()

#   def load_final_model(self):
#     self.idx = len(self.model_files) - 1
#     self.load_model()

#   def load_next_model(self):
#     new_idx = self.idx + 1
#     # Could return something if we're at the first/last model to auto-break
#     # out of external loops?
#     if new_idx >= len(self.model_files):
#       print("Already at last model. Reloading model.")
#     else:
#       self.idx = new_idx

#     self.load_model()

#   def load_prev_model(self):
#     new_idx = self.idx - 1
#     if new_idx <= 0:
#       print("Already at first model. Reloading model.")
#     else:
#       self.idx = new_idx

#     self.load_model()

#   def load_model_files(self):
#     # Parse files for identifier string, earliest and latest model
#     # Assumes all model files have structure "identifier.XX.full"
#     # where XX is a digit, often 00 -> 79
#     # Use first file to determine identifier
#     # Assumes there is only one 'identifier' in the directory, which is an
#     # obvious problem!
#     files = glob.glob(os.path.join(self.model_dir, f"*.{self.file_type}"))
#     identifier = files[0].split('/')[-1].split('.')[0]
#     filtered_files = [file_ for file_ in files if identifier in file_]
#     filtered_files.sort()

#     # Save filtered files
#     self.model_files = filtered_files
#     self.id = identifier

#     print(f"ID: {self.id}. {len(filtered_files)} files found.")

#     # Set current model to first model
#     self.idx = 0
#     self.load_model()

#   def print_files(self):
#     print(f"{self.model_dir} contains:")
#     print('\n'.join([f"\tIndex {i}: {file_}" for i,
#                      file_ in enumerate(self.model_files)]))

#   # -------------------------------------------------------------------------
#   # Functions that loop over models
#   # -------------------------------------------------------------------------
#   def get_quantities_over_models(self, keys: List[str]):
#     # Create a list of quantity 'key' across all models
#     num_snapshots = len(self.model_files)
#     model_idx = self.idx  # store reference before iterating
#     self.load_first_model()

#     quantities = []
#     for i in range(num_snapshots):
#       data = self.current_model.get_quantities_over_snapshots(keys)
#       quantities.extend(data)

#       self.load_next_model()

#     self.idx = model_idx
#     self.load_model()  # revert to original snapshot

#     return np.array(quantities)


# -------------------------------------------------------------------------
# Quantity calculation functions
# -------------------------------------------------------------------------


def calculate_kinetic_energy(density: np.ndarray,
                             v1: np.ndarray, v2: np.ndarray, v3: np.ndarray) -> np.ndarray:
  # Calculate kinetic energy (density, default [erg / cm^3])
  # All input arrays must have the same dimensions
  kinetic_energy = 0.5 * density * (v1**2 + v2**2 + v3**2)

  return kinetic_energy


def calculate_momentum(density: np.ndarray,
                       v1: np.ndarray, v2:  np.ndarray, v3: np.ndarray) -> np.ndarray:
  # Calculate momentum (default [g cm^-2 s^-1])
  # All input arrays must have the same dimensions
  momentum = density * np.sqrt(v1**2 + v2**2 + v3**2)

  return momentum


def initialise_grid(x: np.ndarray, y: np.ndarray, z: np.ndarray,
                    plane: str) -> Union[np.ndarray, np.ndarray]:
  # Expects 'x', 'y', 'z' to properly represent coordinate space, meaning
  # the conversion between index-ordering has already been performed
  # Creates a mesh grid as specified
  plane_grid = {
      'xy': np.meshgrid(x, y),
      'xz': np.meshgrid(x, z),
      'yz': np.meshgrid(y, z)
  }

  X_grid, Y_grid = plane_grid[plane]

  return X_grid, Y_grid


def average_data_to_plane(data: np.ndarray, plane: str, is_idl_idx=True):
  # Average a 3D array to a 2D plane by averaging over the axis
  # not in 'plane'.
  # Assumes IDL indexing if 'is_idl_idx' (opposite to Python), i.e., 'data'
  # is ordered {z, y, x} instead of {x, y, z}
  if plane == 'xy':
    # Average over 'z'
    axis = 0 if is_idl_idx else 2

  elif plane == 'xz':
    # Average over 'y'
    axis = 1

  elif plane == 'yz':
    # Average over 'x'
    axis = 2 if is_idl_idx else 0

  else:
    print(f"Error: {plane} is not a valid option. Valid planes are 'xy', \
      'xz', 'yz'.")

    return None

  # Compute average over axis
  avg_data = np.average(data, axis=axis)

  return avg_data

# Refactor: Move to chemistry utilities


def number_density_to_abundance(number_density: float,
                                hydrogen_number_density: float):
  # A_x = log(n_x) - log(n_H) + 12
  return np.log10(number_density) - np.log10(hydrogen_number_density) + 12
