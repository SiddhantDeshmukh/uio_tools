# %%
from uio_utils import initialise_grid, average_data_to_plane
from .quantities import *

import re
import uio
import eosinter
import opta
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from typing import List, Union
from scipy.interpolate import interp1d, PchipInterpolator
# print('\n'.join(plt.style.available))
plt.style.use('standard-scientific')

sigma_sb = 5.67050e-05  # Stefan-Boltzmann constant, [erg cm^-2 s^-1 K^-4]

# -------------------------------------------------------------------------
# Classes for storing plot types
# -------------------------------------------------------------------------


class UIOData():
  # Default behaviour is to analyse plots in a specific directory since I
  # organise separate runs in separate directories
  # Contains Full and Mean file data
  def __init__(self, model_path: str, eos_file='', opta_file='', par_file='') -> None:
    # Properties
    self.dataset = None
    self.box = None
    self.mean_box = None
    self.id = model_path.split('/')[-1].split('.')[0]
    self.box_keys = None
    self.eos_path = eos_file if eos_file else None
    self.opta_path = opta_file if opta_file else None
    self.par_path = par_file if par_file else None
    self.eos = None
    self.opta = None
    self.par = None

    # Model quantities
    self.dimension = None
    self.temperature = None
    self.gravity = None
    self.carbon_enhancement = None
    self.carbon_oxygen_ratio = None
    self.chemical_network = None

    self.model_path = model_path
    self.full = None
    self.mean = None
    self.model_num = ''

    self.snap_idx = 0
    self.first_snap_idx = 0  # unnecessary, it's always 0
    self.final_snap_idx = 0

    self.min_max_dict = {}  # min-max tuple for each key

    # Infer quantities from model name
    self.infer_quantities()

    if self.eos_path:
      self.load_eos()

    if self.opta_path:
      self.load_opta()

    # Load the specified model
    self.load_model(model_path)

    # Save initial and final times
    self.initial_time = self.get_time()
    self.final_snapshot()
    self.final_time = self.get_time()
    self.first_snapshot()

    # Call this upon updating snapshot!
    self.initialise_quantities()

    # Optical depth
    if self.opta:
      self.tau = self['tau']

    self.add_qlmean_quantities()

    # Standard x-z grid
    self.set_grid(*initialise_grid(self.x, self.y, self.z, 'xz'))

  def initialise_quantities(self):
    # Initialise basic quantities for the model
    # x, y, z arrays from full file
    self.x, self.y, self.z = self.get_x_vectors(self.box)
    self.v1, self.v2, self.v3 = self.get_velocity_vectors(self.box)

    if self.eos:
      self.initialise_eos_quantities()

    if self.opta:
      self.initialise_opta_quantities()

    self.initialise_data_dict()

  def initialise_data_dict(self):
    # Initialise data dict
    x, y, z = self.x, self.y, self.z
    v1, v2, v3 = self.v1, self.v2, self.v3
    self.data_dict = {
        # TODO:
        # Turn everything into lambda statements and evaluate on return
        # Position (grid cells)
        'x': lambda: x,
        'xc1': lambda: x,
        'y': lambda: y,
        'xc2': lambda: y,
        'z': lambda: z,
        'xc3': lambda: z,
        'xc_all': lambda: [x, y, z],
        # Basic quantities from box
        'rho': lambda: self.get_box_quantity('rho'),
        'density': lambda: self.get_box_quantity('rho'),
        'ei': lambda: self.get_box_quantity('ei'),
        'energy': lambda: self.get_box_quantity('ei'),
        # 2D grids
        'xy_grid': lambda: np.meshgrid(x, y),
        'xz_grid': lambda: np.meshgrid(x, z),
        'yz_grid': lambda: np.meshgrid(y, z),
        # 2D grids, pick one
        'xyx_grid': lambda: np.meshgrid(x, y)[0],
        'xyy_grid': lambda: np.meshgrid(x, y)[1],
        'xzx_grid': lambda: np.meshgrid(x, z)[0],
        'xzz_grid': lambda: np.meshgrid(x, z)[1],
        'yzy_grid': lambda: np.meshgrid(y, z)[0],
        'yzz_grid': lambda: np.meshgrid(y, z)[1],
        # Velocities
        'vx': lambda: v1,
        'v1': lambda: v1,
        'vy': lambda: v2,
        'v2': lambda: v2,
        'vz': lambda: v3,
        'v3': lambda: v3,
        'v': lambda: [v1, v2, v3],
        'velocities': lambda: [v1, v2, v3],
        'v_squared': lambda: (self.v1**2 + self.v2**2 + self.v3**2),
        'abs_v': lambda: np.sqrt(self.v_squared),
        # Time
        't': self.get_time,
        'time': self.get_time,
        'modeltime': self.get_time,
        'kinetic energy': lambda: calculate_kinetic_energy(self.rho, self.v1, self.v2, self.v3),
        'ke': lambda: calculate_kinetic_energy(self.rho, self.v1, self.v2, self.v3),
        'momentum': lambda: calculate_momentum(self.rho, self.v1, self.v2, self.v3),
        # EOS quantities
        'pressure': lambda: self.P,
        'temperature': lambda: self.T,
        'entropy': lambda: self.S,
        'dpdrho': lambda: self.dPdrho,
        'dpdei': lambda: self.dPdei,
        'dtdei': lambda: self.dTdei,
        # OPTA quantities
        'opacity': lambda: self.kappa,
        'optical depth': lambda: self.tau,
        # Thermodynamic quantities
        'gamma_1': lambda: (self.rho / self.P) * self.dPdrho + (1 / self.rho) * self.dPdei,
        'gamma_3': lambda: 1 + (1 / self.rho) * self.dPdei,
        'grad_t': lambda: (self.gamma_3 - 1) / self.gamma_1,
        'cv_prime': lambda: (self.P / (self.rho * self.T)) * (1 / self.gamma_3 - 1),
        'cp_prime': lambda: (self.P / (self.rho * self.T)) * (self.gamma_1 / (self.gamma_3 - 1)),
        # Hydrodynamic quantities
        'c_s': lambda: np.sqrt(self.gamma_1 * self.P / self.rho),
        'M_cs': lambda: self.abs_v / self.c_s,

    }

  def initialise_eos_quantities(self):
    rho, ei = self.rho, self.ei
    self.P, self.dPdrho, self.dPdei = self.eos.Pall(rho, ei)
    self.T, self.dTdei = self.eos.Tall(rho, ei)
    self.S = self.eos.STP(rho, ei, quantity='e')
    self.dTdrho = (self.T / rho**2) * self.dPdei - \
        (self.P / rho**2) * self.dTdei

  def initialise_opta_quantities(self):
    rho, P, T = self.rho, self.P, self.T
    self.kappa = self.opta.kappa(T, P)
    self.tau = self.opta.tau(rho, P=P, T=T, z=self.z, axis=0)

  def get_box_quantity(self, key: str):
    # Get quantity from 'box'
    return self.box[key].data

  def __getitem__(self, key: str) -> np.ndarray:
    # !!!
    # REFACTOR ALERT: Change all the if/elif/else blocks into a big dict
    # that holds the different quantities to get (how do I do derived?)
    # !!!
    # Get quantity from current model's current snapshot's box
    data = None
    opts = ['+', '-', '*', '/', '^']

    # Check if 'key' has any operators in it
    # Need to find which operators are here!
    if any(opt in key for opt in opts):
      # Split key by operators
      pattern = r'|'.join([f"\{opt}" for opt in opts])
      opt_keys = [o_key.strip() for o_key in re.split(pattern, key)]

      # Replace expressions in 'key'
      for opt_key in opt_keys:
        key = key.replace(opt_key, f"self['{opt_key}']")

      # Evaluate key
      data = eval(key)
    else:
      # Mean quantities
      if key.endswith('xmean'):
        data = self.mean_box[key].data

      else:
        data = self.data_dict[key]
    # return data.squeeze()
    return data

  def __getattr__(self, name):
    return self[name]

  def quantity_from_model_id(self, search_character: str, num_characters: int):
    # Given the character to search for and the number of characters
    # afterwards
    model = self.model_path.split('/')[-1].split('.')[0]
    quantity = model.split(search_character)[1][:num_characters]

    return quantity

  def infer_quantities(self):
    # Model number is after identifier and before file extension
    self.model_num = self.model_path.split('/')[-1].split('.')[1]

    # From the model path, infer temperature, gravity, metallicity,
    # and carbon enhancement, carbon-to-oxygen ratio and chemistry if present
    model = self.model_path.split('/')[-1].split('.')[0]
    dimension = int(re.findall(r"d[0-9]", model)[0][1])
    temperature = float(re.findall(r"t[0-9]+", model)[0][1:] + '00')
    gravity = re.findall(r"g[0-9]+", model)[0][1:]
    gravity = float(f"{gravity[0]}.{gravity[1]}")

    self.dimension = dimension
    self.temperature = temperature
    self.gravity = 10**gravity  # linear scale!

    # Check for carbon enhancement and carbon-to-oxygen ratio
    if 'c' in model:
      carbon_enhancement = model.split('c')[1][:3]
      carbon_enhancement = f"{carbon_enhancement[0]}.{carbon_enhancement[1:]}"
      self.carbon_enhancement = carbon_enhancement

    if 'co' in model:
      carbon_oxygen_ratio = model.split('co')[1][:3]
      carbon_oxygen_ratio = f"{carbon_oxygen_ratio[0]}.{carbon_oxygen_ratio[1:]}"
      self.carbon_oxygen_ratio = carbon_oxygen_ratio

    # Check chemical network
    if '-' in model:
      network = model.split('-')[-1][1:]
      self.chemical_network = network

  # -----------------------------------------------------------------------
  # Setters
  # -----------------------------------------------------------------------

  def set_grid(self, X_grid: np.ndarray, Y_grid: np.ndarray):
    self.X_grid = X_grid
    self.Y_grid = Y_grid

  def set_mean_box(self, box_idx: int):
    self.mean_box = self.mean.dataset[self.snap_idx].box[box_idx]

  def set_dimension(self, dimension: int):
    self.dimension = dimension

  def set_temperature(self, temperature: float):
    self.temperature = temperature

  def set_gravity(self, gravity: float):
    # Conventionally defined on linear scale
    self.gravity = gravity

  def set_carbon_enhancement(self, carbon_enhancement: float):
    self.carbon_enhancement = carbon_enhancement

  def set_carbon_oxygen_ratio(self, carbon_oxygen_ratio: float):
    self.carbon_oxygen_ratio = carbon_oxygen_ratio

  def set_chemical_network(self, network: str):
    self.chemical_network = network

  # -----------------------------------------------------------------------
  # Loaders
  # -----------------------------------------------------------------------

  def load_eos(self):
    # Load the supplied equation of state file
    self.eos = eosinter.EosInter(self.eos_path)

  def load_opta(self):
    # Load supplied opacity table
    self.opta = opta.Opac(self.opta_path)

  def load_model(self, model_path: str):
    model_type = model_path.split('.')[-1]  # either 'full' or 'mean'
    print(f"Loading model at '{model_path}'")

    model = uio.File(model_path)
    # Current (rough) implementation assumes a mean and full file are both
    # present in this directory & are named the same way

    if model_type == 'full':
      self.full = model
      model_path = model_path.replace('.full', '.mean')
      print(f"Loading model at '{model_path}'")

      model = uio.File(model_path)
      self.mean = model

    elif model_type == 'mean':
      self.mean = model
      model_path = model_path.replace('.mean', '.full')
      print(f"Loading model at '{model_path}'")

      model = uio.File(model_path)
      self.full = model

    else:
      print(f"Error: Model of type '{model_type}' not supported, file must\
              be either 'full' or 'mean'.")

    self.first_snap_idx = 0
    self.final_snap_idx = len(self.full.dataset) - 1
    self.model_num = model_path.split(
        '/')[-1].split('.')[1]  # assumes a certain format
    self.update_snapshot(0)

    # If 'mean' model has been loaded, init the qlmean quantities and tau
    if self.mean and self.gravity:
      self.add_qlmean_quantities()

  def add_qlmean_quantities(self):
    # 'gravity' defined on linear scale
    self.set_mean_box(2)
    # z = self.get_x_vectors(convert_km=False)[2].squeeze()
    xcm = self['kapparho_xmean']
    # tau0 = xcm[-1] / self['rho_xmean'][-1] * \
    #     self['p_xmean'][-1] / self.gravity

    # # Reverse 'xcm' to be monotonically decreasing to get correct integral
    # interp = PchipInterpolator(z, xcm[::-1])
    # Integrate kappa-rho over z to get tau
    # self.tau = cumtrapz(self['kapparho_xmean'], z)
    self.frad = self['ferb_xmean']
    # self.tau = tau0 + interp.antiderivative()(z)

  def quantity_at_tau_val(self, quantity: str, tau_val: float):
    # 'tau_val' on linear scale
    tau = np.mean(self.tau.squeeze(), axis=1)
    return interp1d(tau, self[quantity])(tau_val)

  def z_zero_point(self, kind='tau', tau_val=1):
    # set the 'z' zero point
    # 'kind' is one of 'tau' (set based on optical depth,
    # requires a 'tau_val')
    # or 'bottom' (sets the bottom of the grid to be the zero point)
    kinds = {
        'tau': lambda x: self.quantity_at_tau_val('xc3', x),
        'bottom': lambda x: None
    }

    zero_point = kinds[kind](tau_val)
    return zero_point

  def set_z_zero_point(self, zero_point=None):
    # Set the 'z' zero point
    if not zero_point:
      # Use tau=1 as standard
      zero_point = self.z_zero_point()

    self.z = self['z'] - zero_point

  def set_gravity(self, gravity: float):
    self.gravity = gravity

  def print_properties(self):
    # Simple function to print properties of the instance as well as
    # available keys for plotting for this file
    print("=================")
    print("-----------------")
    print("UIOPlot instance properties:")
    print("-----------------")
    print(f"Model ID:\t{self.id}")
    print(f"Located at:\t{self.model_path}")
    print(f"Model number:\t{self.model_num}")
    print(f"Current snapshot:\t{self.snap_idx} of {self.final_snap_idx}")
    print("-----------------")
    print(f"Available keys:")
    print('\n'.join(self.box_keys))
    print("-----------------")

    print("=================")

  def get_key_name_units(self, key: str):
    # For a given key, get the quantity's name and unit as stored in UIO
    # need to add parser like in 'init'
    name = self.box[key].params['n'].split(" ")[-1]
    unit = self.box[key].params['u']

    return name, unit

  def get_min_max_quantity(self, key):
    min_quantity, max_quantity = np.min(
        self[key]), np.max(self[key])

    return min_quantity, max_quantity

  def get_dataset_from_snapshot(self, snap_idx: int):
    return self.full.dataset[snap_idx]

  def get_p_t_profile(self, axis_mean=None):
    # Calculate and return pressure-temperature profile, averaging over axis
    # if specified
    pressure, temperature = self['pressure'], self['temperature']
    if axis_mean:
      pressure = np.mean(pressure.squeeze(), axis=axis_mean)
      temperature = np.mean(temperature.squeeze(), axis=axis_mean)

    return pressure, temperature

  # -------------------------------------------------------------------------
  # Methods for iterating through snapshots
  # -------------------------------------------------------------------------

  def first_snapshot(self):
    self.update_snapshot(0)

  def final_snapshot(self):
    self.update_snapshot(self.final_snap_idx)

  def prev_snapshot(self):
    # Look for previous model and make it the current model
    potential_idx = self.snap_idx - 1
    if potential_idx < self.first_snap_idx:
      potential_idx = self.first_snap_idx

    self.update_snapshot(potential_idx)

  def next_snapshot(self):
    # Look for next model and make it the current model
    potential_idx = self.snap_idx + 1
    if potential_idx >= self.final_snap_idx:
      potential_idx = self.final_snap_idx

    self.update_snapshot(potential_idx)

  def update_snapshot(self, snap_idx: int, box_idx=2):
    # Update 'snap_idx', 'dataset' and 'box' properties
    self.snap_idx = snap_idx
    self.dataset = self.full.dataset[self.snap_idx]
    self.box = self.dataset.box[0]
    self.box_keys = [key for key in self.box.keys()]

    # Mean box
    self.mean_box = self.mean.dataset[self.snap_idx].box[box_idx]
    self.mean_box_keys = [key for key in self.mean_box.keys()]

    # Update quantities
    self.initialise_quantities()

  # -----------------------------------------------------------------------
  # Methods that loop over all snapshots to calculate quantities
  # -----------------------------------------------------------------------

  def average_quantity_over_snapshots(self, key: str, set_min_max=False):
    # Average the specified 'key' from 'self.box' across all snapshots
    # Has kwarg for setting the 'min-max' dict instance for this key so
    # that we don't have to iterate over all the snapshots to do this again
    # (useful for average quantity analyses)
    num_snapshots = self.final_snap_idx + 1
    snap_idx = self.snap_idx  # store reference before iterating

    self.first_snapshot()  # load first snapshot
    avg_quantity = self.box[key].data

    # Calculate min and max quantities as well
    if set_min_max:
      model_min, model_max = np.min(
          self.box[key].data), np.max(self.box[key].data)

    # Loop over snapshots
    for i in range(num_snapshots):
      data = self.box[key].data
      avg_quantity += data

      if set_min_max:
        min_data, max_data = np.min(data), np.max(data)
        if min_data < model_min:
          model_min = min_data
        if max_data > model_max:
          model_max = max_data

      self.next_snapshot()

    avg_quantity /= num_snapshots

    self.update_snapshot(snap_idx)  # revert to original snapshot

    if set_min_max:
      self.min_max_dict[key] = (model_min, model_max)

    return avg_quantity

  def get_quantities_over_snapshots(self, keys: List[str]):
    # Create a list of quantity 'key' across all snapshots
    num_snapshots = self.final_snap_idx + 1
    snap_idx = self.snap_idx  # store reference before iterating
    self.first_snapshot()

    output_quantities = []
    for i in range(num_snapshots):
      if len(keys) == 1:
        quantities = self[keys[0]]
      else:
        quantities = [self[key] for key in keys]
      output_quantities.append(quantities)

      self.next_snapshot()

    self.update_snapshot(snap_idx)  # revert to original snapshot

    return np.array(output_quantities)

  def min_max_quantity_over_snapshots(self, key: str):
    # Get min & max of specified 'key' from 'self.box' across all snapshots
    snap_idx = self.snap_idx  # store reference before iterating

    self.first_snapshot()  # load first snapshot
    model_min, model_max = self.get_min_max_quantity(key)

    # Loop over snapshots
    for i in range(self.final_snap_idx + 1):
      data = self[key]
      min_data, max_data = np.min(data), np.max(data)
      if min_data < model_min:
        model_min = min_data
      if max_data > model_max:
        model_max = max_data

      self.next_snapshot()

    self.update_snapshot(snap_idx)  # revert to original snapshot

    # Set dict and return? Is this confusing because of the side effect?
    self.min_max_dict[key] = (model_min, model_max)
    return model_min, model_max

  # -----------------------------------------------------------------------
  # Data manipulation functions
  # -----------------------------------------------------------------------
  # Change this function to a dict that gets derived quantities, ignore the
  # standard 'else' case because that's obvious
  # Just use the dict to call the correct functions, in the functions, get
  # the quantities necessary instead of passing them in

  # Package below into a separate UIO script (or keep this entire thing
  # as 'uio_utilities')
  # -----------------------------------------------------------------------
  # Convenience functions for box data
  # -----------------------------------------------------------------------

  def get_x_vectors(self, box=None, convert_km=True) -> Union[
          np.ndarray, np.array, np.array]:
    if not box:
      box = self.box
    # Get 'x', 'y', 'z' from box and squeeze empty dimensions
    x = box['xc1'].data.squeeze()
    y = box['xc2'].data.squeeze()
    z = box['xc3'].data.squeeze()

    if convert_km:
      x /= 1e5
      y /= 1e5
      z /= 1e5

    return x, y, z

  def get_velocity_vectors(self, box=None, convert_km=True) -> Union[
          np.ndarray, np.ndarray, np.ndarray]:
    if not box:
      box = self.box
    # Note 3D array indexing is reversed because of Python vs IDL
    v1 = box['v1'].data
    v2 = box['v2'].data
    v3 = box['v3'].data

    if convert_km:
      v1 /= 1e5
      v2 /= 1e5
      v3 /= 1e5

    return v1, v2, v3

  def get_time(self, dataset=None):
    if not dataset:
      dataset = self.dataset
    return dataset['modeltime'].data

  def get_final_time(self, dataset=None):
    snap_idx = self.snap_idx
    self.final_snapshot()
    final_time = self.get_time(dataset=dataset)
    self.update_snapshot(snap_idx)

    return final_time

  def get_time_difference(self, snap_idx1: int, snap_idx2: int):
    # Calculate the time difference between two snapshots
    time1 = self.get_dataset_from_snapshot(snap_idx1)['modeltime'].data
    time2 = self.get_dataset_from_snapshot(snap_idx2)['modeltime'].data

    return time2 - time1

  # -----------------------------------------------------------------------
  # Plotting methods
  # -----------------------------------------------------------------------

  # This should be a separate function in a plotting module
  def plot_heatmap(self, ax, plot_values, log_quantity=False, title=None,
                   plot_type='image', add_cbar=True, cbar=None,
                   cmap='jet', cbar_label='infer', cbar_label_pos='bottom',
                   origin='lower',
                   vlimits=None):
    # Normalise data to range [0, 1] before applying colours
    if log_quantity:
      plot_values = np.log10(plot_values)

    if vlimits:  # should be a tuple of (vmin, vmax)
      vmin, vmax = vlimits
    else:  # determine min,max from plot values
      vmin, vmax = plot_values.min(), plot_values.max()

    norm = Normalize(vmin=vmin, vmax=vmax)

    # Plot a heatmap using X_, Y_grids and a Z datacube 'plot_values'
    if plot_type == 'mesh':
      im = ax.pcolormesh(self.X_grid, self.Y_grid,
                         plot_values, cmap=cmap, norm=norm)

    elif plot_type == 'contour':
      im = ax.contourf(self.X_grid, self.Y_grid,
                       plot_values, cmap=cmap, norm=norm)

    elif plot_type == 'image':
      x_limits = [self.X_grid[0][0], self.X_grid[0][-1]]
      y_limits = [self.Y_grid[0][0], self.Y_grid[-1][0]]
      extent = (x_limits[0], x_limits[1], y_limits[0], y_limits[1])

      im = ax.imshow(plot_values, interpolation='bilinear', origin=origin,
                     cmap=cmap, norm=norm, extent=extent)

    else:
      print(f"Warning: Plot type {plot_type} is not valid. Valid choices are \
        'mesh', 'contour' and 'image'.")
      print("Defaulting to 'mesh'.")
      im = ax.pcolormesh(self.X_grid, self.Y_grid, plot_values, cmap='jet')

    # Set number of ticks
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_locator(MaxNLocator(5))

    if add_cbar:
      # colorbar same height as heatmap
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)

      # Create new colorbar
      if cbar is None:
        cbar = ax.figure.colorbar(im, ax=ax, cax=cax)

      # Set default colorbar ticks
      ticks = np.linspace(vmin, vmax, num=5)
      cbar.set_ticks(ticks)
      cbar.ax.set_yticklabels([f"{t:.2f}" for t in ticks])

      if cbar_label:
        cbar_label_positions = {
            'top': cbar.ax.set_title,
            'right': cbar.set_label,
            'bottom': cbar.ax.set_xlabel,
        }

        if cbar_label_pos in cbar_label_positions.keys():
          cbar_label_positions[cbar_label_pos](cbar_label)
        else:
          print(f"Error: {cbar_label_pos} is not a valid choice."
                "'top', 'right' and 'bottom' are valid choices. Using"
                "'right'")
          cbar_label_positions['right'](cbar_label)

      if title:
        ax.set_title(title)

      return cbar

  def plot_quantity_heatmap(self, ax: plt.Axes, key: str,
                            plane='xz', set_z_zero_to_tau_1=False,
                            plot_tau=False,
                            title=None, cmap='jet', log_quantity=False,
                            xlabel=None, ylabel=None,
                            auto_label_axes=False,
                            origin='lower',
                            add_cbar=True, cbar_label='infer',
                            cbar_label_pos='right',
                            average_snaps=False):
    # Plot quantity in a certain plane as a heatmap. Axis labels are
    # determined by 'plane'

    # Set 'min-max' bounds for a key if they have not yet been set
    set_min_max = False if key in self.min_max_dict else True

    # Get data from box and set 'min-max' bounds for the key
    # Average over all snapshots
    if average_snaps:
      data = self.average_quantity_over_snapshots(
          key, set_min_max=set_min_max)
    else:
      data = self[key]
      if set_min_max:
        self.min_max_quantity_over_snapshots(key)

    # Average over axis not in plane
    if data is not None:
      if len(data.shape) == 3:
        # For now, assume IDL indexing since 'data' has not been converted to
        # Python indexing
        print(f"Averaging data with key {key} in {plane} plane.")
        avg_data = average_data_to_plane(data, plane, is_idl_idx=True)

      else:
        print(
            f"Data with key {key} has shape {data.shape}, skipping plane average.")
        avg_data = data

    else:
      print(f"Error: Could not compute data for key {key}.")
      return None

    # Plot heatmap of data
    # Make grid points for plotting plane
    if set_z_zero_to_tau_1:
      print("Setting z zero point")
      self.set_z_zero_point()

    # Problem: With 'extent', the 'tau' axis gets squashed since it only
    # goes from around -5 to 5, while 'x' will go up to like 5000
    # Plot 'z'  as optical depth
    z = self.tau if plot_tau else self.z

    self.set_grid(*initialise_grid(self.x, self.y, z, plane))

    # if cbar_label == 'infer':
    #   name, unit = self.get_key_name_units(key)
    #   cbar_label = name

    cbar = self.plot_heatmap(ax, avg_data.squeeze(), log_quantity=log_quantity,
                             plot_type='image', origin=origin,
                             title=title, cmap=cmap, add_cbar=add_cbar,
                             cbar_label=cbar_label,
                             cbar_label_pos=cbar_label_pos)

    # Set labels
    if auto_label_axes:
      xlabel, ylabel = [f"{item} [km]" for item in list(plane)]

    if xlabel is not None:
      ax.set_xlabel(xlabel)

    if ylabel is not None:
      ax.set_ylabel(ylabel)

    return cbar

  def plot_keys(self, ax: plt.Axes, key1: str, key2: str,
                title=None, average_snaps=False,
                label=None, colour=None, ls='-',
                xlabel=None, ylabel=None,
                log_key_1=False, log_key_2=False,
                xscale=None, yscale=None, as_abu=False):
    # Plot 'key1' vs 'key2' in the z-direction (no implementation for x,y)
    self.set_mean_box(2)  # set box to 'z3' (z-axis)

    data_list = []
    log_keys = [log_key_1, log_key_2]
    # Set 'min-max' bounds for a key if they have not yet been set
    for key, log_key in zip([key1, key2], log_keys):
      set_min_max = False if key in self.min_max_dict else True

      if key.lower() == 'tau':
        data = np.log10(self.tau)

      elif key.lower() == 'z':
        data = self.z

      elif not key.endswith('_xmean'):
        xkey = key + '_xmean'

      else:
        xkey = key

      # Get data from box and set 'min-max' bounds for the key
      # Average over all snapshots
      if not key.lower() == 'tau' and not key.lower() == 'z':
        if average_snaps:
          data = self.average_quantity_over_snapshots(
              xkey, set_min_max=set_min_max)
        else:
          data = self[key]

        if set_min_max:
          if key in ['temperature', 'pressure', 'entropy']:
            self.min_max_quantity_over_snapshots(key)
          else:
            self.min_max_quantity_over_snapshots(xkey)

        # Convert QUC to abundances if specified
        # Assumes QUC001 is the hydrogen number density
        if key.startswith('quc') and as_abu:
          data = number_density_to_abundance(data, self['quc001_xmean'])

      if log_key:
        data = np.log10(data)

      data_list.append(data.squeeze())

    # Plot key1 data vs key2 data
    if colour:
      ax.plot(data_list[0], data_list[1], c=colour, label=label, ls=ls)
    else:
      ax.plot(data_list[0], data_list[1], label=label, ls=ls)

    if title:
      ax.set_title(title)

    if xlabel:
      ax.set_xlabel(xlabel)
    if ylabel:
      ax.set_ylabel(ylabel)

    if xscale:
      ax.set_xscale(xscale)
    if yscale:
      ax.set_yscale(yscale)

  def z_plot(self, ax: plt.Axes, key: str,
             title=None, average_snaps=False,
             label=None, colour=None, ls='-',
             xlabel=None, ylabel=None,
             xscale='linear', yscale='linear',
             as_abu=False):
    # Convenience method for plotting keys against 'z'; just light wrapper
    # around 'self.plot_keys()'
    if not xlabel:
      xlabel = 'z [km]'
    self.plot_keys(ax, key1='z', key2=key, title=title,
                   average_snaps=average_snaps,
                   label=label, colour=colour, ls=ls,
                   xlabel=xlabel, ylabel=ylabel,
                   xscale=xscale, yscale=yscale, as_abu=as_abu)

  def tau_plot(self, ax: plt.Axes, key: str,
               title=None, average_snaps=False,
               label=None, colour=None, ls='-',
               xlabel=None, ylabel=None,
               xscale='linear', yscale='linear',
               log_quantity=False,
               as_abu=False):
    # Convenience method for plotting keys against 'log_tau'; just light
    # wrapper around 'self.plot_keys()'
    if not ylabel:
      ylabel = r'$\log{\tau}$'

    self.plot_keys(ax, key1=key, key2='tau', title=title,
                   average_snaps=average_snaps,
                   log_key_1=log_quantity, log_key_2=False,
                   label=label, colour=colour, ls=ls,
                   xlabel=xlabel, ylabel=ylabel,
                   xscale=xscale, yscale=yscale, as_abu=as_abu)
