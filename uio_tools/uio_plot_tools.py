from typing import List
from uio_utils import FullData, MeanData, UIOLoader
import matplotlib.pyplot as plt

import numpy as np
# These can be methods of the Loader; for the bounds, would need to iterate
# through all instances to determine (vmin,vmax) and pass these in!

# Do I need this function?


def heatmap_comparison(figname: str):
  # Open 3 UIOLoader instances for 3 separate models
  solar_model_dir = r"../res/cobold-runs/solar-2d"
  mp_solarHe_model_dir = r"../res/cobold-runs/metalpoor-2d-solarHe"
  mp_model_dir = r"../res/cobold-runs/metalpoor-2d"

  solar_loader = UIOLoader(solar_model_dir, file_type='full')
  mp_solarHe_loader = UIOLoader(mp_solarHe_model_dir, file_type='full')
  mp_loader = UIOLoader(mp_model_dir, file_type='full')

  loaders = [solar_loader, mp_solarHe_loader, mp_loader]
  titles = ["[Fe/H]=0.0, M=11.00",
            "[Fe/H]=-2.0, M=11.00", "[Fe/H]=-2.0, M=10.93"]

  fig, axes = plt.subplots(2, 3, figsize=(12, 4))

  for i, row_ax in enumerate(axes):
    for (loader, title, ax) in zip(loaders, titles, row_ax):
      # Top row: Initial state
      if i == 0:
        loader.load_first_model()
        model = loader.current_model
        model.first_snapshot()

      # Bottom row: Final state
      elif i == 1:
        loader.load_final_model()
        model = loader.current_model
        model.final_snapshot()

      quc_key = "quc005"  # CO
      species = model.box[quc_key].params['n'].split(" ")[-1]
      # suptitle = f"t = {model.dataset['modeltime'].data:.2f}"

      model.plot_quantity_heatmap(ax, quc_key, plane='xz',
                                  title=title, cbar_label=species)

      # plt.suptitle(suptitle)

  plt.savefig(figname, bbox_inches="tight")


def plot_initial_final_heatmap(model_dir: str, key: str,
                               ax_initial, ax_final,
                               title_prefix=''):
  # UIOLoader
  loader = UIOLoader(model_dir, file_type='full')

  # Initial state
  loader.load_first_model()
  model = loader.current_model
  species, unit = model.get_key_name_units(key)
  model.first_snapshot()
  title = f"{model.dataset['modeltime'].data:.2f} s"
  title = f"{title_prefix} {title}"
  model.plot_quantity_heatmap(
      ax_initial, key, plane='xz', title=title, auto_label_axes=True,
      cbar_label=species, cbar_label_pos='bottom')

  # Final state
  loader.load_final_model()
  model = loader.current_model
  species, unit = model.get_key_name_units(key)
  model.first_snapshot()  # final snapshot time can be slightly different
  title = f"{model.dataset['modeltime'].data:.2f} s"
  title = f"{title_prefix} {title}"
  model.plot_quantity_heatmap(
      ax_final, key, plane='xz', title=title, xlabel='x [km]',
      cbar_label=species, cbar_label_pos='bottom',
      ylabel=None)

  # Tick params
  ax_final.yaxis.set_tick_params(labelleft=True)


def plot_initial_final_z_trend(model_dir: str, keys: List[str],
                               ax_initial, ax_final, title_prefix=''):
  # UIOLoader
  loader = UIOLoader(model_dir, file_type='mean')

  # Initialise figure texts
  ylabel = r"Number Density [cm$^{-3}$]"

  # Initial state
  loader.load_first_model()
  model = loader.current_model
  model.first_snapshot()
  title = f"{model.dataset['modeltime'].data:.2f} s"
  title = f"{title_prefix} {title}"
  for k, key in enumerate(keys):
    if not key.endswith('_xmean'):
      key += '_xmean'
    species, unit = model.get_key_name_units(key)
    model.z_plot(ax_initial, key, title=title,
                 line_label=species, ylabel=ylabel,
                 xscale='log', yscale='log')

  ax_initial.legend(loc='lower left')

  # Final state
  loader.load_final_model()
  model = loader.current_model
  model.first_snapshot()  # final snapshot time can be slightly different
  title = f"{model.dataset['modeltime'].data:.2f} s"
  title = f"{title_prefix} {title}"
  for k, key in enumerate(keys):
    if not key.endswith('_xmean'):
      key += '_xmean'
    species, unit = model.get_key_name_units(key)
    model.z_plot(ax_final, key, title=title,
                 line_label=species, ylabel=None,
                 xscale='log', yscale='log')

    # Tick params
    ax_final.yaxis.set_tick_params(labelleft=True)

  ax_final.legend(loc='lower left')
