# While this script is meant for quick prototyping, it currently contains
# a whole host of data-viewing and plotting routines that need to be moved
# to other files. Please excuse the large commented blocks!

# %%
import matplotlib.pyplot as plt
import numpy as np
from uio_tools import eosinter
from uio_tools.uio_utils import UIOLoader, initialise_grid
from uio_tools.uio_plot_tools import plot_initial_final_heatmap, plot_initial_final_z_trend

# print('\n'.join(plt.style.available))
plt.style.use('standard-scientific')

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


# BELOW ARE PLOTTING ROUTINES THAT NEED TO GO INTO A DIFFERENT SCRIPT
# %%
# model directories
# solar-type
solar_model_dir = r"../res/cobold-runs/solar"
solar_sven_dir = f"{solar_model_dir}/2d-wedemeyer"
solar_umist12_dir = f"{solar_model_dir}/2d-umist12"
mp_model_dir = r"../res/cobold-runs/metalpoor"
mp_sven_dir = f"{mp_model_dir}/2d-wedemeyer"
mp_solarHe_dir = f"{mp_model_dir}/2d-solarHe"
mp_umist12_dir = f"{mp_model_dir}/2d-umist12"

# giant
giant_mm00_model_dir = r"../res/cobold-runs/d2t50g25mm00"
giant_mm00cu12_model_dir = r"../res/cobold-runs/d2t50g25mm00cu12"
giant_mm20_model_dir = r"../res/cobold-runs/d2t50g25mm20"
giant_mm20cu12_model_dir = r"../res/cobold-runs/d2t50g25mm20cu12"


model_dirs = [solar_sven_dir, solar_umist12_dir,
              mp_sven_dir, mp_umist12_dir,
              giant_mm00cu12_model_dir,
              giant_mm20cu12_model_dir]

# eos_dir = '/home/sdeshmukh/Documents/chemicalAnalysis/res/eos'
# eos_files = ['eos_mm00_l5.eos', 'eos_mm20_l.eos']
# eos_files = [f"{eos_dir}/{file_}" for file_ in eos_files]

quc_keys = ['quc00' + str(i+1) + '_xmean' for i in range(0, 8)]

# %%
# # How do the chemical quantities vary with abundance and height? (mean)
# mean_fig_dir = '../figs/mean-figs'
# fignames = ['solar-variation.pdf', 'metalpoor-variation.pdf']
# fignames = [f"{mean_fig_dir}/{name}" for name in fignames]

# for (model_dir, figname, eos_file) in zip(model_dirs, fignames, eos_files):
#   fig, axes = plt.subplots(1, 2)

#   # UIOLoader
#   loader = UIOLoader(model_dir, file_type='mean', eos_file=eos_file)
#   loader.load_final_model()
#   model = loader.current_model

#   ylabel = r"Number Density [cm$^{-3}$]"

#   for key in quc_keys:
#     species, unit = model.get_key_name_units(key)

#     # Height
#     model.z_plot(axes[0], key, xscale='log',
#                  yscale='log', line_label=species, ylabel=ylabel)
#     axes[0].legend()

#     # Temperature
#     model.plot_keys(axes[1], 'temperature', key, line_label=species,
#                     xlabel='Temperature [K]', ylabel=ylabel,
#                     xscale='log', yscale='log')
#     axes[1].legend()

#   plt.savefig(figname, bbox_inches="tight")
#   plt.close('all')


# %%
# Plot multiple models in multiple figs
fignames = ['solar-co-ch-wedemeyer.pdf',
            'solar-co-ch-umist12.pdf',
            'metalpoor-co-ch-wedemeyer.pdf',
            'metalpoor-co-ch-umist12.pdf',
            'd2t50g25mm00cu12.pdf',
            'd2t50g25mm20cu12.pdf']
keys = ['quc005', 'quc006']  # CO, CH

# Initial-final heatmaps
# One figure for each model
for (model_dir, figname) in zip(model_dirs, fignames):
  # keys on rows, columns: (initial, final)
  fig, axes = plt.subplots(len(keys), 2, figsize=(10, 6), sharey=True)
  for i, key in enumerate(keys):
    title_prefix = "Time = "
    plot_initial_final_heatmap(
        model_dir, key, axes[i][0], axes[i][1], title_prefix=title_prefix)

  plt.savefig(f"../writeup/2d-evolution/{figname}", bbox_inches="tight")
  plt.close("all")

# Initial-final mean figs for all QUC
quc_keys = ['quc00' + str(i+1) + '_xmean' for i in range(0, 8)]
fignames = ['mean-' + name for name in fignames]
for (model_dir, figname) in zip(model_dirs, fignames):
  # keys on rows, columns: (initial, final)
  fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
  title_prefix = "Time = "
  plot_initial_final_z_trend(
      model_dir, quc_keys, axes[0], axes[1], title_prefix=title_prefix)

  plt.savefig(f"../writeup/2d-evolution/{figname}", bbox_inches="tight")
  plt.close("all")

# %%
# # Check derived quantities for models
# fignames = ['solar-quantities.pdf',
#             'metalpoor-quantities.pdf',
#             'metalpoor-solarHe-quantities.pdf']

# eos_dir = '/home/sdeshmukh/Documents/chemicalAnalysis/res/eos'
# eos_files = ['eos_mm00_l5.eos', 'eos_mm20_l.eos', 'eos_mm20_l.eos']
# eos_files = [f"{eos_dir}/{file_}" for file_ in eos_files]

# for (model_dir, figname, eos_file) in zip(model_dirs, fignames, eos_files):
#   fig, axes = plt.subplots(2, 2)

#   # Load model
#   loader = UIOLoader(model_dir, file_type='full', eos_file=eos_file)
#   loader.load_final_model()
#   model = loader.current_model

#   # Plot density, internal energy
#   model.plot_quantity_heatmap(axes[0][0], 'rho', title="Density")
#   model.plot_quantity_heatmap(axes[0][1], 'ei', title="Internal Energy")

#   # Plot pressure, temperature (EOS quantities)
#   model.plot_quantity_heatmap(axes[1][0], 'pressure', title="Pressure")
#   model.plot_quantity_heatmap(
#       axes[1][1], 'temperature', title="Temperature")

#   plt.savefig(figname, bbox_inches="tight")
#   plt.close('all')

# # %%
# # Plot multiple models in same fig
# # Initial-final heatmaps
# key = 'quc005'  # CO

# fig, axes = plt.subplots(3, 2, figsize=(12, 12))

# for (model_dir, figname, ax) in zip(model_dirs, fignames, axes):
#   plot_initial_final_heatmap(model_dir, key, ax_initial=ax[0], ax_final=ax[1],
#                              title_prefix=f"{figname.replace('.pdf', '')}:")

# plt.savefig('../writeup/2d-evolution/full-co-comparison.pdf',
#             bbox_inches="tight")
# plt.close('all')

# # Initial-final z-plots (mean files)
# keys = ['quc00' + str(i+1) + '_xmean' for i in range(0, 8)]
# fignames = ['mean-' + name for name in fignames]

# fig, axes = plt.subplots(3, 2, figsize=(14, 14))

# for (model_dir, figname, ax) in zip(model_dirs, fignames, axes):
#   plot_initial_final_z_trend(
#       model_dir, keys, ax_initial=ax[0], ax_final=ax[1],
#       title_prefix=f"{figname.replace('.pdf', '')}:")

# plt.savefig('../writeup/2d-evolution/mean-comparison.pdf',
#             bbox_inches="tight")

# plt.close('all')


# # %%
# # Test z-t heatmaps
# # Should be a method of UIOLoader like the other two plots, go through
# # every mean snapshot in chronological order and plot the evolution as a
# # heatmap?
# loader = UIOLoader(solar_model_dir, file_type='mean')
# loader.load_first_model()
# model = loader.current_model

# # Iterate through snapshots, get data and info for grid
# times = []
# cos = []
# chs = []
# z = model.get_z()

# for i in range(model.final_snap_idx - 1):
#   snap_time = model.dataset['modeltime']
#   co = model.get_box_quantity('quc005_xmean')
#   ch = model.get_box_quantity('quc006_xmean')

#   times.append(snap_time)
#   cos.append(co)
#   chs.append(ch)

#   model.next_snapshot()

# fig, axes = plt.subplots(1, 2)

# times = np.array(times)
# X_grid, Y_grid = initialise_grid(z, times, None, 'xy')

# # CO
# im = axes[0].imshow(cos, interpolation='bilinear', origin='lower',
#                     cmap='jet')
# cbar = axes[0].figure.colorbar(im, ax=axes[0])


# # CH
# im = axes[1].imshow(chs, interpolation='bilinear', origin='lower',
#                     cmap='jet')
# cbar = axes[1].figure.colorbar(im, ax=axes[1])

# plt.show()

# # %%
# # Subplots for QUC (CHEM)
# nrows = 4  # 3
# ncols = 2  # 3

# max_idx = 8  # num species, aka num QUC arrays

# # Plot Sven's 2D initial data and then the latest 2D snapshot
# # Colorbar limits have to be the same in both cases for each QUC array
# # plot_type = 'chem-single'  # a CHEM run with 'plustimestep = 1'
# # plot_type = 'chem-initial'  # initial CHEM model with correct abu init
# # plot_type = 'chem-test'  # CHEM run with 'plustimestep=20' to test OpenMP
# # plot_type = '2d-initial'  # 2D solar model snap 00
# # plot_type = '2d-sven'  # Sven's 2D model in x-z plane (y sliced out)
# # plot_type = '2d-final'  # 2D model in x-z plane (y sliced out)
# # plot_type = '2d-relax-initial'  # 2D model in x-z plane
# # plot_type = '2d-relax-final'  # 2D model in x-z plane

# plane = 'xz'

# # Initial plot, set colorbar limits to (min, max) of each array
# cbar_limits = []
# species_names = []

# # uio_single_2d.load_first_model()
# uio_single_2d.load_final_model()
# model = uio_single_2d.current_model
# # plot_type = f"{uio_single_2d.id}-{model.model_num}"

# # Loop over snapshots
# for snap_idx in range(model.final_snap_idx):
#   # Just CO figures
#   quc_key = "quc005"
#   species = model.box[quc_key].params['n'].split(" ")[-1]
#   title = f"t = {model.dataset['modeltime'].data:.2f}"

#   fig, axes = plt.subplots(1, 1, figsize=(8, 8))
#   model.plot_quantity_heatmap(axes, quc_key, plane=plane,
#                               title=title, cbar_label=species)

#   # plt.savefig(
#   #     f'../figs/co-evolution/2d/{uio_single_2d.id}-{model.model_num}-{str(snap_idx).zfill(4)}.png', bbox_inches="tight")
#   plt.close('all')
#   print(f"Done {snap_idx + 1} / {model.final_snap_idx} snapshots.")

#   model.next_snapshot()

# # # All QUC
# # fig_chem, axes_chem = plt.subplots(nrows, ncols, figsize=(18, 18))
# # for i in range(nrows):
# #   for j in range(ncols):
# #     quc_idx = (i * ncols) + j + 1
# #     quc_key = f"quc00{quc_idx}"

# #     species = uio_single_2d.box[quc_key].params['n'].split(" ")[-1]
# #     # unit = box.params['u']
# #     # title = rf"n$_{species}$"
# #     title = r"n$_{\mathrm{" + species + r"}}$"

# #     if quc_idx > max_idx:
# #       continue
# #     uio_single_2d.plot_quantity_heatmap(axes_chem[i][j], quc_key,
# #                                         plane=plane, title=title)
# #     # cbar = plot_quantity_heatmap(axes_chem[i][j], box, quc_key,
# #     #                       plane=plane, title=title)

# # figname = f'{plot_type}-{plane}.png'
# # # plt.savefig(f'../writeup/chem-solar/{figname}', bbox_inches="tight")
# # plt.savefig(f'../figs/{figname}', bbox_inches="tight")

# %%
