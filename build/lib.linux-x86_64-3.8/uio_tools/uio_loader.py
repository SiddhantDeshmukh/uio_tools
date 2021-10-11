from uio_tools.uio_reader import UIOData
import numpy as np
import glob
import os
from typing import List


class UIOLoader():
  # List holding filepaths for UIOPlot() instances
  # Contains helper methods to load next, previous models and holds state
  # of plotting variables to be the same across all models
  def __init__(self, model_directory: str, eos_file='', opta_file='') -> None:
    self.model_dir = model_directory
    self.model_files = None
    self.num_models = None
    self.id = ''
    self.idx = 0  # current idx in 'model_files'
    self.current_model = None
    self.current_model_path = None
    self.eos_file = eos_file if eos_file else None
    self.opta_file = opta_file if opta_file else None

    self.load_model_files()

  def load_model(self):
    # Load the current model
    self.current_model_path = f"{self.model_files[self.idx]}"
    self.current_model = UIOData(
        self.current_model_path, self.eos_file, self.opta_file)

  def load_first_model(self):
    self.idx = 0
    self.load_model()

  def load_final_model(self):
    self.idx = len(self.model_files) - 1
    self.load_model()

  def load_next_model(self):
    new_idx = self.idx + 1
    # Could return something if we're at the first/last model to auto-break
    # out of external loops?
    if new_idx >= len(self.model_files):
      print("Already at last model. Reloading model.")
    else:
      self.idx = new_idx

    self.load_model()

  def load_prev_model(self):
    new_idx = self.idx - 1
    if new_idx <= 0:
      print("Already at first model. Reloading model.")
    else:
      self.idx = new_idx

    self.load_model()

  def load_model_files(self):
    # Parse files for identifier string, earliest and latest model
    # Assumes all model files have structure "identifier.YYY.full"
    # where YYY is a digit, often 00 -> 79
    # Use first file to determine identifier
    # Assumes there is only one 'identifier' in the directory, which is an
    # obvious problem!
    # Store references to all 'full' files; UIOData will automatically load
    # Full and Mean Data
    files = glob.glob(os.path.join(self.model_dir, f"*.full"))
    if not files:
      print(f"Error: No models found in {self.model_dir}")
      return None

    identifier = files[0].split('/')[-1].split('.')[0]
    filtered_files = [file_ for file_ in files if identifier in file_]
    filtered_files.sort()

    # Save filtered files
    self.model_files = filtered_files
    self.num_models = len(self.model_files)
    self.id = identifier

    print(f"ID: {self.id}. {len(filtered_files)} files found.")

    # Set current model to first model
    self.idx = 0
    self.load_model()

  def print_files(self):
    print(f"{self.model_dir} contains:")
    print('\n'.join([f"\tIndex {i}: {file_}" for i,
                     file_ in enumerate(self.model_files)]))

  # -------------------------------------------------------------------------
  # Functions that loop over models
  # -------------------------------------------------------------------------
  def get_quantities_over_models(self, keys: List[str]):
    # Create a list of quantity 'key' across all models
    num_snapshots = len(self.model_files)
    model_idx = self.idx  # store reference before iterating
    self.load_first_model()

    quantities = {}
    for i in range(num_snapshots):
      data = self.current_model.get_quantities_over_snapshots(keys,
                                                              as_arrays=False)
      for key in data.keys():
        if not key in quantities:
          quantities[key] = data[key]
        else:
          # Append on 'snaps' axis
          quantities[key] += (data[key])

      self.load_next_model()

    self.idx = model_idx
    self.load_model()  # revert to original snapshot

    # Convert to arrays
    quantities = {key: np.array(val) for key, val in quantities.items()}

    return quantities

  def average_quantities_over_models(self, keys: List[str]):
    # Get specified quantities across all models and then compute the average
    # in time
    quantities = self.get_quantities_over_models(keys)
    quantities = np.mean(quantities, axis=0)

    return quantities
