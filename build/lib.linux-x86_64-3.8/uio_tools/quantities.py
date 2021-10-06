# Utility script for calculating properties like kinetic energy and momentum
import numpy as np


def calculate_kinetic_energy(density: np.ndarray, v1: np.ndarray,
                             v2: np.ndarray, v3: np.ndarray) -> np.ndarray:
  kinetic_energy = 0.5 * density * (v1**2 + v2**2 + v3**2)
  return kinetic_energy


def calculate_momentum(density: np.ndarray, v1: np.ndarray,
                       v2: np.ndarray, v3: np.ndarray) -> np.ndarray:
  momentum = density * np.sqrt(v1**2 + v2**2 + v3**2)
  return momentum

# ------------------------------------------------------------------------------
# Hydrodynamics quantities
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Thermodynamic quantities
# ------------------------------------------------------------------------------
