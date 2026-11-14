"""Reproducible linear plate benchmarks from the FGM vibration research project."""
from .material import Plate, ceramic_fraction, properties
from .ritz import solve, evaluate_shape, exact_homogeneous_frequency
from .modal import load_u3, mac, compare_to_ritz
__all__ = ['Plate', 'ceramic_fraction', 'properties', 'solve', 'evaluate_shape',
           'exact_homogeneous_frequency', 'load_u3', 'mac', 'compare_to_ritz']
