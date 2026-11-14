"""Shape metrics for matching mode-order-identified Abaqus U3 and Ritz fields."""
import csv
import numpy as np
from .material import Plate
from .ritz import evaluate_shape


def load_u3(path):
    """Read exported Abaqus mode U3; retain node order (one part instance only)."""
    rows = []
    with open(path, newline='', encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            rows.append((r['instance'], int(r['node']), float(r['x']), float(r['y']), float(r['u3'])))
    if not rows:
        raise ValueError('Empty displacement CSV')
    if len({r[0] for r in rows}) != 1:
        raise ValueError('Only single-instance nodal datasets supported; separate instances first')
    if len({r[1] for r in rows}) != len(rows):
        raise ValueError('Duplicate node labels')
    a = np.array([(r[2], r[3], r[4]) for r in rows], dtype=float)
    if not np.isfinite(a).all():
        raise ValueError('Nonfinite nodal values')
    return a[:, 0], a[:, 1], a[:, 2]


def mac(u, v):
    """Unweighted nodal modal assurance criterion (sign invariant)."""
    u = np.asarray(u, dtype=float).ravel()
    v = np.asarray(v, dtype=float).ravel()
    if u.shape != v.shape or not np.isfinite(u).all() or not np.isfinite(v).all():
        raise ValueError('Invalid mode vector shapes/values')
    den = np.dot(u, u)*np.dot(v, v)
    if den <= 0:
        raise ValueError('Zero-norm mode')
    return float((np.dot(u, v)**2)/den)


def shape_residual(u, v):
    """Relative L2 best-fit residual; equals sqrt(1 - MAC) for nonzero data."""
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    if np.dot(v, v) <= 0 or np.dot(u, u) <= 0:
        raise ValueError('Zero-norm mode')
    scale = np.dot(u, v)/np.dot(v, v)
    return float(np.linalg.norm(u-scale*v)/np.linalg.norm(u))


def compare_to_ritz(csv_path, result, plate=Plate()):
    x, y, u = load_u3(csv_path)
    one = np.sin(np.pi*x/plate.a)*np.sin(np.pi*y/plate.b)
    enriched = evaluate_shape(result, x, y, plate)
    return {
        'nodes': len(x),
        'single_sine_mac': mac(u, one),
        'full_ritz_mac': mac(u, enriched),
        'single_sine_residual': shape_residual(u, one),
        'full_ritz_residual': shape_residual(u, enriched),
    }
