# Numerical verification record

## Homogeneous ceramic plate

The classical Ritz solution matches the homogeneous simply supported Kirchhoff–Love benchmark: **594.1807407444 Hz**.

| Abaqus S4R mesh | First transverse-bending frequency | Mode number |
|---|---:|---:|
| 20 × 20 | 578.76 Hz | 2 |
| 40 × 40 | 576.98 Hz | 2 |
| 80 × 80 | 576.50 Hz | 2 |

The 40→80 refinement changes this frequency by approximately −0.083%.

## Thickness comparison, 40 × 40 mesh

| Thickness (m) | Kirchhoff reference (Hz) | Recorded S4R (Hz) | Bending mode number |
|---|---:|---:|---:|
| 0.04 | 594.18 | 576.98 | 2 |
| 0.02 | 297.09 | 293.64 | 1 |
| 0.01 | 148.55 | 147.94 | 1 |

The deviations decrease approximately 2.89% → 1.16% → 0.41% as the plate gets thinner. This is consistent with a model-form difference, but does not isolate shear and rotary-inertia contributions.

## Graded plate benchmarks

| Gradation (α, β) | Ritz 10 × 10 (Hz) | Recorded Abaqus S4R, 40 × 40 (Hz) |
|---|---:|---:|
| (0,0) | 594.18074 | 576.98 |
| (1,1) | 402.57036 | 391.20 |
| (5,5) | 316.53602 | 305.98 |
| (1,5) | 334.65096 | 323.58 |
| (5,1) | 334.65096 | 323.58 |

The equal (1,5) and (5,1) frequencies provide a coordinate-exchange symmetry check for the square plate.

## Nodal U3 shape verification

| Gradation | Single-sine MAC | Full Ritz 10 × 10 MAC |
|---|---:|---:|
| (0,0) | 0.99998717 | 0.99998717 |
| (1,1) | 0.98689444 | 0.99998642 |
| (1,5) | 0.99391837 | 0.99997379 |

The tabulated MAC values are calculated from the committed Abaqus U3 displacement CSVs and linear Ritz eigenvectors. Reproduce them using `python scripts/compare_fem.py`; the published datasets are evaluated by the automated regression tests.

## Scope of automated checks

Automated tests cover material gradation, homogeneous reference frequencies, graded Ritz benchmarks, coordinate exchange, thickness scaling, recorded U3 modal correlation, MAC algebra, and input-deck generation. Abaqus solver execution takes place separately in a licensed environment.
