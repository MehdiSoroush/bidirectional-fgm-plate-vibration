"""Plot recorded linear frequency benchmarks from tabulated data.

The plots here require no proprietary Abaqus installation.
The FEM numbers were recorded from earlier Abaqus 2021 runs,
not recalculated by this script.
"""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "reference"
OUTPUT = ROOT / "results" / "figures"


def rows(name):
    with (DATA / name).open(newline="") as file:
        return list(csv.DictReader(file))


def save(name):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plt.grid(alpha=0.20)
    plt.tight_layout()
    plt.savefig(OUTPUT / name, dpi=180)
    plt.close()


def main():
    sweep = rows("symmetric_gradation_sweep.csv")
    p = [float(r["alpha"]) for r in sweep]
    f = [float(r["frequency_hz"]) for r in sweep]
    plt.figure(figsize=(8, 4.8))
    plt.plot(p, f, "o-", linewidth=2)
    plt.xlabel("Symmetric gradation exponent p")
    plt.ylabel("First transverse bending frequency (Hz)")
    plt.title("Recorded Abaqus S4R graded-plate sweep, 40 x 40")
    save("frequency_vs_gradation.png")

    mesh = rows("homogeneous_mesh.csv")
    n = [int(r["mesh_n"]) for r in mesh]
    frequencies = [float(r["frequency_hz"]) for r in mesh]
    plt.figure(figsize=(8, 4.8))
    plt.plot(n, frequencies, "o-", linewidth=2)
    plt.xlabel("Elements per plate edge")
    plt.ylabel("First transverse bending frequency (Hz)")
    plt.title("Recorded Abaqus mesh refinement: homogeneous ceramic")
    save("homogeneous_mesh_convergence.png")

    thickness = rows("thickness_study.csv")
    h = [1000 * float(r["thickness_m"]) for r in thickness]
    classical = [float(r["ritz_cpt_hz"]) for r in thickness]
    shell = [float(r["abaqus_s4r_hz"]) for r in thickness]
    plt.figure(figsize=(8, 4.8))
    plt.plot(h, classical, "o-", label="Classical Kirchhoff plate")
    plt.plot(h, shell, "s--", label="Recorded Abaqus S4R")
    plt.xlabel("Thickness (mm)")
    plt.ylabel("First transverse bending frequency (Hz)")
    plt.title("Thickness check")
    plt.legend()
    save("thickness_comparison.png")


if __name__ == "__main__":
    main()
