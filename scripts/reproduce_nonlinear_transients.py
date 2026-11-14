"""Reconstruct the undamped cubic responses and free-frequency curve."""
from pathlib import Path
import csv
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from fgm_plate.analytical_single_mode import source_coefficients
from fgm_plate.nonlinear import DuffingModel, exact_free_angular_frequency

DATA=ROOT/"data/reference/nonlinear_source_coefficients.csv"
OUT=ROOT/"results/figures"


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    with DATA.open(newline="") as f:
        ref={r["case"]:r for r in csv.DictReader(f)}
    fig,ax=plt.subplots(figsize=(8.5,4.6))
    for name,label in [
        ("s00",r"$(\alpha,\beta)=(0,0)$"),
        ("s10",r"$(\alpha,\beta)=(1,0)$"),
        ("s11",r"$(\alpha,\beta)=(1,1)$"),
    ]:
        row=ref[name]
        c=source_coefficients(float(row["alpha"]),float(row["beta"]),
            slenderness=float(row["slenderness"]),
            aspect_ratio=float(row["aspect_ratio"]))
        model=DuffingModel(c.p,c.q,forcing_amplitude=5e6,forcing_frequency=1500)
        t,x,_=model.integrate(t_end=.008,q0=.02,v0=0,num=1601)
        ax.plot(t*1e3,x,label=label,lw=1.65)
    ax.set(xlabel="Time (ms)",ylabel="Generalized displacement q (m)",
           title="Reconstructed nonlinear transient responses")
    ax.grid(alpha=.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT/"nonlinear_transient_reconstruction.png",dpi=200)
    plt.close(fig)

    c=source_coefficients(0,0)
    amplitudes=np.linspace(0,.02,81)
    normalized=np.array([
        exact_free_angular_frequency(c.p,c.q,float(a))/np.sqrt(c.p)
        for a in amplitudes
    ])
    fig,ax=plt.subplots(figsize=(8.2,4.5))
    ax.plot(amplitudes*1e3,normalized,lw=1.8)
    ax.set(xlabel="Free-vibration amplitude (mm)",
           ylabel=r"Frequency ratio $\omega(A)/\omega_0$",
           title="Source-model cubic hardening at fixed square-plate geometry")
    ax.text(.02,.94,
        r"Gradations (0,0), (1,0), (1,1) coincide: $Q/P=822.65625$ m$^{-2}$",
        transform=ax.transAxes,va="top",fontsize=9.8)
    ax.grid(alpha=.25)
    fig.tight_layout()
    fig.savefig(OUT/"nonlinear_free_frequency.png",dpi=200)
    plt.close(fig)

    print("Reconstruction figures generated in",OUT)
    print("Single-mode nonlinear oscillator responses computed.")


if __name__=="__main__":
    main()
