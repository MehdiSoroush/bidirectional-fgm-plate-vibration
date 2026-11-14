# Abaqus 2021 linear modal benchmark

The Abaqus scripts use the Python 2.7 runtime distributed with Abaqus 2021. The input-deck generator also runs under modern Python, while the U3 exporter needs the Abaqus ODB modules.

Generate a 40 × 40 mesh graded plate input deck:

    abaqus python abaqus/generate_abaqus2021.py --n 40 --alpha 1 --beta 1 --h 0.04 --outdir inputs

Or generate all five reference gradations:

    abaqus python abaqus/generate_abaqus2021.py --all --n 40 --outdir inputs

Submit the generated .inp file using the licensed Abaqus solver and inspect U3 contours to identify the first true *transverse bending* mode. At the baseline h=0.04 m, the recorded bending mode was Mode 2, not Mode 1.

Export a selected U3 eigenmode to CSV:

    abaqus python abaqus/export_u3_abaqus2021.py Job.odb 2 mode_11_n40.csv

The exported CSV has columns instance,node,x,y,z,u3. The modal-comparison script **scripts/compare_fem.py** reads the displacement datasets in **data/abaqus_u3/**.

**Boundary conditions:** U3=0 along all edges, minimal in-plane pins to remove rigid-body modes. These restraints are intended for the documented *linear* vibration benchmark. They should not be reused unchanged for nonlinear von Kármán membrane modelling.

**Numerical scope:** The generator writes input decks for linear modal analysis. Additional solver results require independent verification, and S4R shell kinematics differ from the Kirchhoff–Love analytical formulation.
