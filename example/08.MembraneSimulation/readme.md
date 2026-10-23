This example shows how to add an implicit lipid force field to simulate a protein

normal membrane protein simulation:

    (source ../../source.sh)
    python 0.normal.run.py

channel membrane protein simulation:

    python 1.channel.run.py

a membrane of fixed curvature (radius 120 A):

    python 3.fixed_curvature.run.py

a curvature that evolves with the proteins (`--curvature-changer-interval`), three 1rkl copies held by nails
(nail-xyz.dat, radius 250 A to start; nail-y.dat, radius 140 A to start):

    python 4.curvature_dynamics1.run.py
    python 5.curvature_dynamics2.run.py

TODO:
show the identification of protein-lipid and protein-water interfaces

