Simulation using PHAST of the system described in
   Tang, D.H., E.O. Frind, and E.A. Sudicky. 1981. Contaminant transport in
   fractured porous media: Analytical solution for a single fracture,
   Water Resources Research 17, 555-564.
   https://doi.org/10.1029/WR017i003p00555

The chemistry database is "%INSTALLDIR%\database\phast.dat"
The "real" input files are
   Tang.chem.dat
   Tang.trans.dat
all other files are created by Phast4windows, except
files with names starting with an underscore (_).

As explained in the PHAST manual, nodes define `elements` and `cells`.
An element is a rectangular prism defined by exactly eight nodes that
are located at the element corners.  Every node is contained within
a single cell. Cell faces bisect the distance between adjacent nodes.

Porous-media properties must be defined for each _element_, including
porosity, hydraulic conductivity, specific storage, and dispersivity.
Porous-media properties for a _cell_ are defined by aggregating the
properties of the active elements that are contained in the cell.
When simulating flow and solute transport in a fracture and in the
rock matrix, we have to define nodes so that the cells have either
(a) fracture properties, or (b) matrix properties.

In the case of Tang et al (1981) we use a non-uniform
Z-grid with the following nodes:

   0    5.5e-005    6.5e-005    6.7e-005    0.0002
   0.0006    0.0012    0.002    0.003

This gives the following cell borders:

   0       2.75e-5  6.0e-5   6.6e-5   0.0001335
   0.0004  0.0009   0.0016   0.0025   0.003

The fracture half-width is 60 um in Tang et al (1981).
So the first two cells are within the fracture.
In the input file the MEDIA properties are defined for nodes with
    0 <= Z <= 6.6e-5
with (for example) porosity = 1 for the first two cells,
and porosity 0.35 for cells with Z-borders >= 6.6e-4.
The single cell between Z= 6.e-5 and 6.6e-5 has intermediate
properties between those of the fracture and the rock matrix.
This does not affect the results.

The hydraulic conductivity is calculated from the given
groundwater velocity, in this case v=0.75 m/day, the porosity
which is set eps=1 in the fracture, and a fictive head
difference Delta_h = 0.01 m in L=0.8 m (the length of the model):

  K_x = (v eps) / (Delta_h / L) = 0.0006944 m/s

