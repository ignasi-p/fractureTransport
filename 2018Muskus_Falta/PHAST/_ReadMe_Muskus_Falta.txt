Simulation using PHAST of the system described in
   Muskus and Falta 2018. Semi-analytical method for matrix diffusion in heterogeneous
   and fractured systems with parent-daughter reactions.
   Journal of Contaminant Hydrology, 218, 94-109.
   https://doi.org/10.1016/j.jconhyd.2018.10.002.

The chemistry database is "%INSTALLDIR%\database\phast.dat"
The "real" input files are
   Muskus_Falta.chem.dat
   Muskus_Falta.dat
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

In the case of Muskus_Falta (2018) we use a non-uniform
Z-grid with the following nodes:

   0    4.8e-005    5.2e-005    5.8e-005    0.0166667
   0.05    0.1    0.16667    0.25

This gives the following cell borders:

   0        2.4e-5  5.0e-4    5.5e-4    0.008359
   0.03333  0.075   0.1333    0.20833   0.25

The fracture half-width is 50 um in Muskus_Falta (2018).
So the first two cells are within the fracture.
In the input file the MEDIA properties are defined for nodes with
    0 <= Z <= 5.3e-4
to set (for example) porosity = 1 for the first two cells,
and porosity 0.1 for cells with Z-borders >= 5.3e-4.
The single cell between Z= 5.e-4 and 5.5e-4 has intermediate
properties between those of the fracture and the rock matrix.
This does not affect the results.

The hydraulic conductivity is calculated from the given
groundwater velocity, in this case v=100 m/yr, the porosity
which is set eps=1 in the fracture, and a fictive head
difference Delta_h = 0.01 m in L=50 m (the length of the model):

  K_x = (v eps) / (Delta_h / L) = 0.015844 m/s
