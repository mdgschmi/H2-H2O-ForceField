# H2-H2O-ForceField
This code calculates the energies and gradients between a point-like parahydrogen molecule confined in an external clathrate hydrate ("water cage") for use with MMTK.

Found in this directory are:
- H2H2OFF.py : This code is directly called by the main python driver code
- MMTK_h2_h2o.c : This code directly calculates the energies and gradients
- watercom.20 : Data file of the center of masses for the water molecules of the 5^12 ("small") clathrate hydrate
- rotmat.20 : Data file of the rotation matrices for the water molecules of the 5^12 ("small") clathrate hydrate

Missing (Data Files too large)
pot - Potential Energy Surface for H2O-pH2 as a function of radius (r), theta (t), chi (c)
gr, gc, gt - Gradients for H2O-pH2 dV/dr, dV/dt, dV/dc

Please request them at: mdgschmi@uwaterloo.ca
