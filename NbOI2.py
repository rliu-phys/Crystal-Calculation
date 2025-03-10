import matplotlib.pyplot as plt
import xrayutilities as xu
import numpy as np
import pandas as pd


Nb = xu.materials.elements.Nb
O = xu.materials.elements.O
I = xu.materials.elements.I

# Define the NbOI2 crystal structure with correct atomic positions
NbOI2 = xu.materials.Crystal(
    "NbOI2", 
    xu.materials.SGLattice(
        5,  # Space group (monoclinic, unique b-axis)
        15.1882, 3.9329, 7.5232,  # Lattice parameters (a, b, c in Å)
        105.404,  # Monoclinic beta angle (only required angle)
        atoms=["Nb", "I", "I", "O"],  # Atomic species
        pos=[
            (0.50079, 0.2905, 0.7895),  # Nb1
            (0.37130, 0.2319, 0.43185),  # I1
            (0.34503, 0.2378, 0.91445),  # I2
            (0.5002, 0.755, 0.7857)  # O1
        ]
    )
)

# #Display unit cell (optional visualization)
# f = plt.figure()
# NbOI2.show_unitcell(fig=f, subplot=111)
# plt.title('NbOI2 Monoclinic Structure')
# plt.show()

# 2S+2D goniometer
qconv=xu.QConversion(('x-','z+'), ('z+', 'x-'), (0, 1, 0))


# [0,1,0] surface normal
#hxrd = xu.Experiment(NbOI2.Q(1,2,0), NbOI2.Q(28,-1,0), qconv=qconv, en= 9500)
hxrd = xu.Experiment(NbOI2.Q(0,1,0), NbOI2.Q(8,0,-1), qconv=qconv, en= 9500)
#hxrd = xu.Experiment(NbOI2.Q(0,1,0), NbOI2.Q(8,-1,0), qconv=qconv, en= 9500)

#hxrd = xu.Experiment(NbOI2.Q(0,1,0), NbOI2.Q(56,0,15), qconv=qconv, en= 9500)
hxrd = xu.Experiment(NbOI2.Q(0,1,0), NbOI2.Q(10,0,-19), qconv=qconv, en= 9500)
hxrd = xu.Experiment(NbOI2.Q(0,-1,0), NbOI2.Q(1,0,0), qconv=qconv, en=9500)
 
hkl=(4,0,0)
q_material = NbOI2.Q(hkl)
q_laboratory = hxrd.Transform(q_material) # transform
print(q_laboratory)

print(f"NbOI2: hkl {hkl}, qvec {np.round(q_material, 5)}") 
print(f"Lattice plane distance: {NbOI2.planeDistance(hkl):.4f}")

#### determine the goniometer angles with the correct geometry restrictions
# tell bounds of angles / (min,max) pair or fixed value for all motors.
# maximum of three free motors! here the first goniometer angle is fixed.

# Phi, Theta, delta, nu

bounds = ((0), (0,100), (-30,60), (0,60))
ang, qerror, errcode = xu.Q2AngFit(q_laboratory, hxrd, bounds)
print(f"err: {errcode:8.3f} ({qerror:.3f}) angles (Phi, Theta, delta, nu)= {np.round(ang, 5)}")
# check that qerror is small!!
print("sanity check with back-transformation (hkl): ",
      np.round(hxrd.Ang2HKL(*ang,mat=NbOI2),5))

print(ang)
#print("501L (hkl): ",
#      np.round(hxrd.Ang2HKL(*ang,mat=Ga2O3),5))


