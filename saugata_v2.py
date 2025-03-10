import matplotlib.pyplot as plt
import xrayutilities as xu
import numpy as np



Ga = xu.materials.elements.Ga
O = xu.materials.elements.O


Ga2O3 = xu.materials.Crystal("GaO", xu.materials.SGLattice(12, 12.2140, 3.0371, 5.7981, 103.83, atoms=["Ga", "Ga", "O", "O", "O"], pos=[(0.0905, 0, 0.7946), (0.15866, 0.5, 0.31402), (0.1645, 0, 0.1098), (0.1733, 0, 0.5632), (-0.0041, 0.5, 0.2566)]))

"""f = plt.figure()
Ga2O3.show_unitcell(fig=f, subplot=111)
plt.title('Ga2O3 monoclinic')"""

# 2S+2D goniometer
qconv=xu.QConversion(('x-','z+'), ('z+', 'x-'), (0, 1, 0))


# [0,1,0] surface normal
#hxrd = xu.Experiment(Ga2O3.Q(1, -2, 0), Ga2O3.Q(32, 1, 0), qconv=qconv, en= 10400)
#hxrd = xu.Experiment(Ga2O3.Q(-1, -2, 0), Ga2O3.Q(32, -1, 0), qconv=qconv, en= 10400)
#hxrd = xu.Experiment(Ga2O3.Q(0, -1, 0), Ga2O3.Q(1, 0, 0), qconv=qconv, en= 10400)
#hxrd = xu.Experiment(Ga2O3.Q(0, 1, 0), Ga2O3.Q(1, 0, 0), qconv=qconv, en= 10400)
#hxrd = xu.Experiment(Ga2O3.Q(-1, 2, 0), Ga2O3.Q(32, 1, 0), qconv=qconv, en= 10400) # right edge
hxrd = xu.Experiment(Ga2O3.Q(1, 2, 0), Ga2O3.Q(32, -1, 0), qconv=qconv, en= 10400) # left edge

"""sub= xu.materials.Crystal("Ga2O3", xu.materials.SGLattice(12, 12.2140, 3.0371, 5.7981, 103.83, atoms=["Ga", "Ga", "O", "O", "O"], pos=[(0.0905, 0, 0.7946), (0.15866, 0.5, 0.31402), (0.1645, 0, 0.1098), (0.1733, 0, 0.5632), (-0.0041, 0.5, 0.2566)]))
hsub = xu.HXRD(sub.Q(0, 0, 1), sub.Q(1, 0, 0))
ax, h = xu.materials.show_reciprocal_space_plane(sub, hsub, ttmax=160)"""

hkl=(2, 0, 0)
q_material = Ga2O3.Q(hkl)
q_laboratory = hxrd.Transform(q_material) # transform
print(q_laboratory)

print(f"Ga2O3: hkl {hkl}, qvec {np.round(q_material, 5)}") 
print(f"Lattice plane distance: {Ga2O3.planeDistance(hkl):.4f}")

#### determine the goniometer angles with the correct geometry restrictions
# tell bounds of angles / (min,max) pair or fixed value for all motors.
# maximum of three free motors! here the first goniometer angle is fixed.

# Phi, Theta, delta, nu

bounds = ((0), (-30,30), (-30,60), (0,40))
ang, qerror, errcode = xu.Q2AngFit(q_laboratory, hxrd, bounds)
print(f"err: {errcode:8.3f} ({qerror:.3f}) angles (Phi, Theta, delta, nu)= {np.round(ang, 5)}")
# check that qerror is small!!
print("sanity check with back-transformation (hkl): ",
      np.round(hxrd.Ang2HKL(*ang,mat=Ga2O3),5))

print(ang)
#print("501L (hkl): ",
#      np.round(hxrd.Ang2HKL(*ang,mat=Ga2O3),5))


