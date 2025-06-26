import matplotlib.pyplot as plt
import xrayutilities as xu
import numpy as np



W = xu.materials.elements.W
Se = xu.materials.elements.Se


MoS2 = xu.materials.Crystal.fromCIF("MoS2/sd_0526530.cif")


# 2S+2D goniometer
qconv=xu.QConversion(('x-','z+'), ('z+', 'x-'), (0, 1, 0))

fout = open("MoS2/MoS2.txt", "w")

hxrd = xu.Experiment(MoS2.Q(1, 1, 0), MoS2.Q(0, 0, 1), qconv=qconv, en= 9659) 

for h in range(-5, 6, 1):
    for k in range(-5, 6, 1):
        for l in range(0, 10):
            if l%2 and not ((h+2*k)%3):
                continue
            hkl=(h, k, l)
            #print(hkl)
            q_material = MoS2.Q(hkl)
            q_laboratory = hxrd.Transform(q_material) # transform

            bounds = ((0,130), (0), (0,40), (0,60))
            ang, qerror, errcode = xu.Q2AngFit(q_laboratory, hxrd, bounds);
            if qerror<1e-3:
                h1, k1, l1 = np.round(hxrd.Ang2HKL(*ang,mat=MoS2),5)
                fout.write("{0:6.0f} {1:6.0f} {2:6.0f} {3:6.1f} {4:6.1f} {5:6.1f} {6:7.3f} {7:7.3f} {8:7.3f} {9:7.3f}\n".format(h,k,l,h1,k1,l1,ang[0],ang[1],ang[2],ang[3]))
            #print(f"err: {errcode:8.3f} ({qerror:.3f}) angles (Phi, Theta, delta, nu)= {np.round(ang, 5)}")
            ## check that qerror is small!!
            #print("sanity check with back-transformation (hkl): ",



"""
hkl=(0, 0, 8)
q_material = MoS2.Q(hkl)
q_laboratory = hxrd.Transform(q_material) # transform
print(q_laboratory)

print(f"MoS2: hkl {hkl}, qvec {np.round(q_material, 5)}") 
print(f"Lattice plane distance: {MoS2.planeDistance(hkl):.4f}")

#### determine the goniometer angles with the correct geometry restrictions
# tell bounds of angles / (min,max) pair or fixed value for all motors.
# maximum of three free motors! here the first goniometer angle is fixed.

# Theta, phi, delta, nu

bounds = ((0,50), (0), (0,0), (0,60))
ang, qerror, errcode = xu.Q2AngFit(q_laboratory, hxrd, bounds)
print(f"err: {errcode:8.3f} ({qerror:.3f}) angles (Phi, Theta, delta, nu)= {np.round(ang, 5)}")
# check that qerror is small!!
print("sanity check with back-transformation (hkl): ",
      np.round(hxrd.Ang2HKL(*ang,mat=MoS2),5))

print(ang)
"""