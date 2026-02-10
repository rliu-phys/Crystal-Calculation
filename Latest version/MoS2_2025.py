import matplotlib.pyplot as plt
import xrayutilities as xu
import numpy as np



Mo = xu.materials.elements.Mo
S = xu.materials.elements.S

MoS2 = xu.materials.Crystal("MoS2-3R", xu.materials.SGLattice(166, 3.1585273, 18.37, atoms=["Mo","Mo", "Mo", "S","S","S","S","S","S"], pos=[(0.00000, 0.00000, 0.00000),
                                                                                                                                       (0.00480, 0.00480, 0.00750),
                                                                                                                                       (0.00240, 0.00000, 0.00000),
                                                                                                                                       (0.00000, 0.00000, 0.25160),
                                                                                                                                       (0.00380, 0.00380, 0.00510),
                                                                                                                                       (0.00190, 0.00000, 0.00000),
                                                                                                                                       (0.00000, 0.00000, 0.41510),
                                                                                                                                       (0.00380, 0.00380, 0.00510),
                                                                                                                                       (0.00190, 0.00000, 0.00000)])) #space group, lattice parameter, position of atoms


# 2S+2D goniometer
qconv=xu.QConversion(('y+','x-','z+'), ('y+', 'x-'), (0, 0, 1))

#fout = open("WSe2.txt", "w")

hxrd = xu.Experiment(MoS2.Q(1, 1, 0), MoS2.Q(0, 0, 1), qconv=qconv, en= 9500, sampleor="x+")

"""

hkl=(0, 0, 12)
q_material = MoS2.Q(hkl)
q_laboratory = hxrd.Transform(q_material) # transform
print(q_laboratory)


bounds = ((0, 130), (0), (-5,5), (0,60), (1.11))
ang, qerror, errcode = xu.Q2AngFit(q_laboratory, hxrd, bounds)
print(f"err: {errcode:8.3f} ({qerror:.3f}) angles (Theta, Phi, nu, del)= {np.round(ang, 5)}")
# check that qerror is small!!
print("sanity check with back-transformation (hkl): ",
      np.round(hxrd.Ang2HKL(*ang,mat=MoS2),5))

print(ang)



d2 = 1000

for d_th in np.arange(0, 1, 0.01):
    for d_phi in np.arange(-5, -4, 0.01):
        ang = [114+d_th,  d_phi, 1.30202, 48.7, -3.81]
        hkl = np.round(hxrd.Ang2HKL(*ang,mat=MoS2),5)
        d3 = (hkl[0]-1)**2+(hkl[1]-1)**2+(hkl[2]-0)**2
        if d3 < d2:
            d1 = hkl
            d2 = d3
            d4 = [d_th, d_phi]
print(d2,d1, d4)

"""



hkl=(1, 0, 1)
q_material = MoS2.Q(hkl)
q_laboratory = hxrd.Transform(q_material) # transform
print(q_laboratory)


bounds = ((0, 130), (-4.61), (1.30202), (0,60), (11.21))
ang, qerror, errcode = xu.Q2AngFit(q_laboratory, hxrd, bounds)
print(f"err: {errcode:8.3f} ({qerror:.3f}) angles (Theta, Phi, nu, del)= {np.round(ang, 5)}")
# check that qerror is small!!
print("sanity check with back-transformation (hkl): ",
      np.round(hxrd.Ang2HKL(*ang,mat=MoS2),5))

print(ang)



"""

# sannity checks

hkl=(0, 0, 12)
q_material = MoS2.Q(hkl)
q_laboratory = hxrd.Transform(q_material) # transform
#print(q_laboratory)


bounds = ((0, 130), (-4.61), (1.30202), (0,60), (1.11))
ang, qerror, errcode = xu.Q2AngFit(q_laboratory, hxrd, bounds)
#print(f"err: {errcode:8.3f} ({qerror:.3f}) angles (Theta, Phi, nu, del)= {np.round(ang, 5)}")
# check that qerror is small!!
print("sanity check with back-transformation (hkl): ",
      np.round(hxrd.Ang2HKL(*ang,mat=MoS2),5))


hkl=(1, 1, 0)
q_material = MoS2.Q(hkl)
q_laboratory = hxrd.Transform(q_material) # transform
#print(q_laboratory)


bounds = ((0, 130), (-4.61), (1.30202), (48.7), (-3.81))
ang, qerror, errcode = xu.Q2AngFit(q_laboratory, hxrd, bounds)
#print(f"err: {errcode:8.3f} ({qerror:.3f}) angles (Theta, Phi, nu, del)= {np.round(ang, 5)}")
# check that qerror is small!!
print("sanity check with back-transformation (hkl): ",
      np.round(hxrd.Ang2HKL(*ang,mat=MoS2),5))
print(ang)
"""
