import matplotlib.pyplot as plt
import xrayutilities as xu
import numpy as np



AlN = xu.materials.Crystal.fromCIF("AlScN.cif")


# 2S+2D goniometer
qconv=xu.QConversion(('y+','x-','z+'), ('y+', 'x-'), (0, 0, 1))



hxrd = xu.Experiment(AlN.Q(1, 0, 0), AlN.Q(0, 0, 1), qconv=qconv, en= 9500, sampleor="x+")

hkl=(1, 0, 3)
q_material = AlN.Q(hkl)
q_laboratory = hxrd.Transform(q_material) # transform
print(q_laboratory)

arrr = np.zeros((360,5))

for iii in range(360):

    bounds = ((0, 60), (iii-180), (0), (0,60), (0,30))
    ang, qerror, errcode = xu.Q2AngFit(q_laboratory, hxrd, bounds)
    if qerror < 0.1:
        arrr[iii] = ang
    #print(f"err: {errcode:8.3f} ({qerror:.3f}) angles (Theta, Phi, nu, del)= {np.round(ang, 5)}")
    # check that qerror is small!!
    #print("sanity check with back-transformation (hkl): ",
      #np.round(hxrd.Ang2HKL(*ang,mat=AlN),5))



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
