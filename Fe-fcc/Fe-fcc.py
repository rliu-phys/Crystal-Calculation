import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xrayutilities as xu

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

reflections = [(h, k, l) for h in range(-9, 10) 
                           for k in range(-9, 10) 
                           for l in range(-9, 10) 
                           if (h, k, l) != (0, 0, 0)]
Fe = xu.materials.Crystal.fromCIF("Fe-fcc/fE_fcc.cif")

# Setup goniometer conversion and experimental geometry.
qconv = xu.QConversion(('x-', 'z+'), ('z+', 'x-'), (0, 1, 0))
hxrd = xu.Experiment(Fe.Q(1,0,0), Fe.Q(0,0,1), qconv=qconv, en=10500)

# Define goniometer angle bounds: (theta, phi fixed, delta, nu)
bounds = ((0,150), (0), (0,40), (0,70))

results = []  # to store calculated parameters

for idx, hkl in enumerate(reflections):
    # Compute Q-vector in the crystal frame
    q_mat = Fe.Q(hkl)
    # Transform Q-vector to laboratory frame
    q_lab = hxrd.Transform(q_mat)
    # Calculate lattice plane distance (d-spacing)
    d_spacing = Fe.planeDistance(hkl)
    # Fit goniometer angles
    ang, qerror, errcode = xu.Q2AngFit(q_lab, hxrd, bounds)
    # Back-transform to verify the reflection indices
    hkl_back = hxrd.Ang2HKL(*ang, mat=Fe)
    
    if qerror < 0.1:
        results.append({
            "h": hkl[0],
            "k": hkl[1],
            "l": hkl[2],
            "d_spacing (Å)": round(d_spacing, 4),
            "theta (deg)": round(ang[0], 4),
            "phi (deg)": round(ang[1], 4),
            "delta (deg)": round(ang[2], 4),
            "nu (deg)": round(ang[3], 4),
            "qerror": round(qerror, 4),
            "errcode": errcode,
            "hkl_back": str(np.round(hkl_back, 4)),
            # "I": intensity_from_file[idx]
        })
# filtered_results = [r for r in results if r["nu (deg)"] > r["theta (deg)"] and r["nu (deg)"] > 0]
filtered_results = results
df = pd.DataFrame(filtered_results)
print(df)
df.to_csv("Fe-fcc/Fe_reflections_filtered.txt", sep="\t", index=False)

hkl=(2, 0, 0.05)
q_material = Fe.Q(hkl)
q_laboratory = hxrd.Transform(q_material) # transform
print(q_laboratory)

print(f"Fe: hkl {hkl}, qvec {np.round(q_material, 5)}") 
print(f"Lattice plane distance: {Fe.planeDistance(hkl):.4f}")

ang, qerror, errcode = xu.Q2AngFit(q_laboratory, hxrd, bounds)
print(f"err: {errcode:8.3f} ({qerror:.3f}) angles (Theta, phi, delta, nu)= {np.round(ang, 5)}")
print("sanity check with back-transformation (hkl): ",
      np.round(hxrd.Ang2HKL(*ang,mat=Fe),5))