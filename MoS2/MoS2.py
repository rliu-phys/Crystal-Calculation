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

MoS2_3R = xu.materials.Crystal.fromCIF("MoS2/sd_0526530.cif")

# Setup goniometer conversion and experimental geometry.
qconv = xu.QConversion(('x-', 'z+'), ('z+', 'x-'), (0, 1, 0))
hxrd = xu.Experiment(MoS2_3R.Q(1,1,0), MoS2_3R.Q(0,0,1), qconv=qconv, en=9659)

# Define goniometer angle bounds: (theta, phi fixed, delta, nu)
bounds = ((0,130), (0), (0,40), (0,90))

results = []  # to store calculated parameters

for idx, hkl in enumerate(reflections):
    # Compute Q-vector in the crystal frame
    q_mat = MoS2_3R.Q(hkl)
    # Transform Q-vector to laboratory frame
    q_lab = hxrd.Transform(q_mat)
    # Calculate lattice plane distance (d-spacing)
    d_spacing = MoS2_3R.planeDistance(hkl)
    # Fit goniometer angles
    ang, qerror, errcode = xu.Q2AngFit(q_lab, hxrd, bounds)
    # Back-transform to verify the reflection indices
    hkl_back = hxrd.Ang2HKL(*ang, mat=MoS2_3R)
    
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
df.to_csv("MoS2/MoS2_3R_reflections_filtered.txt", sep="\t", index=False)
