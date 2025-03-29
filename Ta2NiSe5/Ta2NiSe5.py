import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xrayutilities as xu

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
# List of reflections (Miller indices) – top 100 list (example)

# # Read the reflection data from the text file (assumes whitespace-delimited)
# df_ref = pd.read_csv("Al2O3/Al2O3.txt", delim_whitespace=True)

# # Extract the reflection array and the intensity (I) from the file
# # (Assuming the file columns are labeled as "h", "k", "l", ..., and "I")
# reflections = df_ref[["h", "k", "l"]].values.tolist()
# intensity_from_file = df_ref["2θ"].tolist()

# # define reflections by yourself
# reflections = [
#     (1, -1, 4),
# ]
reflections = [(h, k, l) for h in range(-9, 10) 
                           for k in range(-9, 10) 
                           for l in range(-9, 10) 
                           if (h, k, l) != (0, 0, 0)]

# reflections = [
#     (1, 1, 6),
# ]

Ta2NiSe5 = xu.materials.Crystal(
    "Ta2NiSe5",
    xu.materials.SGLattice(
        15,                   # space group number for C2/c
        3.496, 12.829, 15.641, 90.53,   # lattice parameters: a, b, c, beta (monoclinic; α=γ=90° are implicit)
        atoms=["Ta", "Ta", "Ni", "Se", "Se", "Se"],
        pos=[
            (0.0,    0.167, 0.25),    # Ta1 (Wyckoff 4e)
            (0.0,    0.668, 0.25),    # Ta2 (Wyckoff 4e)
            (0.0,    0.5,   0.0),     # Ni  (Wyckoff 4e)
            (0.277,  0.103, 0.119),   # Se1 (e.g. Wyckoff 8f)
            (0.277,  0.103, 0.619),   # Se2 (symmetry partner of Se1)
            (0.0,    0.5,   0.75)     # Se3 (e.g. Wyckoff 4e)
        ]
    )
)


# Setup goniometer conversion and experimental geometry.
qconv = xu.QConversion(('x-', 'z+'), ('z+', 'x-'), (0, 1, 0))
hxrd = xu.Experiment(Ta2NiSe5.Q(0,0,1), Ta2NiSe5.Q(0,1,0), qconv=qconv, en=9500)

# Define goniometer angle bounds: (theta, phi fixed, delta, nu)
bounds = ((0,120), (0), (0,30), (0,60))

results = []  # to store calculated parameters

for idx, hkl in enumerate(reflections):
    # Compute Q-vector in the crystal frame
    q_mat = Ta2NiSe5.Q(hkl)
    # Transform Q-vector to laboratory frame
    q_lab = hxrd.Transform(q_mat)
    # Calculate lattice plane distance (d-spacing)
    d_spacing = Ta2NiSe5.planeDistance(hkl)
    # Fit goniometer angles
    ang, qerror, errcode = xu.Q2AngFit(q_lab, hxrd, bounds)
    # Back-transform to verify the reflection indices
    hkl_back = hxrd.Ang2HKL(*ang, mat=Ta2NiSe5)
    
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
df.to_csv("Ta2NiSe5/Ta2NiSe5_reflections_filtered.txt", sep="\t", index=False)
