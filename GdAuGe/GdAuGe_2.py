import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xrayutilities as xu

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
# List of reflections (Miller indices) – top 100 list (example)

# Read the reflection data from the text file (assumes whitespace-delimited)
df_ref = pd.read_csv("GdAuGe/GdAuGe_pristine.txt", delim_whitespace=True)

# Extract the reflection array and the intensity (I) from the file
# (Assuming the file columns are labeled as "h", "k", "l", ..., and "I")
reflections = df_ref[["h", "k", "l"]].values.tolist()
intensity_from_file = df_ref["2θ"].tolist()

# # define reflections by yourself
# reflections = [
#     (1, -1, 4),
# ]
reflections = [(h, k, l) for h in range(-6, 7) 
                           for k in range(-6, 7) 
                           for l in range(-6, 7) 
                           if (h, k, l) != (0, 0, 0)]

GdAuGe = xu.materials.Crystal(
    "GdAuGe",
    xu.materials.SGLattice(
        186,           # space group P6₃mc
        4.4264, 7.4211,    # lattice parameters a and c (hexagonal: b=a, γ=120° implicit)
        atoms=["Gd", "Gd", "Au", "Au", "Ge", "Ge"],
        pos=[
            (0.0, 0.0, 0.0),       # Gd at 2a
            (0.0, 0.0, 0.5),       # symmetry-related Gd
            (1/3, 2/3, 0.633),     # Au at 2b
            (2/3, 1/3, 0.133),     # symmetry-related Au
            (1/3, 2/3, 0.125),     # Ge at 2b
            (2/3, 1/3, 0.625)      # symmetry-related Ge
        ]
    )
)

# Setup goniometer conversion and experimental geometry.
qconv = xu.QConversion(('x-', 'z+'), ('z+', 'x-'), (0, 1, 0))
hxrd = xu.Experiment(GdAuGe.Q(1,0,0), GdAuGe.Q(0,0,1), qconv=qconv, en=9500)

# Define goniometer angle bounds: (theta, phi fixed, delta, nu)
bounds = ((0,120), (27), (0,30), (0,60))

results = []  # to store calculated parameters

for idx, hkl in enumerate(reflections):
    # Compute Q-vector in the crystal frame
    q_mat = GdAuGe.Q(hkl)
    # Transform Q-vector to laboratory frame
    q_lab = hxrd.Transform(q_mat)
    # Calculate lattice plane distance (d-spacing)
    d_spacing = GdAuGe.planeDistance(hkl)
    # Fit goniometer angles
    ang, qerror, errcode = xu.Q2AngFit(q_lab, hxrd, bounds)
    # Back-transform to verify the reflection indices
    hkl_back = hxrd.Ang2HKL(*ang, mat=GdAuGe)
    
    if qerror < 0.01:
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
filtered_results = [r for r in results if r["nu (deg)"] > r["theta (deg)"] and r["nu (deg)"] > 0]

df = pd.DataFrame(filtered_results)
print(df)
df.to_csv("GdAuGe/GdAuGe_reflections_filtered.txt", sep="\t", index=False)
