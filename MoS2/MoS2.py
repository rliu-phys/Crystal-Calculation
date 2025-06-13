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

# Al2O3 = xu.materials.Crystal(
#     "Al2O3",
#     xu.materials.SGLattice(
#         167,             # space group R-3c (No. 167)
#         4.758, 12.991,   # lattice parameters: a and c (hexagonal: b = a, γ = 120° implicit)
#         atoms=["Al", "Al", "O", "O", "O"],
#         pos=[
#             (0.0, 0.0, 0.352),   # Al (Wyckoff 12c, independent position)
#             (0.0, 0.0, 0.852),   # symmetry-related Al (0.352 + 0.5)
#             (0.306, 0.0, 0.25),   # O (Wyckoff 18e, independent position)
#             (0.0, 0.306, 0.75),   # one symmetry equivalent for O
#             (-0.306, -0.306, 0.25)  # another symmetry equivalent for O
#         ]
#     )
# )
MoS2_3R = xu.materials.Crystal.fromCIF("MoS2/MoS2_mp-1434_computed.cif")

MoS2_3R = xu.materials.Crystal(
    "MoS2",
    xu.materials.SGLattice(
        194,               # Space group P6₃/mmc
        3.161, 12.295,     # a, c (in Å)
    ),
    atoms=["Mo", "S"],
    pos=[
        (0.0000, 0.0000, 0.2500),  # Mo at 2c
        (0.3333, 0.6667, 0.6210),  # S at 4f (z ≈ 0.6210)
    ]
)

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
filtered_results = [r for r in results if r["nu (deg)"] > r["theta (deg)"] and r["nu (deg)"] > 0]
# filtered_results = results
df = pd.DataFrame(filtered_results)
print(df)
df.to_csv("MoS2/MoS2_3R_reflections_filtered.txt", sep="\t", index=False)
