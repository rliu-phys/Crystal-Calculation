import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xrayutilities as xu

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

GdAuGe = xu.materials.Crystal(
    "GdAuGe",
    xu.materials.SGLattice(
        186,           # space group number for P6₃mc
        4.4264, 7.4211,    # lattice parameters a and c (hexagonal, so b = a and γ = 120° are implicit)
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


# Setup goniometer conversion and experimental geometry
qconv = xu.QConversion(('x-', 'z+'), ('z+', 'x-'), (0, 1, 0))
hxrd = xu.Experiment(GdAuGe.Q(1,0,0), GdAuGe.Q(0,0,1), qconv=qconv, en=9500)
# If needed, you could alternatively use:
# hxrd = xu.Experiment(GdAuGe.Q(0,1,0), GdAuGe.Q(28,-1,0), qconv=qconv, en=9500)

# List of reflections (Miller indices)
reflections = [(0,0,4)]
reflections = [
    (0, 0, 4),
    (0, 0, 2),
    (0, 0, 6),
    (1, 0, 4),
    (0, 1, 4),
    (1, 0, 2),
    (0, 1, 2),
    (1, 1, 0),
    (2, 0, 0),
    (1, -1, 0),
    (1, 1, 2),
    (0, 0, 8),
    (1, 0, 6),
    (0, 1, 6),
    (2, 0, 2),
    (1, 1, 4),
    (2, 1, 0),
    (1, -1, 2),
    (0, 2, 0),
    (0, 0, 10),
    (2, 0, 4),
    (1, 0, 8),
    (0, 1, 8),
    (2, 2, 0),
    (3, 0, 0),
    (1, 1, 6),
    (1, 0, 10),
    (0, 1, 10),
    (1, -1, 4),
    (2, 0, 6),
    (2, 1, 2),
    (1, 2, 0),
    (1, 1, -2),
    (0, 2, 2),
    (0, 0, 12),
    (1, 0, 12),
    (0, 1, 12),
    (2, 0, 8),
    (1, 2, 2),
    (2, 1, 4),
    (3, 0, 2),
    (2, 0, -2),
    (1, -1, 6),
    (0, 2, 4),
    (1, 0, -2),
    (0, 0, -4),
    (0, 0, -2),
    (0, 0, -6),
    (1, 0, -4),
    (0, 1, -4),
    (1, 0, -2),
    (0, 1, -2),
    (1, 1, 0),    # duplicate reflection may occur in the intensity ranking
    (2, 0, 0),
    (1, -1, 0),
    (1, 1, -2),
    (0, 0, -8),
    (1, 0, -6),
    (0, 1, -6),
    (2, 0, -2),
    (1, 1, -4),
    (2, 1, 0),
    (1, -1, -2),
    (0, 2, 0),
    (0, 0, -10),
    (2, 0, -4),
    (1, 0, -8),
    (0, 1, -8),
    (2, 2, 0),
    (3, 0, 0),
    (1, 1, -6),
    (1, 0, -10),
    (0, 1, -10),
    (1, -1, -4),
    (2, 0, -6),
    (2, 1, -2),
    (1, 2, 0),
    (1, 1, 2),
    (0, 2, -2),
    (0, 0, 14),
    (1, 0, 14),
    (0, 1, 14),
    (2, 0, 10),
    (1, 2, 4),
    (2, 1, 6),
    (3, 0, 4),
    (2, 0, 0),
    (1, -1, 8),
    (0, 2, 6),
    (1, 0, 16),
    (0, 1, 16),
    (2, 0, 12),
    (1, 1, 8),
    (2, 2, 4),
    (3, 0, 6),
    (1, 2, 6),
    (0, 0, -12),
    (1, 0, -12),
    (0, 1, -12),
    (1, 1, -8)
]



# Define goniometer angle bounds: (theta, phi fixed, delta, nu)
bounds = ((20,25), (0,60), (0,30), (0,60))

results = []  # to store calculated parameters

for hkl in reflections:
    # Compute Q-vector in the crystal (material) frame
    q_mat = GdAuGe.Q(hkl)
    # Transform Q-vector to laboratory frame
    q_lab = hxrd.Transform(q_mat)
    # Calculate lattice plane distance (d-spacing)
    d_spacing = GdAuGe.planeDistance(hkl)
    # Fit goniometer angles (phi, theta, delta, nu)
    ang, qerror, errcode = xu.Q2AngFit(q_lab, hxrd, bounds)
    # Back-transform to verify the reflection indices
    hkl_back = hxrd.Ang2HKL(*ang, mat=GdAuGe)
    
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
        "hkl_back": str(np.round(hkl_back, 4))
    })

# Convert results to a DataFrame and output as a text file
df = pd.DataFrame(results)
print(df)

# Save table into a tab-delimited text file
df.to_csv("GdAuGe_reflections.txt", sep="\t", index=False)

# # Calculate the hkl
# hkl_back = hxrd.Ang2HKL(22.7, -10, 20, 31, mat=GdAuGe)
# print(hkl_back)