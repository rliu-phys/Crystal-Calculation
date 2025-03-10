import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xrayutilities as xu

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

# Define elements and crystal structure for NbOI2
Nb = xu.materials.elements.Nb
O = xu.materials.elements.O
I = xu.materials.elements.I

NbOI2 = xu.materials.Crystal(
    "NbOI2", 
    xu.materials.SGLattice(
        5,  # Space group (monoclinic, unique b-axis)
        15.1882, 3.9329, 7.5232,  # Lattice parameters (a, b, c in Å)
        105.404,  # Monoclinic beta angle (in degrees)
        atoms=["Nb", "I", "I", "O"],
        pos=[
            (0.50079, 0.2905, 0.7895),  # Nb1
            (0.37130, 0.2319, 0.43185),  # I1
            (0.34503, 0.2378, 0.91445),  # I2
            (0.5002, 0.755, 0.7857)  # O1
        ]
    )
)

# Setup goniometer conversion and experimental geometry
qconv = xu.QConversion(('x-', 'z+'), ('z+', 'x-'), (0, 1, 0))
hxrd = xu.Experiment(NbOI2.Q(0,0,-1), NbOI2.Q(1,0,0), qconv=qconv, en=9500)
# If needed, you could alternatively use:
# hxrd = xu.Experiment(NbOI2.Q(0,1,0), NbOI2.Q(28,-1,0), qconv=qconv, en=9500)

# List of reflections (Miller indices)
reflections = [
    (1, 1, 0),
    (2, 0, 2),
    (1, -1, 0),
    (4, 0, -2),
    (2, 0, 0),
    (5, 1, -2),
    (5, -1, -2),
    (6, 0, 0),
    (3, 1, 2),
    (1, -1, -2),
    (3, -1, 2),
    (1, 1, -2),
    (2, 0, -4),
    (3, -1, 0),
    (3, 1, 0),
    (4, 0, 0),
    (7, 1, 0),
    (7, -1, 0),
    (0, 2, 0),
    (0, -2, 0),
    (6, 0, -2),
    (2, 2, 2),
    (4, 2, -2),
    (2, -2, 2),
    (4, -2, -2),
    (8, 0, 0),
    (1, 1, 2),
    (1, -1, 2),
    (6, 0, -1),
    (4, 0, 1),
    (2, 0, -2),
    (3, 1, -4),
    (3, -1, -4),
    (4, 0, 2),
    (0, 0, 3),
    (3, 1, -2),
    (3, -1, -2),
    (2, 2, -4),
    (2, -2, -4),
    (6, 2, 0),
    (5, 1, 1),
    (4, 0, -3),
    (5, -1, 1),
    (6, -2, 0),
    (6, 0, -4),
    (1, -1, 4),
    (9, 1, -4),
    (1, 1, -4),
    (9, -1, -4),
    (1, 1, 4),
    (1, -1, -4),
    (4, 0, 4),
    (4, -2, 0),
    (0, 0, 1),
    (10, 0, -3),
    (5, -1, -4),
    (5, 1, -4),
    (5, 1, 4),
    (5, -1, 4),
    (2, 0, -3),
    (5, 1, -1),
    (8, 2, 0),
    (10, 0, -4),
    (5, -1, -1),
    (12, 0, -2),
    (8, 0, -4),
    (8, -2, 0),
    (7, 1, 3),
    (0, 0, 2),
    (4, 2, 0),
    (7, -1, 3),
    (5, 3, -2),
    (8, 0, -2),
    (11, 1, -2),
    (6, 2, -2),
    (0, 0, 6),
    (2, 0, 4),
    (1, 1, 3),
    (11, -1, -2),
    (2, 0, 5),
    (6, 0, 1),
    (3, -1, -6),
    (3, 1, -6),
    (1, 3, 0),
    (1, -1, 3),
    (5, -3, -2),
    (6, -2, -2),
    (8, 0, -5),
    (1, -3, -2),
    (7, 1, -6),
    (7, -1, -6),
    (8, 0, 3),
    (9, 1, 2),
    (2, -2, -2),
    (3, 3, 2),
    (9, -1, 2),
    (8, 0, 1),
    (6, 0, 2),
    (6, 0, 3),
    (1, -3, 0),
    (10, 0, 2),
    (8, 0, -8),
    (7, 3, 0),
    (1, 3, -2),
    (6, 0, -6),
    (3, -3, 2),
    (11, 1, -3),
    (6, -2, -4),
    (11, -1, -3),
    (15, 1, -5),
    (9, -1, -1),
    (3, -3, 0),
    (12, 0, 1),
    (10, 0, 5),
    (9, 1, -1),
    (15, -1, -5),
    (7, -3, 0),
    (4, 0, -6),
    (14, 0, 3),
    (1, 1, -8),
    (6, 2, -4),
    (1, -1, -8),
    (3, 1, 5),
    (3, 1, -3),
    (6, 2, -1),
    (3, -1, 5),
    (8, 0, 2),
    (3, -1, -3),
    (6, -2, -1),
    (2, 2, -2),
    (3, 3, 0),
    (13, 1, -7),
    (2, 4, -4),
    (17, 1, -3),
    (5, 1, 7),
    (4, 2, 2),
    (4, 2, 4),
    (2, 0, 3),
    (10, 0, -1),
    (17, -1, -3),
    (13, 1, -1),
    (13, -1, -7),
    (12, 2, -2),
    (13, -1, -1),
    (0, 2, 6),
    (5, -1, 7),
    (8, -2, -2),
    (10, 2, -3),
    (0, -2, 6),
    (2, 4, 2),
    (4, -2, 2),
    (0, -2, 2),
    (7, 1, -1),
    (2, -2, 0),
    (4, -2, 4),
    (10, 2, -4),
    (7, -1, -1),
    (16, 0, -5),
    (2, -4, -4),
    (12, -2, -2),
    (9, 1, 5),
    (10, -2, -3),
    (4, 4, -2),
    (1, 1, -3),
    (7, 3, -6),
    (12, 0, -7),
    (6, 0, 4),
    (9, 3, -4),
    (7, -3, -6),
    (10, -2, -4),
    (9, -1, 5),
    (10, 0, -6),
    (1, -1, -3),
    (3, -3, -6),
    (10, 0, -2),
    (6, -2, 2),
    (3, 3, -4),
    (3, 3, -6),
    (5, -1, 2),
    (2, -4, 2),
    (4, 2, 1),
    (1, 1, 6),
    (9, -3, -4),
    (1, -1, 6),
    (0, 2, 3),
    (2, -2, 4),
    (7, 1, -8),
    (12, 0, -1),
    (12, 2, 1),
    (7, -1, -8),
    (11, 1, 1),
    (4, -2, 1),
    (1, -3, 4),
    (3, -3, -4),
    (4, -4, -2),
    (6, 0, 7),
    (12, -2, 1),
    (6, 2, -6),
    (8, 2, -5)
]


# Define goniometer angle bounds: (phi fixed, theta, delta, nu)
bounds = ((0,120), (0), (0,60), (0,60))

results = []  # to store calculated parameters

for hkl in reflections:
    # Compute Q-vector in the crystal (material) frame
    q_mat = NbOI2.Q(hkl)
    # Transform Q-vector to laboratory frame
    q_lab = hxrd.Transform(q_mat)
    # Calculate lattice plane distance (d-spacing)
    d_spacing = NbOI2.planeDistance(hkl)
    # Fit goniometer angles (phi, theta, delta, nu)
    ang, qerror, errcode = xu.Q2AngFit(q_lab, hxrd, bounds)
    # Back-transform to verify the reflection indices
    hkl_back = hxrd.Ang2HKL(*ang, mat=NbOI2)
    
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
df.to_csv("NbOI2_reflections.txt", sep="\t", index=False)