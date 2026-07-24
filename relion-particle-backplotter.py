from pathlib import Path

import matplotlib.pyplot as plt
import mrcfile
import numpy as np
import starfile

## change the below information
STAR_FILE = r"pathway to your run_it025_data.star file"
MICROGRAPH = r"pathway to your .mrc micrograph"
CLASSES = [1, 2, 5, 6] #choose which class number you want to display
TUBE_ID = 4 #choose which specific fibril you want to display
OUTPUT = r"pathway to your output folder"



MARKER_SIZE = 30
ALPHA = 0.8
DOWNSAMPLE = 1.0
particles = starfile.read(STAR_FILE)

if isinstance(particles, dict):
    if "particles" in particles:
        particles = particles["particles"]
    else:
        particles = list(particles.values())[0]
        
if "rlnMicrographName" in particles.columns:

    micrograph_name = Path(MICROGRAPH).name

    particles = particles[
        particles["rlnMicrographName"].apply(
            lambda x: Path(str(x)).name == micrograph_name
        )
    ]
    
if TUBE_ID is not None:

    particles = particles[
        particles["rlnHelicalTubeID"] == TUBE_ID
    ]

particles = particles[
    particles["rlnClassNumber"].isin(CLASSES)
]

print()

print("Particles remaining:", len(particles))

print()

for c in CLASSES:

    print(
        f"Class {c}:",
        np.sum(particles["rlnClassNumber"] == c)
    )

print()

with mrcfile.open(MICROGRAPH) as mrc:

    image = np.squeeze(mrc.data)

vmin, vmax = np.percentile(image, [1, 99])

fig, ax = plt.subplots(figsize=(12, 12))

ax.imshow(
    image,
    cmap="gray",
    origin="upper",
    vmin=vmin,
    vmax=vmax,
)

colours = [
    "red",
    "dodgerblue",
    "limegreen",
    "gold",
    "magenta",
    "cyan",
    "orange",
    "white",
    "purple",
    "brown",
]

for i, class_number in enumerate(CLASSES):

    subset = particles[
        particles["rlnClassNumber"] == class_number
    ]

    ax.scatter(
        subset["rlnCoordinateX"] * DOWNSAMPLE,
        subset["rlnCoordinateY"] * DOWNSAMPLE,
        s=MARKER_SIZE,
        color=colours[i % len(colours)],
        alpha=ALPHA,
        edgecolors="black",
        linewidths=0.3,
        label=f"Class {class_number}",
    )

ax.set_title("RELION particle classes")

ax.legend()

plt.tight_layout()

plt.savefig(
    OUTPUT,
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print()

print("Finished!")

print("Image saved to:")

print(OUTPUT)
