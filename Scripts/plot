"""Density visualization of side-chain angles in tripeptide (XRX) helices,
classified based on size group of the preceding residue.

Equivalent to plot_angles.R but implemented using matplotlib.
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


input_file = sys.argv[1] if len(sys.argv) > 1 else "final/angles.tsv"
output_img = sys.argv[2] if len(sys.argv) > 2 else "final/angle_density_plot.png"

SIZE_GROUPS = ["Tiny", "Small", "Intermediate", "Large", "Bulky"]

COLOR_MAP = {
    "Tiny": "#e6e2d3",
    "Small": "#f2c879",
    "Intermediate": "#f28e5c",
    "Large": "#e4572e",
    "Bulky": "#b11226",
}

data = pd.read_csv(input_file, sep="\t")

# normalize angle to [-180, 180]
data["angle"] = ((data["angle"] + 180) % 360) - 180

# keep only valid size groups
data = data[data["size_class"].isin(SIZE_GROUPS)]

total_count = len(data)

fig, axis = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor("#bdbdbd")
axis.set_facecolor("#bdbdbd")

x_vals = np.linspace(-180, 180, 1000)

for group in SIZE_GROUPS:
    subset_angles = data[data["size_class"] == group]["angle"].values

    if len(subset_angles) < 2:
        continue

    kde_estimator = gaussian_kde(subset_angles, bw_method="silverman")
    density = kde_estimator(x_vals)

    axis.plot(
        x_vals,
        density,
        color=COLOR_MAP[group],
        linewidth=1.6,
        label=group
    )

axis.set_xlim(-180, 180)
axis.set_xticks(np.arange(-180, 181, 50))

# vertical reference lines
for val in np.arange(-180, 181, 50):
    axis.axvline(val, linestyle=":", linewidth=0.5, alpha=0.7, color="blue")

axis.set_xlabel(r"Angle between adjacent C-$\alpha$ $\rightarrow$ Centroid vectors [°]")
axis.set_ylabel("Normalized Frequency [A.U.]")
axis.set_title(f"Tripeptide (XRX) Helix Distribution (n = {total_count})")

legend = axis.legend(loc="upper left", facecolor="white", edgecolor="black")
legend.get_frame().set_linewidth(0.8)

axis.spines["top"].set_visible(False)
axis.spines["right"].set_visible(False)

os.makedirs(os.path.dirname(output_img) or ".", exist_ok=True)

plt.tight_layout()
plt.savefig(output_img, dpi=300, facecolor=fig.get_facecolor())

print(f"Saved output: {output_img}")
print(f"Total samples: {total_count}")
