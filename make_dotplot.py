import sys
import matplotlib
matplotlib.use('Agg') # Forces headless rendering on the cluster
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import pandas as pd

input_file = sys.argv[1]
output_image = sys.argv[2]

print("Loading data...")
df = pd.read_csv(input_file, sep='\t', header=None, 
                 usecols=[2, 3, 6, 7, 8, 9], 
                 names=['pident', 'length', 'ystart', 'yend', 'xstart', 'xend'])

print("Filtering matches...")
df = df[(df['pident'] >= 95) & (df['length'] >= 25000)]

max_xlen = df[['xstart', 'xend']].max().max()
max_ylen = df[['ystart', 'yend']].max().max()

print("Building plot... (this might take a few seconds)")
segments = [
    [(row.xstart, row.ystart), (row.xend, row.yend)]
    for row in df.itertuples(index=False)
]

fig, ax = plt.subplots(figsize=(6, 6), dpi=300)
lc = LineCollection(segments, colors="black", linewidths=0.5)
ax.add_collection(lc)

plt.margins(0)
if max_ylen > 0:
    plt.ylim(0, max_ylen)
if max_xlen > 0:
    plt.xlim(0, max_xlen)

plt.xlabel("Reference P. infestans genome (bp)")
plt.ylabel("Assembled P. infestans genome (bp)")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

plt.savefig(output_image)
print(f"Success! Image saved as {output_image}")
