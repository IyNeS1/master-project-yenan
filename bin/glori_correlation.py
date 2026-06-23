import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--modkit", required=True)
parser.add_argument("--glori", required=True)
parser.add_argument("--min_coverage", type=int, default=10)
args = parser.parse_args()

# Read modkit
modkit_cols = ["chrom", "start", "end", "mod", "score", "strand",
               "start2", "end2", "color", "coverage", "freq",
               "mod_count", "unmod_count", "c1", "c2", "c3", "c4", "c5"]
modkit = pd.read_csv(args.modkit, sep="\t", header=None, names=modkit_cols)
modkit = modkit[(modkit["mod"] == "a") & (modkit["coverage"] >= args.min_coverage)]
modkit["site"] = modkit["chrom"] + "_" + modkit["start"].astype(str) + "_" + modkit["strand"]

# Read GLORI
glori = pd.read_csv(args.glori, sep="\t")
glori["start_0based"] = glori["Sites"] - 1
glori["site"] = glori["Chr"] + "_" + glori["start_0based"].astype(str) + "_" + glori["Strand"]

# Merge on shared sites
merged = pd.merge(
    modkit[["site", "freq"]].rename(columns={"freq": "pipeline_freq"}),
    glori[["site", "Ratio"]].rename(columns={"Ratio": "glori_ratio"}),
    on="site"
)
# Convert GLORI ratio (0-1) to percentage (0-100)
merged["glori_pct"] = merged["glori_ratio"] * 100

print(f"Overlapping sites: {len(merged)}")

# Spearman correlation
rho, pval = stats.spearmanr(merged["pipeline_freq"], merged["glori_pct"])
print(f"Spearman correlation: {rho:.4f} (p = {pval:.2e})")

# Pearson correlation
r, pval_p = stats.pearsonr(merged["pipeline_freq"], merged["glori_pct"])
print(f"Pearson correlation: {r:.4f} (p = {pval_p:.2e})")

# Scatter plot
plt.figure(figsize=(8, 8))
plt.scatter(merged["glori_pct"], merged["pipeline_freq"], alpha=0.05, s=1)
plt.plot([0, 100], [0, 100], "r--", label="y = x")
plt.xlabel("GLORI modification frequency (%)")
plt.ylabel("Pipeline modification frequency (%)")
plt.title(f"Pipeline vs GLORI (n={len(merged)}, Spearman rho={rho:.3f})")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("glori_correlation.png", dpi=150)
print("Plot saved: glori_correlation.png")
