import pandas as pd
import matplotlib.pyplot as plt
import pysam
import argparse
import re
import os

# Accept command line arguments
parser = argparse.ArgumentParser(description="Calibration sweep for m6A frequency threshold")
parser.add_argument("--m6a", required=True, help="Path to modkit all_modifications.bed (unfiltered)")
parser.add_argument("--ground_truth", required=True, help="Path to ground truth file (ONT BED or GLORI)")
parser.add_argument("--reference", default=None, help="Path to reference fasta (required for DRACH filter with GLORI)")
parser.add_argument("--min_coverage", type=int, default=10, help="Minimum coverage filter (default: 10)")
parser.add_argument("--threshold_min", type=int, default=0, help="Minimum threshold percentage (default: 0)")
parser.add_argument("--threshold_max", type=int, default=100, help="Maximum threshold percentage (default: 100)")
parser.add_argument("--threshold_step", type=int, default=5, help="Threshold step size (default: 5)")
args = parser.parse_args()

# ── DRACH motif check ──
# D=A/G/T, R=A/G, A, C, H=A/C/T (DNA version of DRACH)
DRACH_PATTERN = re.compile(r'[AGT][AG]AC[ACT]')

def is_drach(chrom, pos, strand, fasta):
    """Check if position is in DRACH context. pos is 0-based."""
    try:
        seq = fasta.fetch(chrom, pos - 2, pos + 3).upper()
        if len(seq) != 5:
            return False
        if strand == "-":
            comp = str.maketrans("ACGT", "TGCA")
            seq = seq.translate(comp)[::-1]
        return bool(DRACH_PATTERN.match(seq))
    except:
        return False

# ── Read modkit bedMethyl (unfiltered) ──
modkit_cols = ["chrom", "start", "end", "mod", "score", "strand",
               "start2", "end2", "color", "coverage", "freq",
               "mod_count", "unmod_count", "c1", "c2", "c3", "c4", "c5"]

m6a_all = pd.read_csv(args.m6a, sep="\t", header=None, names=modkit_cols)
m6a_all = m6a_all[(m6a_all["mod"] == "a") & (m6a_all["coverage"] >= args.min_coverage)]
m6a_all["site"] = m6a_all["chrom"] + "_" + m6a_all["start"].astype(str) + "_" + m6a_all["strand"]
print(f"ONT m6A sites (coverage >= {args.min_coverage}): {len(m6a_all)}")

# ── Auto-detect ground truth format ──
gt_path = args.ground_truth
with open(gt_path, "r") as f:
    first_line = f.readline()

if "Chr" in first_line and "Sites" in first_line:
    # GLORI format
    print(f"Detected GLORI format: {gt_path}")
    glori = pd.read_csv(gt_path, sep="\t")
    glori["start_0based"] = glori["Sites"].astype(int) - 1
    if not glori["Chr"].iloc[0].startswith("chr"):
        glori["Chr"] = "chr" + glori["Chr"].astype(str)
    glori["site"] = glori["Chr"] + "_" + glori["start_0based"].astype(str) + "_" + glori["Strand"].astype(str)
    print(f"Total GLORI sites: {len(glori)}")

    # Overlap filter: restrict ground truth to ONT-covered sites only
    ont_sites = set(m6a_all["site"])
    glori_overlapping = glori[glori["site"].isin(ont_sites)]
    print(f"Overlapping sites (ONT covered): {len(glori_overlapping)}")

    # DRACH filter if reference provided
    if args.reference:
        print("Applying DRACH motif filter...")
        fasta = pysam.FastaFile(args.reference)
        drach_mask = glori_overlapping.apply(
            lambda row: is_drach(row["Chr"], row["start_0based"], row["Strand"], fasta), axis=1
        )
        glori_overlapping = glori_overlapping[drach_mask]
        fasta.close()
        print(f"DRACH-filtered overlapping sites: {len(glori_overlapping)}")

    gt_sites = set(glori_overlapping["site"])

else:
    # ONT BED format
    print(f"Detected ONT BED format: {gt_path}")
    gt_df = pd.read_csv(gt_path, sep="\t", header=None,
                         names=["chrom", "start", "end", "mod", "score", "strand"])
    gt_df = gt_df[gt_df["mod"] == "a"]
    gt_df["site"] = gt_df["chrom"] + "_" + gt_df["start"].astype(str) + "_" + gt_df["strand"].fillna("").astype(str)
    gt_sites = set(gt_df["site"])

print(f"Ground truth sites for evaluation: {len(gt_sites)}")

# ── Run calibration sweep ──
thresholds = list(range(args.threshold_min, args.threshold_max + 1, args.threshold_step))
results = []

for threshold in thresholds:
    filtered = m6a_all[m6a_all["freq"] > threshold].copy()
    m6a_sites = set(filtered["site"])

    # For GLORI, only count sites within the overlapping ground truth set
    TP = len(m6a_sites & gt_sites)
    FP = len(m6a_sites - gt_sites)
    FN = len(gt_sites - m6a_sites)

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall    = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1        = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    results.append({
        "threshold": threshold,
        "called": len(m6a_sites),
        "TP": TP, "FP": FP, "FN": FN,
        "Precision": precision,
        "Recall": recall,
        "F1": f1
    })
    print(f"Threshold: {threshold}% | Called: {len(m6a_sites)} | TP: {TP} | FP: {FP} | FN: {FN} | Precision: {precision:.3f} | Recall: {recall:.3f} | F1: {f1:.3f}")

# Convert to dataframe
df = pd.DataFrame(results)

# Find best threshold
best = df.loc[df["F1"].idxmax()]
print(f"\nBest threshold: {best['threshold']:.0f}% | F1: {best['F1']:.3f} | Precision: {best['Precision']:.3f} | Recall: {best['Recall']:.3f}")

# Plot 1: All metrics on one graph
plt.figure(figsize=(10, 6))
plt.plot(df["threshold"], df["Precision"], label="Precision", marker="o", markersize=4)
plt.plot(df["threshold"], df["Recall"], label="Recall", marker="o", markersize=4)
plt.plot(df["threshold"], df["F1"], label="F1 Score", marker="o", markersize=4)
plt.axvline(x=best["threshold"], color="gray", linestyle="--", label=f"Best threshold ({best['threshold']:.0f}%)")
plt.xlabel("Frequency threshold (%)")
plt.ylabel("Score")
plt.title("Calibration sweep: Precision, Recall and F1 Score vs frequency threshold")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("calibration_sweep_combined.png", dpi=150)

# Plot 2: Separate subplots for each metric
fig, axes = plt.subplots(3, 1, figsize=(10, 12))

axes[0].plot(df["threshold"], df["Precision"], color="blue", marker="o", markersize=4)
axes[0].axvline(x=best["threshold"], color="gray", linestyle="--")
axes[0].set_ylabel("Precision")
axes[0].set_title("Precision vs frequency threshold")
axes[0].grid(True)

axes[1].plot(df["threshold"], df["Recall"], color="orange", marker="o", markersize=4)
axes[1].axvline(x=best["threshold"], color="gray", linestyle="--")
axes[1].set_ylabel("Recall")
axes[1].set_title("Recall vs frequency threshold")
axes[1].grid(True)

axes[2].plot(df["threshold"], df["F1"], color="green", marker="o", markersize=4)
axes[2].axvline(x=best["threshold"], color="gray", linestyle="--", label=f"Best threshold ({best['threshold']:.0f}%)")
axes[2].set_xlabel("Frequency threshold (%)")
axes[2].set_ylabel("F1 Score")
axes[2].set_title("F1 Score vs frequency threshold")
axes[2].legend()
axes[2].grid(True)

plt.tight_layout()
plt.savefig("calibration_sweep_separate.png", dpi=150)

# Save results to CSV
df.to_csv("calibration_sweep_results.csv", index=False)
print("Plots saved as calibration_sweep_combined.png and calibration_sweep_separate.png")
print("Results saved as calibration_sweep_results.csv")
