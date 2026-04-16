import pandas as pd
import matplotlib.pyplot as plt
import argparse

# Accept command line arguments
parser = argparse.ArgumentParser(description="Calibration sweep for m6A frequency threshold")
parser.add_argument("--m6a", required=True, help="Path to modkit output .bed file")
parser.add_argument("--ground_truth", required=True, help="Path to ground truth .bed file")
args = parser.parse_args()

# Read files
modkit_cols = ["chrom", "start", "end", "mod", "score", "strand",
               "start2", "end2", "color", "coverage", "freq",
               "mod_count", "unmod_count", "c1", "c2", "c3", "c4", "c5"]

m6a_all = pd.read_csv(args.m6a, sep="\t", header=None, names=modkit_cols)
ground_truth = pd.read_csv(args.ground_truth, sep="\t", header=None,
                            names=["chrom", "start", "end", "mod", "score", "strand"])

# Filter m6A only
m6a_all = m6a_all[m6a_all["mod"] == "a"]
ground_truth = ground_truth[ground_truth["mod"] == "a"]

# Create ground truth site identifiers
ground_truth["site"] = ground_truth["chrom"].astype(str) + "_" + ground_truth["start"].astype(str) + "_" + ground_truth["strand"].fillna("").astype(str)
gt_sites = set(ground_truth["site"])

# Run calibration sweep from 70% to 100% at 1% intervals
thresholds = list(range(70, 101))
results = []

for threshold in thresholds:
    filtered = m6a_all[m6a_all["freq"] > threshold].copy()
    filtered["site"] = filtered["chrom"].astype(str) + "_" + filtered["start"].astype(str) + "_" + filtered["strand"].fillna("").astype(str)
    m6a_sites = set(filtered["site"])

    TP = len(m6a_sites & gt_sites)
    FP = len(m6a_sites - gt_sites)
    FN = len(gt_sites - m6a_sites)

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall    = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1        = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    results.append({
        "threshold": threshold,
        "TP": TP, "FP": FP, "FN": FN,
        "Precision": precision,
        "Recall": recall,
        "F1": f1
    })
    print(f"Threshold: {threshold}% | Precision: {precision:.3f} | Recall: {recall:.3f} | F1: {f1:.3f}")

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
plt.show()

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
plt.show()

print("Plots saved as calibration_sweep_combined.png and calibration_sweep_separate.png")