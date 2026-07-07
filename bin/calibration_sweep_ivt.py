import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

# Accept command line arguments
parser = argparse.ArgumentParser(description="Calibration sweep for m6A frequency threshold, with IVT negative control support")
parser.add_argument("--m6a", required=True, help="Path to modkit output .bed file")
parser.add_argument("--ground_truth", required=True, help="Path to ground truth file (ONT BED, GLORI CSV, or empty file for IVT negative control)")
parser.add_argument("--threshold_min", type=int, default=5, help="Minimum threshold percentage (default: 5)")
parser.add_argument("--threshold_max", type=int, default=100, help="Maximum threshold percentage (default: 100)")
parser.add_argument("--min_coverage", type=int, default=10, help="Minimum coverage per site (default: 10)")
args = parser.parse_args()

# Read modkit bedMethyl
modkit_cols = ["chrom", "start", "end", "mod", "score", "strand",
               "start2", "end2", "color", "coverage", "freq",
               "mod_count", "unmod_count", "c1", "c2", "c3", "c4", "c5"]

m6a_all = pd.read_csv(args.m6a, sep="\t", header=None, names=modkit_cols)
m6a_all = m6a_all[m6a_all["mod"] == "a"]

# Coverage filter
m6a_all = m6a_all[m6a_all["coverage"] >= args.min_coverage]
print(f"Sites after coverage filter (min_coverage={args.min_coverage}): {len(m6a_all)}")

# Auto-detect ground truth format, and handle empty file for IVT negative control
gt_path = args.ground_truth
gt_columns = ["chrom", "start", "strand"]

if os.path.getsize(gt_path) == 0:
    print(f"Ground truth file is empty (IVT negative control mode): {gt_path}")
    gt_df = pd.DataFrame(columns=gt_columns)
else:
    with open(gt_path, "r") as f:
        first_line = f.readline()

    if "Chr" in first_line and "Sites" in first_line:
        # GLORI format: tab-separated with header
        print(f"Detected GLORI format: {gt_path}")
        gt_raw = pd.read_csv(gt_path, sep="\t")
        # GLORI Sites is 1-based; bedMethyl start is 0-based
        gt_df = pd.DataFrame({
            "chrom": gt_raw["Chr"].astype(str),
            "start": gt_raw["Sites"].astype(int) - 1,
            "strand": gt_raw["Strand"].astype(str)
        })
        if not gt_df["chrom"].iloc[0].startswith("chr"):
            gt_df["chrom"] = "chr" + gt_df["chrom"]
    else:
        # ONT BED format: tab-separated, no header, 6 columns
        print(f"Detected ONT BED format: {gt_path}")
        gt_df = pd.read_csv(gt_path, sep="\t", header=None,
                             names=["chrom", "start", "end", "mod", "score", "strand"])
        gt_df = gt_df[gt_df["mod"] == "a"]
        gt_df = gt_df[gt_columns]

# Create ground truth site identifiers
gt_df["site"] = gt_df["chrom"].astype(str) + "_" + gt_df["start"].astype(str) + "_" + gt_df["strand"].fillna("").astype(str)
gt_sites = set(gt_df["site"])
print(f"Ground truth sites: {len(gt_sites)}")

# Run calibration sweep
thresholds = list(range(args.threshold_min, args.threshold_max + 1, 5))
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
    print(f"Threshold: {threshold}% | TP: {TP} | FP: {FP} | FN: {FN} | Precision: {precision:.3f} | Recall: {recall:.3f} | F1: {f1:.3f}")

df = pd.DataFrame(results)

# For IVT negative control, F1 based "best threshold" is meaningless since TP is always 0
# Report FP count trend instead
if len(gt_sites) == 0:
    print("\nIVT negative control mode: no ground truth sites, F1 is not meaningful.")
    print("FP counts across thresholds:")
    print(df[["threshold", "FP"]].to_string(index=False))
else:
    best = df.loc[df["F1"].idxmax()]
    print(f"\nBest threshold: {best['threshold']:.0f}% | F1: {best['F1']:.3f} | Precision: {best['Precision']:.3f} | Recall: {best['Recall']:.3f}")

# Plot 1: All metrics on one graph
plt.figure(figsize=(10, 6))
plt.plot(df["threshold"], df["Precision"], label="Precision", marker="o", markersize=4)
plt.plot(df["threshold"], df["Recall"], label="Recall", marker="o", markersize=4)
plt.plot(df["threshold"], df["F1"], label="F1 Score", marker="o", markersize=4)
plt.xlabel("Frequency threshold (%)")
plt.ylabel("Score")
plt.title("Calibration sweep: Precision, Recall and F1 Score vs frequency threshold")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("calibration_sweep_ivt_combined.png", dpi=150)

# Plot 2: FP count trend, useful for IVT negative control
plt.figure(figsize=(10, 6))
plt.plot(df["threshold"], df["FP"], color="red", marker="o", markersize=4)
plt.xlabel("Frequency threshold (%)")
plt.ylabel("False positive count")
plt.title("False positive count vs frequency threshold (IVT negative control)")
plt.grid(True)
plt.tight_layout()
plt.savefig("calibration_sweep_ivt_fp_trend.png", dpi=150)

df.to_csv("calibration_sweep_ivt_results.csv", index=False)
print("Plots saved as calibration_sweep_ivt_combined.png and calibration_sweep_ivt_fp_trend.png")
print("Results saved as calibration_sweep_ivt_results.csv")
