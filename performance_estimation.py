import pandas as pd
import argparse

# Accept command line arguments
parser = argparse.ArgumentParser(description="Performance estimation for m6A pipeline")
parser.add_argument("--m6a", required=True, help="Path to m6A sample modkit output .bed file")
parser.add_argument("--control", required=True, help="Path to control sample modkit output .bed file")
parser.add_argument("--ground_truth", required=True, help="Path to ground truth .bed file")
args = parser.parse_args()

# Read files
modkit_cols = ["chrom", "start", "end", "mod", "score", "strand",
               "start2", "end2", "color", "coverage", "freq",
               "mod_count", "unmod_count", "c1", "c2", "c3", "c4", "c5"]

m6a = pd.read_csv(args.m6a, sep="\t", header=None, names=modkit_cols)
control = pd.read_csv(args.control, sep="\t", header=None, names=modkit_cols)
ground_truth = pd.read_csv(args.ground_truth, sep="\t", header=None,
                            names=["chrom", "start", "end", "mod", "score", "strand"])

# Keep only m6A sites (mod == "a")
m6a = m6a[m6a["mod"] == "a"]
control = control[control["mod"] == "a"]
ground_truth = ground_truth[ground_truth["mod"] == "a"]

# Create unique site identifier using chrom + start + strand
m6a["site"] = m6a["chrom"].astype(str) + "_" + m6a["start"].astype(str) + "_" + m6a["strand"].fillna("").astype(str)
control["site"] = control["chrom"].astype(str) + "_" + control["start"].astype(str) + "_" + control["strand"].fillna("").astype(str)
ground_truth["site"] = ground_truth["chrom"].astype(str) + "_" + ground_truth["start"].astype(str) + "_" + ground_truth["strand"].fillna("").astype(str)

gt_sites = set(ground_truth["site"])
m6a_sites = set(m6a["site"])
control_sites = set(control["site"])

# Calculate TP, FP, FN for m6A sample
TP = len(m6a_sites & gt_sites)
FP = len(m6a_sites - gt_sites)
FN = len(gt_sites - m6a_sites)

precision = TP / (TP + FP) if (TP + FP) > 0 else 0
recall    = TP / (TP + FN) if (TP + FN) > 0 else 0
f1        = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

print("=== m6A sample (positive control) ===")
print(f"TP: {TP}, FP: {FP}, FN: {FN}")
print(f"Precision: {precision:.3f}")
print(f"Recall:    {recall:.3f}")
print(f"F1 Score:  {f1:.3f}")

# Control sample should ideally have no m6A sites detected
print("\n=== Control sample (negative control) ===")
print(f"Sites detected: {len(control_sites)}")
print(f"False positives (sites in ground truth): {len(control_sites & gt_sites)}")