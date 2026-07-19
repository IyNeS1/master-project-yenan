#!/usr/bin/env python3
"""
Plot false positive trend for IVT negative control
and compare with WT calibration sweep results.
"""

import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

def main():
    if len(sys.argv) < 4:
        print("Usage: python fp_trend_plot.py <ivt_csv> <wt_csv> <output_dir>")
        sys.exit(1)

    ivt_csv = sys.argv[1]
    wt_csv = sys.argv[2]
    output_dir = sys.argv[3]

    os.makedirs(output_dir, exist_ok=True)

    ivt = pd.read_csv(ivt_csv)
    wt = pd.read_csv(wt_csv)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: IVT FP trend
    ax1 = axes[0]
    ax1.plot(ivt["threshold"], ivt["FP"], marker="o", color="#d62728", linewidth=2)
    ax1.set_xlabel("Modification Frequency Threshold (%)", fontsize=12)
    ax1.set_ylabel("Number of False Positives", fontsize=12)
    ax1.set_title("A. IVT Negative Control: FP Trend", fontsize=13, fontweight="bold")
    ax1.set_xticks(range(0, 105, 10))
    ax1.grid(True, alpha=0.3)
    ax1.axvspan(15, 20, alpha=0.15, color="orange", label="Elbow region (15-20%)")
    ax1.legend(fontsize=10)

    # Panel B: FP comparison IVT vs WT
    ax2 = axes[1]
    ax2.plot(ivt["threshold"], ivt["FP"], marker="o", color="#d62728",
             linewidth=2, label="IVT (negative control)")
    ax2.plot(wt["threshold"], wt["FP"], marker="s", color="#1f77b4",
             linewidth=2, label="WT (experimental)")
    ax2.set_xlabel("Modification Frequency Threshold (%)", fontsize=12)
    ax2.set_ylabel("Number of False Positives", fontsize=12)
    ax2.set_title("B. FP Comparison: IVT vs WT", fontsize=13, fontweight="bold")
    ax2.set_xticks(range(0, 105, 10))
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)

    plt.tight_layout()
    outpath = os.path.join(output_dir, "fp_trend_ivt_vs_wt.png")
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    print(f"Saved: {outpath}")
    plt.close()

if __name__ == "__main__":
    main()
