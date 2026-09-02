"""
Standalone test script for calibration_sweep.py Option 2 changes.
Tests two things independently:
1. Coverage filtering logic
2. Empty ground truth file handling logic
Does NOT modify or import from calibration_sweep.py.
"""

import pandas as pd
import os
import io

# ---------- Part 1: Test coverage filtering ----------

def apply_coverage_filter(df, min_coverage):
    """Filter dataframe rows by minimum coverage threshold."""
    return df[df["coverage"] >= min_coverage]


def test_coverage_filter():
    print("=== Test 1: Coverage filter ===")
    # synthetic bedMethyl-like data
    data = {
        "chrom": ["chr1", "chr1", "chr1", "chr1"],
        "start": [100, 200, 300, 400],
        "strand": ["+", "+", "+", "+"],
        "coverage": [5, 15, 9, 100],
        "freq": [90, 85, 95, 82],
    }
    df = pd.DataFrame(data)
    print("Before filter:")
    print(df)

    filtered = apply_coverage_filter(df, min_coverage=10)
    print("After filter (min_coverage=10):")
    print(filtered)

    expected_rows = 2  # only coverage 15 and 100 should survive
    assert len(filtered) == expected_rows, f"Expected {expected_rows} rows, got {len(filtered)}"
    print("PASS: coverage filter works as expected\n")


# ---------- Part 2: Test empty ground truth handling ----------

def read_ground_truth_safe(gt_path, columns):
    """
    Read ground truth file. If file is empty, return an empty dataframe
    with the same columns instead of raising an error.
    """
    if os.path.getsize(gt_path) == 0:
        print(f"Ground truth file is empty: {gt_path}")
        return pd.DataFrame(columns=columns)
    try:
        return pd.read_csv(gt_path, sep="\t", header=None, names=columns)
    except pd.errors.EmptyDataError:
        print(f"Ground truth file has no parseable data: {gt_path}")
        return pd.DataFrame(columns=columns)


def test_empty_ground_truth():
    print("=== Test 2: Empty ground truth handling ===")
    columns = ["chrom", "start", "strand"]

    # create a real empty temp file to simulate IVT control ground truth
    empty_path = "test_empty_ground_truth.tsv"
    with open(empty_path, "w") as f:
        pass  # write nothing, truly empty file

    gt_df = read_ground_truth_safe(empty_path, columns)
    print("Resulting ground truth dataframe:")
    print(gt_df)
    assert len(gt_df) == 0, "Expected empty ground truth dataframe"
    print("PASS: empty file handled without crash\n")

    os.remove(empty_path)

    # now test that TP/FP/FN logic behaves correctly with empty ground truth
    print("=== Test 3: TP/FP/FN with empty ground truth ===")
    called_sites = {"chr1_100_+", "chr1_200_+", "chr1_300_+"}
    gt_sites = set(gt_df["chrom"].astype(str) + "_" + gt_df["start"].astype(str) + "_" + gt_df["strand"].astype(str))

    TP = len(called_sites & gt_sites)
    FP = len(called_sites - gt_sites)
    FN = len(gt_sites - called_sites)

    print(f"TP: {TP}, FP: {FP}, FN: {FN}")
    assert TP == 0, "TP should be 0 when ground truth is empty"
    assert FP == len(called_sites), "All called sites should be FP when ground truth is empty"
    assert FN == 0, "FN should be 0 when ground truth is empty"
    print("PASS: all called sites correctly counted as FP\n")


if __name__ == "__main__":
    test_coverage_filter()
    test_empty_ground_truth()
    print("All tests passed.")