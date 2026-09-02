#!/usr/bin/env python3
"""Generate thesis-quality IVT false positive trend plot."""

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.size'] = 14

thresholds = [5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100]
fp_counts = [46390,10492,5446,2076,1204,652,432,230,171,102,77,54,40,27,18,14,14,9,6,0]

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(thresholds, fp_counts, 'o-', color='#D62728', markersize=5, linewidth=1.5)
ax.axvline(x=15, color='grey', linestyle='--', linewidth=1, label='15% threshold')
ax.set_xlabel('Frequency threshold (%)')
ax.set_ylabel('False positive count')
ax.legend(frameon=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('/data/projects/punim2974/master-project-yenan/figures/fp_trend_ivt_thesis.png', dpi=300, bbox_inches='tight')
plt.savefig('/data/projects/punim2974/master-project-yenan/figures/fp_trend_ivt_thesis.pdf', bbox_inches='tight')
print('Saved to figures/fp_trend_ivt_thesis.png and .pdf')
