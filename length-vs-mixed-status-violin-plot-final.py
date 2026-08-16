import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

INPUT_FILE = # the path to the csv file we generated with the lengh-versus-mixed-state.py script uploaded on github called fibril_length_vs_mixed.csv
MIXED_COLOR = '#3fae7a'
NONMIXED_COLOR = '#9c9c9c'


df = pd.read_csv(INPUT_FILE)

mixed_lengths = df.loc[df['is_mixed'], 'total_particles'].values
nonmixed_lengths = df.loc[~df['is_mixed'], 'total_particles'].values
stat, p_value = mannwhitneyu(mixed_lengths, nonmixed_lengths, alternative='two-sided')
p_label = '<0.0001' if p_value < 0.0001 else f'{p_value:.4g}'

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 15
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'

fig, ax = plt.subplots(figsize=(4.5, 6))

parts = ax.violinplot([mixed_lengths, nonmixed_lengths], positions=[1, 2],
                       showmeans=False, showmedians=False, showextrema=False, widths=0.7)

violin_colors = [MIXED_COLOR, NONMIXED_COLOR]
for pc, color in zip(parts['bodies'], violin_colors):
    pc.set_facecolor(color)
    pc.set_edgecolor(color)
    pc.set_linewidth(0.5)
    pc.set_alpha(0.9)

mean_mixed = np.mean(mixed_lengths)
mean_nonmixed = np.mean(nonmixed_lengths)
ax.hlines(mean_mixed, 1 - 0.25, 1 + 0.25, color='black', linewidth=2.5, zorder=5)
ax.hlines(mean_nonmixed, 2 - 0.25, 2 + 0.25, color='black', linewidth=2.5, zorder=5)

y_max = max(mixed_lengths.max(), nonmixed_lengths.max())
bracket_y = y_max * 1.05
ax.plot([1, 1, 2, 2], [bracket_y, bracket_y * 1.02, bracket_y * 1.02, bracket_y],
        color='black', linewidth=1.5)
ax.text(1.5, bracket_y * 1.03, p_label, ha='center', va='bottom', fontsize=15, fontweight='bold')

ax.set_xticks([1, 2])
ax.set_xticklabels(['Mixed', 'Non-mixed'])
ax.set_ylabel('Fibril length (particle count)', fontweight='bold')
ax.set_ylim(0, bracket_y * 1.15)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params(width=2)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontweight('bold')

plt.tight_layout()
plt.savefig('fibril_length_scatter.png', dpi=300)
plt.show()
print('Saved: fibril_length_scatter.png')
print(f'Mann-Whitney p-value: {p_value:.4g}')
print(f'\nMean fibril length (particle count):')
print(f'  Mixed:     {mean_mixed:.2f}')
print(f'  Non-mixed: {mean_nonmixed:.2f}')
print(f'\nMedian fibril length (particle count):')
print(f'  Mixed:     {np.median(mixed_lengths):.2f}')
print(f'  Non-mixed: {np.median(nonmixed_lengths):.2f}')

