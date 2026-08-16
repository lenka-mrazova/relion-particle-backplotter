import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, pointbiserialr

INPUT_FILE = # the path to the csv that has the particle counts per 3D class for each fibril ID 
MICROGRAPH_COL = 'rlnMicrographName'
TUBE_ID_COL = 'rlnHelicalTubeID'
DOMINANT_CLASSES = ['Class_1']
MINOR_CLASSES = ['Class_2']
ALL_CLASS_COLS = [f'Class_{i}' for i in range(1, 13)]   # Class_1 .. Class_12
MIN_PARTICLES = 5


df = pd.read_csv(INPUT_FILE)
df['fibril_id'] = df[MICROGRAPH_COL].astype(str) + '_' + df[TUBE_ID_COL].astype(str)

df['dominant_count'] = df[DOMINANT_CLASSES].sum(axis=1)
df['minor_count'] = df[MINOR_CLASSES].sum(axis=1)

df['total_particles'] = df[ALL_CLASS_COLS].sum(axis=1)

df['is_mixed'] = ((df['dominant_count'] >= MIN_PARTICLES) &
                   (df['minor_count'] >= MIN_PARTICLES))

mixed_lengths = df.loc[df['is_mixed'], 'total_particles']
not_mixed_lengths = df.loc[~df['is_mixed'], 'total_particles']

print(f'Mixed fibrils: {df["is_mixed"].sum()}')
print(f'Non-mixed fibrils: {(~df["is_mixed"]).sum()}\n')

print('Fibril length (total particles) summary:')
print(f'  Mixed     - median: {mixed_lengths.median():.0f}, mean: {mixed_lengths.mean():.1f}')
print(f'  Non-mixed - median: {not_mixed_lengths.median():.0f}, mean: {not_mixed_lengths.mean():.1f}\n')

stat, p_value = mannwhitneyu(mixed_lengths, not_mixed_lengths, alternative='greater')
print(f'Mann-Whitney U test (H1: mixed fibrils are longer):')
print(f'  p-value: {p_value:.4g}')

corr, corr_p = pointbiserialr(df['is_mixed'].astype(int), df['total_particles'])
print(f'\nPoint-biserial correlation (mixed status vs. fibril length):')
print(f'  r = {corr:.3f}, p = {corr_p:.4g}')

df['length_quartile'] = pd.qcut(df['total_particles'], 4, labels=['Q1 (shortest)', 'Q2', 'Q3', 'Q4 (longest)'])
pct_mixed_by_quartile = df.groupby('length_quartile')['is_mixed'].mean() * 100
print(f'\n% of fibrils that are mixed, by length quartile:')
print(pct_mixed_by_quartile.round(2))

df.to_csv('fibril_length_vs_mixed.csv', index=False)
print('\nFull data saved to fibril_length_vs_mixed.csv')
