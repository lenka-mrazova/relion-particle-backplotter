import pandas as pd


INPUT_FILE = #path to the csv file which has all fibrils classified as dominant only, minor only, mixed or neither
MICROGRAPH_COL = 'rlnMicrographName'
TUBE_ID_COL = 'rlnHelicalTubeID'
DOMINANT_CLASSES = ['Class_1']
MINOR_CLASSES = ['Class_2']
MIN_PARTICLES = 5


df = pd.read_csv(INPUT_FILE)
df['fibril_id'] = df[MICROGRAPH_COL].astype(str) + '_' + df[TUBE_ID_COL].astype(str)
df['dominant_count'] = df[DOMINANT_CLASSES].sum(axis=1)
df['minor_count'] = df[MINOR_CLASSES].sum(axis=1)

total_dominant = df['dominant_count'].sum()
total_minor = df['minor_count'].sum()
total_both = total_dominant + total_minor

print('=== Dataset-wide conformer proportions (all Class_1/Class_2 particles) ===')
print(f'Dominant particles: {total_dominant}')
print(f'Minor particles:    {total_minor}')
print(f'Dominant: {100*total_dominant/total_both:.1f}%   Minor: {100*total_minor/total_both:.1f}%')
print(f'Ratio (dominant:minor): {total_dominant/total_minor:.1f} : 1\n')


df['is_mixed'] = ((df['dominant_count'] >= MIN_PARTICLES) &
                   (df['minor_count'] >= MIN_PARTICLES))
mixed = df[df['is_mixed']].copy()

mixed['dom_fraction'] = mixed['dominant_count'] / (mixed['dominant_count'] + mixed['minor_count'])

print(f'=== Balance within mixed fibrils (n = {len(mixed)}) ===')
print(f'Mean dominant fraction per mixed fibril:   {100*mixed["dom_fraction"].mean():.1f}%')
print(f'Median dominant fraction per mixed fibril: {100*mixed["dom_fraction"].median():.1f}%')
print(f'  (mean minor fraction ~ {100*(1-mixed["dom_fraction"].mean()):.1f}%)\n')


majority_dom = (mixed['dom_fraction'] > 0.6).sum()
majority_min = (mixed['dom_fraction'] < 0.4).sum()
roughly_even = ((mixed['dom_fraction'] >= 0.4) & (mixed['dom_fraction'] <= 0.6)).sum()

print(f'Of the {len(mixed)} mixed fibrils:')
print(f'  Majority dominant (>60% dominant): {majority_dom} ({100*majority_dom/len(mixed):.1f}%)')
print(f'  Roughly even (40-60% dominant):    {roughly_even} ({100*roughly_even/len(mixed):.1f}%)')
print(f'  Majority minor (<40% dominant):    {majority_min} ({100*majority_min/len(mixed):.1f}%)')
