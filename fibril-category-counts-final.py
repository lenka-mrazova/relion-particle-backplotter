
import pandas as pd
import matplotlib.pyplot as plt


INPUT_FILE = # path to a csv file that counts the number of particles per fibril that go into each 3D class
MICROGRAPH_COL = 'rlnMicrographName'
TUBE_ID_COL = 'rlnHelicalTubeID'
DOMINANT_CLASSES = ['Class_1']
MINOR_CLASSES = ['Class_2']
MIN_PARTICLES = 5

DOMINANT_COLOR = '#65c8e3'
MINOR_COLOR = '#81ce6d'
MIXED_COLOR = '#3fae7a'   
NEITHER_COLOR = '#bfbfbf'


df = pd.read_csv(INPUT_FILE)
df['fibril_id'] = df[MICROGRAPH_COL].astype(str) + '_' + df[TUBE_ID_COL].astype(str)
df['dominant_count'] = df[DOMINANT_CLASSES].sum(axis=1)
df['minor_count'] = df[MINOR_CLASSES].sum(axis=1)

def classify_fibril(row):
    has_dom = row['dominant_count'] >= MIN_PARTICLES
    has_min = row['minor_count'] >= MIN_PARTICLES
    if has_dom and has_min:
        return 'Mixed'
    elif has_dom:
        return 'Dominant only'
    elif has_min:
        return 'Minor only'
    else:
        return 'Neither'

df['category'] = df.apply(classify_fibril, axis=1)
n_total = len(df)
counts = df['category'].value_counts()

categories = ['Mixed', 'Minor only', 'Dominant only']
colors = [MIXED_COLOR, MINOR_COLOR, DOMINANT_COLOR]
percentages = [100 * counts.get(c, 0) / n_total for c in categories]

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 18
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'

fig, ax = plt.subplots(figsize=(5, 5))
bars = ax.barh(categories, percentages, color=colors, edgecolor='black', linewidth=0)

ax.set_xlabel('Fibrils (%)', fontweight='bold')
ax.set_xlim(0, 30)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params(width=2)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontweight('bold')

plt.tight_layout()
plt.savefig('fibril_category_barh.png', dpi=2000)
plt.show()
print('Saved: fibril_category_barh.png')
for c, p in zip(categories, percentages):
    print(f'{c}: {p:.2f}%')
