import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mrcfile


STAR_FILE = # path to the .star file 
MICROGRAPH_COL = 'rlnMicrographName'
TUBE_ID_COL = 'rlnHelicalTubeID'
COORD_X_COL = 'rlnCoordinateX'
COORD_Y_COL = 'rlnCoordinateY'
CLASS_COL = 'rlnClassNumber'
DOMINANT_CLASS_NUM = 1
MINOR_CLASS_NUM = 2
MIN_PARTICLES = 5

TARGET_MICROGRAPH = # the name of the micrograph.mrc  


MRC_FOLDER = # folder where the above micrograph is located 

DOMINANT_COLOR = # colour for plotting the dominant conformer
MINOR_COLOR = # colour for plotting the minor conformer

def read_star_particles(path):
    with open(path) as f:
        lines = f.readlines()

    header_cols = []
    data_start = None
    in_loop = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('loop_'):
            in_loop = True
            header_cols = []
            continue
        if in_loop and stripped.startswith('_rln'):
            col_name = stripped.split()[0][1:]
            header_cols.append(col_name)
            continue
        if in_loop and stripped and not stripped.startswith('_') and not stripped.startswith('#'):
            if MICROGRAPH_COL in header_cols and CLASS_COL in header_cols:
                data_start = i
                break
            else:
                in_loop = False

    if data_start is None:
        raise ValueError('Could not find a data block containing the required columns.')

    rows = []
    for line in lines[data_start:]:
        stripped = line.strip()
        if not stripped or stripped.startswith('data_') or stripped.startswith('loop_'):
            break
        rows.append(stripped.split())

    return pd.DataFrame(rows, columns=header_cols)


print('Reading star file...')
df = read_star_particles(STAR_FILE)

df[COORD_X_COL] = pd.to_numeric(df[COORD_X_COL], errors='coerce')
df[COORD_Y_COL] = pd.to_numeric(df[COORD_Y_COL], errors='coerce')
df[CLASS_COL] = pd.to_numeric(df[CLASS_COL], errors='coerce')
df['fibril_id'] = df[MICROGRAPH_COL].astype(str) + '_' + df[TUBE_ID_COL].astype(str)

fibril_counts = df.groupby('fibril_id')[CLASS_COL].agg(
    dominant_count=lambda x: (x == DOMINANT_CLASS_NUM).sum(),
    minor_count=lambda x: (x == MINOR_CLASS_NUM).sum()
).reset_index()
fibril_counts['is_mixed'] = ((fibril_counts['dominant_count'] >= MIN_PARTICLES) &
                              (fibril_counts['minor_count'] >= MIN_PARTICLES))
mixed_fibril_ids = set(fibril_counts.loc[fibril_counts['is_mixed'], 'fibril_id'])
print(f'Mixed fibrils in full dataset: {len(mixed_fibril_ids)}')

mic_df = df[df[MICROGRAPH_COL].str.contains(TARGET_MICROGRAPH, na=False)].copy()
mic_df = mic_df[mic_df['fibril_id'].isin(mixed_fibril_ids)]
mic_df = mic_df[mic_df[CLASS_COL].isin([DOMINANT_CLASS_NUM, MINOR_CLASS_NUM])]

if len(mic_df) == 0:
    raise ValueError(f'No mixed-fibril dominant/minor particles found on micrograph "{TARGET_MICROGRAPH}". '
                      f'Try a different micrograph, or check MIN_PARTICLES / class numbers.')

print(f'Particles plotted on this micrograph (mixed fibrils, dominant/minor only): {len(mic_df)}')

mrc_path = MRC_FOLDER + TARGET_MICROGRAPH.split('/')[-1]
print(f'Loading micrograph image: {mrc_path}')

with mrcfile.open(mrc_path, permissive=True) as mrc:
    image = mrc.data.astype(np.float32)
    if image.ndim == 3:
        image = image[0]

dominant_df = mic_df[mic_df[CLASS_COL] == DOMINANT_CLASS_NUM]
minor_df = mic_df[mic_df[CLASS_COL] == MINOR_CLASS_NUM]

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 15

plt.figure(figsize=(5, 5))
vmin, vmax = np.percentile(image, [1, 99])
plt.imshow(image, cmap='gray', vmin=vmin, vmax=vmax, origin='upper')

plt.scatter(dominant_df[COORD_X_COL], dominant_df[COORD_Y_COL],
            c=DOMINANT_COLOR, marker='o', s=40, edgecolors='black', linewidths=0.3,
            label='Dominant')
plt.scatter(minor_df[COORD_X_COL], minor_df[COORD_Y_COL],
            c=MINOR_COLOR, marker='^', s=40, edgecolors='black', linewidths=0.3,
            label='Minor')

plt.axis('off')
plt.legend(loc='upper right', fontsize=20, frameon=True)
plt.tight_layout()
plt.savefig('micrograph_particle_overlay_full.png', dpi=4000)
plt.show()

print('\nSaved: micrograph_particle_overlay_full.png')
