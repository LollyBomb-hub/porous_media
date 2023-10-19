from scipy.interpolate import CloughTocher2DInterpolator
import pickle
import pandas as pd
from pandas import DataFrame
import os.path as path
import json
import numpy as np

CONFIGURATION_FILE_PATH: str = "./configuration.json"
RESULT_FOLDER_PATH: str = r".\cmake-build-release-visual-studio\results"

CONCENTRATION: str = path.join(RESULT_FOLDER_PATH, "conc.dat")
REST: str = path.join(RESULT_FOLDER_PATH, "rest.dat")

c_df: DataFrame = pd.read_csv(CONCENTRATION, sep=';', skiprows=0, header=None, names=['ix', 'iy', 'x', 'y', 'c'],
                              index_col=False)
s_df: DataFrame = pd.read_csv(REST, sep=';', skiprows=0, header=None, names=['ix', 'iy', 'x', 'y', 's'],
                              index_col=False)

c_df = c_df.merge(s_df, on=['ix', 'iy', 'x', 'y'])

del s_df

c_df.to_csv(path_or_buf='results.csv', sep=';', header=True)


assert path.exists(CONFIGURATION_FILE_PATH), f"Configuration file path does not exists! {CONFIGURATION_FILE_PATH}"

with open(CONFIGURATION_FILE_PATH) as configuration_file:
    CONFIGURATION: dict = json.load(configuration_file)

c_df = c_df.sort_values(['ix', 'iy'])
c_df.replace("-nan(ind)", np.nan, inplace=True)
c_df.c = c_df.c.astype(float).fillna(-CONFIGURATION['C_0_t'])
c_df.s = c_df.s.astype(float).fillna(-CONFIGURATION['s_m'])
print(c_df)

interp_c = CloughTocher2DInterpolator(np.stack([c_df.x.values, c_df.y.values], -1), c_df.c.values)
interp_s = CloughTocher2DInterpolator(np.stack([c_df.x.values, c_df.y.values], -1), c_df.s.values)

print("Dumping")

pickle.dump(interp_c, open("interp_c", mode='wb'))
pickle.dump(interp_s, open("interp_s", mode='wb'))
