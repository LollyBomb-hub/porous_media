from scipy.interpolate import CloughTocher2DInterpolator
import pickle
import pandas as pd
import os.path as path
import json
from pandas import DataFrame
import numpy as np

CONFIGURATION_FILE_PATH: str = "C:/Users/a.pesterev/CLionProjects/porous_media/configuration.json"
RESULT_FOLDER_PATH: str = r"C:\Users\a.pesterev\CLionProjects\porous_media\cmake-build-release-visual-studio\results"

assert path.exists(CONFIGURATION_FILE_PATH), f"Configuration file path does not exists! {CONFIGURATION_FILE_PATH}"

with open(CONFIGURATION_FILE_PATH) as configuration_file:
    CONFIGURATION: dict = json.load(configuration_file)

CONCENTRATION: str = path.join(RESULT_FOLDER_PATH, "conc.dat")
REST: str = path.join(RESULT_FOLDER_PATH, "rest.dat")

c_df: DataFrame = pd.read_csv(CONCENTRATION, sep=';', skiprows=0, header=None, names=['ix', 'iy', 'x', 'y', 'c'],
                              index_col=False)
print(c_df)
s_df: DataFrame = pd.read_csv(REST, sep=';', skiprows=0, header=None, names=['ix', 'iy', 'x', 'y', 's'],
                              index_col=False)

c_df = c_df.merge(s_df, on=['ix', 'iy', 'x', 'y'])

del s_df

c_df = c_df.sort_values(['ix', 'iy'])

interp_c = CloughTocher2DInterpolator(np.stack([c_df.x.values, c_df.y.values], -1), c_df.c.values)
interp_s = CloughTocher2DInterpolator(np.stack([c_df.x.values, c_df.y.values], -1), c_df.s.values)

print("Dumping")

pickle.dump(interp_c, open("interp_c", mode='wb'))
pickle.dump(interp_s, open("interp_s", mode='wb'))
