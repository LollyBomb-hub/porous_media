import contextlib
import glob
import json
import os
import os.path as path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from pandas import DataFrame

CONFIGURATION_FILE_PATH: str = "C:/Users/a.pesterev/CLionProjects/porous_media/configuration.json"
RESULT_FOLDER_PATH: str = r"C:\Users\a.pesterev\CLionProjects\porous_media\cmake-build-release-visual-studio\results"

assert path.exists(CONFIGURATION_FILE_PATH), f"Configuration file path does not exists! {CONFIGURATION_FILE_PATH}"

with open(CONFIGURATION_FILE_PATH) as configuration_file:
    CONFIGURATION: dict = json.load(configuration_file)

RESULT_FILE: str = path.join(RESULT_FOLDER_PATH, "result101_400001.csv")

df: DataFrame = pd.read_csv(RESULT_FILE, sep=';', skiprows=4, index_col=False,
                            dtype={'c': np.float64, 's': np.float64, 'x': np.float64, 't': np.float64}, na_values='-nan(ind)')

ANALYTIC_SOLUTION_IS_KNOWN = False


def plot_with_condition(passed_df: pd.DataFrame, condition_instance: object, **kwargs: object) -> None:
    """
    :param kwargs:
    :rtype: None
    :param passed_df:
    :type condition_instance: object
    """
    fixed_x_df = passed_df[condition_instance]
    # print(fixed_x_df)
    fixed_x_df.plot(**kwargs)


def make_gif_from_pngs(png_files_regex: str | list[str], output_file: str, duration: int = 200):
    # use exit stack to automatically close opened images
    with contextlib.ExitStack() as stack:
        # lazily load images
        if png_files_regex is str:
            imgs = (stack.enter_context(Image.open(f))
                    for f in glob.glob(png_files_regex))
        elif isinstance(png_files_regex, list):
            imgs = (stack.enter_context(Image.open(f))
                    for f in png_files_regex)
        else:
            raise ValueError(f"Unexpected type of param {png_files_regex}, {type(png_files_regex)}!")

        # extract  first image from iterator
        img = next(imgs)

        # https://pillow.readthedocs.io/en/stable/handbook/image-file-formats.html#gif
        img.save(fp=output_file, format='GIF', append_images=imgs,
                 save_all=True, duration=duration, loop=0)


mult = 500
RDX: float = CONFIGURATION['dx']
DX: float = CONFIGURATION['dx'] * mult
COUNT_X_DOTS: int = int(CONFIGURATION['N'] / mult)

print(df.columns)

print(df.head(5))

plotted_xs: list[int] = list(map(lambda el: el * mult, range(COUNT_X_DOTS)))

subplot_numbers: list[int] = [111 for _ in range(len(plotted_xs))]

C_COORDINATE = 'c(i, j)'

output_filenames: list[str] = list()

if not os.path.exists("./images/"):
    os.mkdir("./images")

plotted_x_i: int
subplot_number: int
for plotted_x_i, subplot_number in zip(plotted_xs, subplot_numbers):
    plotted_x: float = plotted_x_i * RDX
    ax: plt.Axes = plt.subplot(subplot_number)
    condition = df.i == plotted_x_i
    args: dict[str, str | object] = dict(
        kind="line",
        x="t",
        title=f"При x = {str(round(plotted_x, len(str(mult))))}",
        ax=ax,
        xlabel="Время, t"
    )
    plot_with_condition(df, condition, label='Вычисленная концентрация', y=C_COORDINATE, linestyle='-', **args)
    if ANALYTIC_SOLUTION_IS_KNOWN:
        plot_with_condition(df, condition, y="c_e", label='Аналитическое решение для концентрации', linestyle='-.',
                            **args)
    plot_with_condition(df, condition, y="s", label='Вычисленный осадок', linestyle='-', **args)
    if ANALYTIC_SOLUTION_IS_KNOWN:
        plot_with_condition(df, condition, y="s_e", label='Аналитическое решение для осадка', linestyle='-.', **args)
    plt.legend()
    filename: str = f'images/x={int(plotted_x_i / mult)}.png'
    output_filenames.append(filename)
    plt.savefig(filename)
    ax.cla()
    plt.clf()

make_gif_from_pngs(output_filenames, 'output200.gif')
make_gif_from_pngs(output_filenames, 'output50.gif', duration=50)

# plotted_ts: list[int | float] = [0.5, 1, 2, 5, 10, 20, 50, 100]
#
# for plotted_t, subplot_number in zip(plotted_ts, subplot_numbers):
#     ax: plt.Axes = plt.subplot(subplot_number)
#     condition = df.t == plotted_t
#     args: dict[str, str | object] = dict(
#         kind="line",
#         x="x",
#         title=f"При t = {str(plotted_t)}",
#         ax=ax,
#         xlabel="x"
#     )
#     plot_with_condition(df, condition, label='Вычисленная концентрация', y="c", linestyle='-', **args)
#     if ANALYTIC_SOLUTION_IS_KNOWN:
#         plot_with_condition(df, condition, y="c_e", label='Аналитическое решение для концентрации', linestyle='-.',
#                             **args)
#     plot_with_condition(df, condition, y="s", label='Вычисленный осадок', linestyle='-', **args)
#     if ANALYTIC_SOLUTION_IS_KNOWN:
#         plot_with_condition(df, condition, y="s_e", label='Аналитическое решение для осадка', linestyle='-.', **args)
#     plt.legend()
#     plt.show()
