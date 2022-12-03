import json
import os.path as path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame

CONFIGURATION_FILE_PATH: str = "C:/Users/a.pesterev/CLionProjects/porous_media/configuration.json"
RESULT_FOLDER_PATH: str = "C:/Users/a.pesterev/CLionProjects/porous_media/cmake-build-release/results"

assert path.exists(CONFIGURATION_FILE_PATH), f"Configuration file path does not exists! {CONFIGURATION_FILE_PATH}"

with open(CONFIGURATION_FILE_PATH) as configuration_file:
    CONFIGURATION: object = json.load(configuration_file)

RESULT_FILE: str = path.join(RESULT_FOLDER_PATH, "result.csv")

df: DataFrame = pd.read_csv(RESULT_FILE, sep=';', skiprows=4, index_col=False,
                            dtype={'c': np.float64, 's': np.float64, 'x': np.float64, 't': np.float64})

print(df)

ANALYTIC_SOLUTION_IS_KNOWN = False


def plot_with_condition(passed_df: pd.DataFrame, condition_instance: object, **kwargs: object) -> None:
    """
    :param kwargs:
    :rtype: None
    :param passed_df:
    :type condition_instance: object
    """
    fixed_x_df = passed_df[condition_instance]
    print(fixed_x_df)
    fixed_x_df.plot(**kwargs)


plotted_xs: list[int | float] = np.arange(0, 0.25, 0.01)
subplot_numbers: list[int] = [111 for _ in range(len(plotted_xs))]

plotted_x: int | float
subplot_number: int
for plotted_x, subplot_number in zip(plotted_xs, subplot_numbers):
    ax: plt.Axes = plt.subplot(subplot_number)
    condition = df.x == plotted_x
    args: dict[str, str | object] = dict(
        kind="line",
        x="t",
        title=f"При x = {str(plotted_x)}",
        ax=ax,
        xlabel="Время, t"
    )
    plot_with_condition(df, condition, label='Вычисленная концентрация', y="c", linestyle='-', **args)
    if ANALYTIC_SOLUTION_IS_KNOWN:
        plot_with_condition(df, condition, y="c_e", label='Аналитическое решение для концентрации', linestyle='-.',
                            **args)
    plot_with_condition(df, condition, y="s", label='Вычисленный осадок', linestyle='-', **args)
    if ANALYTIC_SOLUTION_IS_KNOWN:
        plot_with_condition(df, condition, y="s_e", label='Аналитическое решение для осадка', linestyle='-.', **args)
    plt.legend()
    plt.show()

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
