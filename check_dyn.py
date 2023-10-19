import json
import os
import os.path as path
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from mayavi import mlab
from pandas import DataFrame


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


CONFIGURATION_FILE_PATH: str = "./configuration.json"
RESULT_FOLDER_PATH: str = r".\cmake-build-release-visual-studio\results"

assert path.exists(CONFIGURATION_FILE_PATH), f"Configuration file path does not exists! {CONFIGURATION_FILE_PATH}"

with open(CONFIGURATION_FILE_PATH) as configuration_file:
    CONFIGURATION: dict = json.load(configuration_file)

CONCENTRATION: str = path.join(RESULT_FOLDER_PATH, "conc.dat")
REST: str = path.join(RESULT_FOLDER_PATH, "rest.dat")

c_df: DataFrame = pd.read_csv(CONCENTRATION, sep=';', skiprows=0, header=None, names=['ix', 'iy', 'x', 'y', 'c'],
                              index_col=False)
s_df: DataFrame = pd.read_csv(REST, sep=';', skiprows=0, header=None, names=['ix', 'iy', 'x', 'y', 's'],
                              index_col=False)

c_df = c_df.merge(s_df, on=['ix', 'iy', 'x', 'y'])

del s_df

c_df = c_df.sort_values(['ix', 'iy'])

c_df.replace("-nan(ind)", np.nan, inplace=True)
# c_df.c = c_df.c.astype(float).fillna(CONFIGURATION['C_0_t'])
# c_df.s = c_df.s.astype(float).fillna(CONFIGURATION['s_m'])

# Interpolating procedure
print(c_df.x, c_df.c + c_df.s)


def test_plot_iy(dataframe, iy, max_iy=1999):
    _ix = 0
    _x = np.zeros((max_iy - iy + 1))
    _c = np.zeros((max_iy - iy + 1))
    _s = np.zeros((max_iy - iy + 1))
    for i in range(iy, max_iy + 1):
        _df = dataframe.loc[(dataframe.ix == _ix) & (dataframe.iy == i)]
        _c[_ix] = _df.c
        _s[_ix] = _df.s
        _x[_ix] = _df.x
        _ix += 1

    plt.title(f"iy = {iy}")
    plt.plot(_x, _c)
    plt.plot(_x, _s)
    plt.show()


def test_plot_ix(dataframe, x):
    _ix = 0
    _df = dataframe[abs(dataframe.x - x) <= 1e-6]

    plt.title(f"x = {x}")
    plt.scatter(_df.y, _df.c)
    plt.scatter(_df.y, _df.s)
    plt.show()


# test_plot_ix(c_df, 0.2)
#
#
# for i in range(20, 41, 10):
#     test_plot_iy(c_df, i)


PLOT_3D = False
MPL = False

MAX_X = max(c_df.x.values)

DPI = 600

PLOTTING_ARGS_C = {
    "linestyle": "--",
    "color": "blue",
    "linewidth": 3.
}

PLOTTING_ARGS_S = {
    "linestyle": "-",
    "color": "red",
    "linewidth": 1.7
}

LEGEND_SETTINGS = {
    "fontsize": 'large'
}

print(f"MAX_X = {MAX_X}")

test_dots = np.stack([c_df.x.values, c_df.y.values], -1)

print("Got test dots")

with open("interp_c", mode='rb') as f:
    interp_c = pickle.load(f)

with open("interp_s", mode='rb') as f:
    interp_s = pickle.load(f)

print("Loaded")

# print(max(c_df.c.values - interp_c(test_dots)))

if PLOT_3D:
    if MPL:

        interpolated_c = interp_c(test_dots)
        interpolated_s = interp_s(test_dots)

        print("Interpolated")

        print(max(c_df.c.values - interpolated_c))
        print(max(c_df.s.values - interpolated_s))

        # fig = plt.figure()
        # ax = fig.add_subplot(1, 3, 1, projection='3d')
        # interp_ax = fig.add_subplot(1, 3, 2, projection='3d')
        # difference = fig.add_subplot(1, 3, 3, projection='3d')
        #
        # X = c_df.x.values
        # # print(X)
        # Y = c_df.y.values
        # # print(Y)
        # Z_C = c_df.c.values
        # # print(Z)
        #
        # # print(max(interpolated_c - Z))
        #
        # interp_ax.plot_trisurf(X, Y, interpolated_c, cmap=cm.jet, linewidth=0.2)
        # ax.plot_trisurf(X, Y, Z_C, cmap=cm.coolwarm, linewidth=0.2)
        # difference.plot_trisurf(X, Y, Z_C - interpolated_c, cmap=cm.bwr, linewidth=0.2)
        #
        # difference.set_xlabel('X')
        # difference.set_ylabel('T')
        # difference.set_zlabel('C - C(interp)')
        # interp_ax.set_xlabel('X')
        # interp_ax.set_ylabel('T')
        # interp_ax.set_zlabel('C(interp)')
        # ax.set_xlabel('X')
        # ax.set_ylabel('T')
        # ax.set_zlabel('C')
        # plt.show()
        #
        # plt.cla()
        # plt.clf()

        fig = plt.figure()
        ax = fig.add_subplot(1, 3, 1, projection='3d')
        interp_ax = fig.add_subplot(1, 3, 2, projection='3d')
        difference = fig.add_subplot(1, 3, 3, projection='3d')

        X = c_df.x.values
        # print(X)
        Y = c_df.y.values
        # print(Y)
        Z_C = c_df.s.values
        # print(Z)

        # print(max(interpolated_c - Z))

        interp_ax.plot_trisurf(X, Y, interpolated_s, cmap=cm.jet, linewidth=0.2)
        ax.plot_trisurf(X, Y, Z_C, cmap=cm.coolwarm, linewidth=0.2)
        difference.plot_trisurf(X, Y, Z_C - interpolated_s, cmap=cm.bwr, linewidth=0.2)

        difference.set_xlabel('X')
        difference.set_ylabel('T')
        difference.set_zlabel('S - S(interp)')
        interp_ax.set_xlabel('X')
        interp_ax.set_ylabel('T')
        interp_ax.set_zlabel('S(interp)')
        ax.set_xlabel('X')
        ax.set_ylabel('T')
        ax.set_zlabel('S')
        plt.show()

        plt.clf()
        plt.cla()
    else:
        X = c_df.x.values
        # print(X)
        Y = c_df.y.values
        # print(Y)
        Z_C = c_df.c.values
        # print(Z)
        Z_S = c_df.s.values

        # print(X, Y, Z_C)

        # print(max(interpolated_c - Z))

        # mlab.plot3d(X, Y, interpolated_c, interpolated_c, name='C(interp)', colormap='Spectral')
        mlab.plot3d(X, Y, Z_C, Z_C, name='C', colormap='Spectral')
        # mlab.plot3d(X, Y, Z_C - interpolated_c, 'C - C(interp)', colormap='Spectral')
        mlab.axes(xlabel='X', ylabel='Y', zlabel='Z', ranges=[min(X), max(X), min(Y), max(Y), min(Z_C), max(Z_C)])

        # mlab.plot3d(X, Y, interpolated_s, interpolated_s, name='S(interp)', colormap='Spectral')
        mlab.plot3d(X, Y, Z_S, Z_S, name='S', colormap='Spectral')
        # mlab.plot3d(X, Y, Z_S + Z_C, name='their_sum', colormap='Spectral')
        # mlab.plot3d(X, Y, interp_c(X, Y) - interp_s(X, Y) / 0.5, 'C - S', colormap='Spectral')
        # mlab.plot3d(X, Y, Z_S - interpolated_s, 'C - C(interp)', colormap='Spectral')
        mlab.axes(xlabel='X', ylabel='Y', zlabel='Z', ranges=[min(X), max(X), min(Y), max(Y), min(Z_S), max(Z_S)])
        mlab.show()

subplot_number = 111
XMAX_T = CONFIGURATION['max_t']
plotted_t = list(np.arange(0, XMAX_T, 0.1))

XMIN_T = 0
YMIN = 0
YMAX = 1

os.makedirs('images/t/', exist_ok=True)

for plotted_t_v in plotted_t:
    print(plotted_t_v)
    ax: plt.Axes = plt.subplot(subplot_number)
    condition = c_df.y == plotted_t_v
    inner_df = c_df[condition]
    # print(inner_df)
    x_v = inner_df.x.values.tolist()
    if len(x_v) != 0:
        max_x_v = max(x_v)
        c_v = list(inner_df.c.values)
        s_v = list(inner_df.s.values)
        plt.title(f"t = {str(plotted_t_v)}")

        x_zeros = np.arange(max_x_v, MAX_X, 0.1)

        zeros_y = np.zeros(len(x_zeros)).tolist()

        x_v += x_zeros.tolist()
        c_v += zeros_y
        s_v += zeros_y

        plt.xlim((XMIN_T, XMAX_T))
        plt.ylim((YMIN, YMAX))
        plt.plot(x_v, c_v, label='Suspended concentration C', **PLOTTING_ARGS_C)
        plt.plot(x_v, s_v, label='Retained concentration S', **PLOTTING_ARGS_S)

        plt.ylabel("C,S", rotation=0, y=1, horizontalalignment='right')
        plt.xlabel("x", x=1, horizontalalignment='right')

        plt.legend(**LEGEND_SETTINGS)
        fig = plt.gcf()
        fig.set_size_inches(6, 6)

        plt.savefig(f"./images/t/t({plotted_t.index(plotted_t_v):02d}) = {str(plotted_t_v)}.tiff", dpi=DPI)

        plt.clf()
        plt.cla()

MAX_T = CONFIGURATION['max_t']
dt = CONFIGURATION['dt']

t = np.arange(0, MAX_T, dt)

# for erosion!
plotted_xs = list(np.arange(0, MAX_X, 0.1))
# no erosion
# plotted_xs = [0.1, 0.25, 0.35, 0.5, 0.65, 0.75, 0.85, 1.]

initial = c_df[c_df.ix == 0]

XMIN_X = 0
XMAX_X = CONFIGURATION['max_t']

os.makedirs('images/x/', exist_ok=True)

subplot_number = 111
ax: plt.Axes = plt.subplot(subplot_number)
plt.title(f"x = 0")
ax.plot(initial.y, initial.c, label='Suspended concentration C', **PLOTTING_ARGS_C)
ax.plot(initial.y, initial.s, label='Retained concentration S', **PLOTTING_ARGS_S)
plt.legend(**LEGEND_SETTINGS)
plt.xlim((XMIN_X, XMAX_X))
plt.ylim((YMIN, YMAX))
plt.ylabel("C,S", rotation=0, y=1., horizontalalignment='right')
plt.xlabel("t", x=1., horizontalalignment='right')
fig = plt.gcf()
fig.set_size_inches(6, 6)
plt.savefig("./images/x/x(00) = 0.tiff", dpi=DPI)
plt.cla()

DELTA_0 = CONFIGURATION['delta_0']
S_M = CONFIGURATION['s_m']
F0 = CONFIGURATION['F0']

print(len(t))


def ds_dt(S, C):
    def Delta(S):
        return (DELTA_0 * (S_M - S))**0.5

    def F(S):
        return F0 * S

    return Delta(S) * C - F(S)


for plotted_x in plotted_xs:
    subplot_number = 111
    ax: plt.Axes = plt.subplot(subplot_number)
    plotted_c = list()
    plotted_s = list()
    print(plotted_x)
    for t_dot in range(len(t)):
        # Временной срез
        time_slice = c_df[c_df.iy == t_dot]
        if len(time_slice.x.values) == 0:
            m_x = -1
        else:
            m_x = max(time_slice.x.values)
        if m_x >= plotted_x:
            print(m_x, t[t_dot])
            from_t = t[t_dot:]
            with_x = plotted_x * np.ones(len(from_t))
            dots = np.stack([with_x, from_t], -1)
            for c_calc in interp_c(dots):
                plotted_c.append(c_calc)
            for s_calc in interp_s(dots):
                plotted_s.append(s_calc)
            break
        else:
            plotted_c.append(0.)
            plotted_s.append(0.)

    print(np.stack([t, plotted_c, plotted_s], -1))
    fig = plt.gcf()
    fig.set_size_inches(6, 6)
    plt.title(f"x = {plotted_x}")
    plt.xlim((XMIN_X, XMAX_X))
    plt.ylim((YMIN, YMAX))
    plt.ylabel("C,S", rotation=0, y=1., horizontalalignment='right')
    plt.xlabel("t", x=1., horizontalalignment='right')
    ax.plot(t, plotted_c, label='Suspended concentration C', **PLOTTING_ARGS_C)
    ax.plot(t, plotted_s, label='Retained concentration S', **PLOTTING_ARGS_S)
    plt.legend(**LEGEND_SETTINGS)
    plt.savefig(f"./images/x/x({(plotted_xs.index(plotted_x) + 1):02d}) = {plotted_x}.tiff", dpi=DPI)
    # plt.show()
    plt.cla()

# condition = c_df.ix == plotted_x
# args: dict[str, str | object] = dict(
#     kind="line",
#     x="y",
#     title=f"При x = {str(plotted_x)}",
#     ax=ax,
#     xlabel="t"
# )
# plot_with_condition(c_df, condition, label='Вычисленная концентрация', y="c", linestyle='-', **args)
# plot_with_condition(c_df, condition, y="s", label='Вычисленный осадок', linestyle='-', **args)
#
# plt.show()
