import decimal
import json
import os.path as path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame

with open("configuration.json", mode='r') as f:
    config = json.load(f)

Decimal = decimal.Decimal

decimal.getcontext().prec = 100

RESULT_FOLDER_PATH: str = r"C:\Users\a.pesterev\CLionProjects\porous_media\cmake-build-release-visual-studio\results"

CONCENTRATION: str = path.join(RESULT_FOLDER_PATH, "conc.dat")
REST: str = path.join(RESULT_FOLDER_PATH, "rest.dat")

c_df: DataFrame = pd.read_csv(CONCENTRATION, sep=';', skiprows=0, header=None, names=['ix', 'iy', 'x', 'y', 'c'],
                              index_col=False)
print(c_df)
s_df: DataFrame = pd.read_csv(REST, sep=';', skiprows=0, header=None, names=['ix', 'iy', 'x', 'y', 's'],
                              index_col=False)

MAX_T = config['max_t']
DT = config['dt']

print("Loaded")

lambda_p = config['delta_0']
k = config['F0']
S_max = config['s_m']

Sm = lambda_p * S_max / (lambda_p + k)


def Sasymp(x, t):
    return Decimal(Sm) * (Decimal(1) - Decimal(-(lambda_p + k) * Decimal(t - (1 + Sm) * x)).exp())


def Casymp(x, t):
    return Decimal(1) - (Decimal(1) - Decimal(-(lambda_p + k) * Sm * x).exp()) * Decimal(Decimal(-(lambda_p + k)) * Decimal(t - (1 + Sm) * x)).exp()


# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1, projection='3d')
# M_T = 1
# RDT = DT
# tos = np.arange(0., M_T, RDT)
# avs = list()
# for t0 in tos:
#     T = np.arange(t0, 4 * M_T, RDT)
#     T0 = t0 * np.ones(len(T))
#     X = T - t0
#     a = Sm - interp_s(X, T) / interp_c(X, T) - Sm * np.exp(-lambda_p * X + k * T) + np.exp(k * t0)
#     avs.append((T, a))
#     ax.scatter(T0, T, a, color='steelblue')
#     # if len(a) >= 4:
#     #     print(t0, a[-1], a[-2], a[-3], a[-4])
#     #     d1 = (-(a[-1] - a[-2]) / RDT)
#     #     d2 = (-(a[-2] - a[-3]) / RDT)
#     #     d3 = (-(a[-3] - a[-4]) / RDT)
#     #     print(
#     #         d1 / a[-1],
#     #         d2 / a[-2],
#     #         d3 / a[-3]
#     #     )
#     #     dd1 = (d1 - d2) / RDT
#     #     dd2 = (d2 - d3) / RDT
#     #     print(dd1, dd2)
#     #     print((dd1 - dd2) / RDT)
#
# ax.set_xlabel('t0')
# ax.set_ylabel('t')
# plt.show()
#
# tas = list()
# ras = list()
#
# for i in range(len(avs)):
#     tas.append(tos[i])
#     ras.append(avs[i][1][0])
#     # plt.cla()
#     # plt.clf()
#     # plt.plot(avs[i][0], avs[i][1])
#     # plt.title(f't0 = {tos[i]}')
#     # plt.show()
#
# plt.cla()
# plt.clf()
# plt.plot(tas, ras)
# plt.title(f't0 = t')
# plt.show()

def get_calced_C(ix, iy):
    return c_df.loc[(c_df['ix'] == ix) & (c_df['iy'] == iy)]['c'].values[0]


def get_calced_S(ix, iy):
    return s_df.loc[(s_df['ix'] == ix) & (s_df['iy'] == iy)]['s'].values[0]


tested_xs = [int(0.1 / DT), int(0.5 / DT), int(1. / DT), int(2. / DT), int(5. / DT), int(10. / DT)]

for i in range(len(tested_xs)):
    # cf = open(f"results_x={tested_xs[i] * DT}_dt={DT}.csv", "w")
    # cf.write("t;c;s;c_asymp;s_asymp;expc\n")
    dc = None
    ds = None
    expp = None
    ts = range(tested_xs[i], int(MAX_T / DT))
    vc = list()
    vs = list()
    c_asymp = list()
    s_asymp = list()
    dltc = list()
    dlts = list()
    v1 = list()
    v2 = list()
    for t in ts:
        expc = Decimal(2 * lambda_p * t * DT).exp()
        c_v = Decimal(float(get_calced_C(tested_xs[i], t)))
        s_v = Decimal(float(get_calced_S(tested_xs[i], t)))
        vc.append(c_v)
        vs.append(s_v)
        c_asymp_v = Decimal(Casymp(tested_xs[i] * DT, t * DT))
        c_asymp.append(c_asymp_v)
        s_asymp_v = Decimal(Sasymp(tested_xs[i] * DT, t * DT))
        s_asymp.append(s_asymp_v)
        v1.append(c_v - c_asymp_v)
        v2.append(s_v - s_asymp_v)
        dltc.append((c_v - c_asymp_v) * expc)
        dlts.append((s_v - s_asymp_v) * expc)
        # cf.write(f"{t * DT};{c_v};{s_v};{c_asymp_v};{s_asymp_v};{expc}\n")
    # cf.close()
    # plt.cla()
    # plt.clf()
    # plt.plot(ts, vc, label='numeric')
    # plt.plot(ts, c_asymp, label='asymp')
    # plt.legend()
    # plt.show()
    # plt.cla()
    # plt.clf()
    # plt.plot(ts, vs, label='numeric')
    # plt.plot(ts, s_asymp, label='asymp')
    # plt.legend()
    # plt.show()
    plt.cla()
    plt.clf()
    plt.title(f"x = {tested_xs[i] * DT}")
    plt.plot(np.array(ts) * DT, dltc, label='delta_C')
    plt.plot(np.array(ts) * DT, v1, label='delta_C2')
    plt.plot(np.array(ts) * DT, dlts, label='delta_S')
    plt.plot(np.array(ts) * DT, v2, label='delta_S2')
    # plt.plot(ts, characteristic, label='characteristic')
    plt.legend()
    plt.savefig(f"results_x={tested_xs[i] * DT}_dt={DT}.png")
