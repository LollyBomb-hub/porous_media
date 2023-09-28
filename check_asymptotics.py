import json
import pickle

import matplotlib.pyplot as plt
import numpy as np

with open("configuration.json", mode='r') as f:
    config = json.load(f)

with open("interp_c", mode='rb') as f:
    interp_c = pickle.load(f)

with open("interp_s", mode='rb') as f:
    interp_s = pickle.load(f)

MAX_T = config['max_t']
DT = config['dt']

print("Loaded")

lambda_p = config['delta_0']
k = config['F0']
S_max = config['s_m']

Sm = lambda_p * S_max / (lambda_p + k)


def Sasymp(x, t):
    return Sm * (1 - np.exp(-(lambda_p + k) * (t - (1 + Sm) * x)))


def Casymp(x, t):
    return 1 - (1 - np.exp(-(lambda_p + k) * Sm * x)) * np.exp(-(lambda_p + k) * (t - (1 + Sm) * x))


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


tested_xs = [0.1, 0.5, 1., 2., 5., 10.]

for i in range(len(tested_xs)):
    dc = None
    ds = None
    expp = None
    ts = np.arange(tested_xs[i], MAX_T, DT)
    vc = list()
    vs = list()
    c_asymp = list()
    s_asymp = list()
    dltc = list()
    dlts = list()
    v1 = list()
    v2 = list()
    for t in ts:
        expc = float(np.exp(2 * lambda_p * t))
        c_v = float(interp_c(tested_xs[i], t))
        s_v = float(interp_s(tested_xs[i], t))
        vc.append(c_v)
        vs.append(s_v)
        c_asymp_v = float(Casymp(tested_xs[i], t))
        c_asymp.append(c_asymp_v)
        s_asymp_v = float(Sasymp(tested_xs[i], t))
        s_asymp.append(s_asymp_v)
        v1.append(c_v - c_asymp_v)
        v2.append(s_v - s_asymp_v)
        if expp is None:
            expp = expc
        if dc is None or dc == 0.:
            dc = (c_v - c_asymp_v)
            dltc.append(dc * expc)
        else:
            k2 = expc / expp
            cd = c_v - c_asymp_v
            k1 = cd / dc
            dc_v = dltc[-1] * k1 * k2
            dltc.append(dc_v)
            dc = cd
        if ds is None or ds == 0.:
            ds = (s_v - s_asymp_v)
            dlts.append(ds * expc)
        else:
            k2 = expc / expp
            sd = s_v - s_asymp_v
            k1 = sd / ds
            ds_v = dlts[-1] * k1 * k2
            dlts.append(ds_v)
            ds = sd
    plt.cla()
    plt.clf()
    plt.plot(ts, vc, label='numeric')
    plt.plot(ts, c_asymp, label='asymp')
    plt.legend()
    plt.show()
    plt.cla()
    plt.clf()
    plt.plot(ts, vs, label='numeric')
    plt.plot(ts, s_asymp, label='asymp')
    plt.legend()
    plt.show()
    plt.cla()
    plt.clf()
    # plt.plot(ts, dltc, label='delta_C')
    plt.plot(ts, v1, label='delta_C2')
    # plt.plot(ts, dlts, label='delta_S')
    plt.plot(ts, v2, label='delta_S2')
    # plt.plot(ts, characteristic, label='characteristic')
    plt.legend()
    plt.show()
