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


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
M_T = 1
RDT = DT
tos = np.arange(0., M_T, RDT)
avs = list()
for t0 in tos:
    T = np.arange(t0, 4 * M_T, RDT)
    T0 = t0 * np.ones(len(T))
    X = T - t0
    a = Sm - interp_s(X, T) / interp_c(X, T) - Sm * np.exp(-lambda_p * X + k * T) + np.exp(k * t0)
    avs.append((T, a))
    ax.scatter(T0, T, a, color='steelblue')
    # if len(a) >= 4:
    #     print(t0, a[-1], a[-2], a[-3], a[-4])
    #     d1 = (-(a[-1] - a[-2]) / RDT)
    #     d2 = (-(a[-2] - a[-3]) / RDT)
    #     d3 = (-(a[-3] - a[-4]) / RDT)
    #     print(
    #         d1 / a[-1],
    #         d2 / a[-2],
    #         d3 / a[-3]
    #     )
    #     dd1 = (d1 - d2) / RDT
    #     dd2 = (d2 - d3) / RDT
    #     print(dd1, dd2)
    #     print((dd1 - dd2) / RDT)

ax.set_xlabel('t0')
ax.set_ylabel('t')
plt.show()

tas = list()
ras = list()

for i in range(len(avs)):
    tas.append(tos[i])
    ras.append(avs[i][1][0])
    # plt.cla()
    # plt.clf()
    # plt.plot(avs[i][0], avs[i][1])
    # plt.title(f't0 = {tos[i]}')
    # plt.show()

plt.cla()
plt.clf()
plt.plot(tas, ras)
plt.title(f't0 = t')
plt.show()
