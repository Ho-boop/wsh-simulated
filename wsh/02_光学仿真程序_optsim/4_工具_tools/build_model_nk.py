# -*- coding: utf-8 -*-
"""
生成 PM6:Y6 与 SPT(磷钨酸钠, NaPW12) 的 n,k 数据（色散模型）。

公开数据库中没有这两种材料的椭偏数据，因此这里用
Kramers-Kronig 自洽的高斯振子模型（与椭偏软件中的 Gaussian oscillator 相同）构造 n,k：

    eps2(E) = A [exp(-((E-E0)/s)^2) - exp(-((E+E0)/s)^2)]
    eps1(E) = (2A/sqrt(pi)) [D((E+E0)/s) - D((E-E0)/s)]      D = Dawson 函数

PM6:Y6：振子能量取自 Y6 (≈1.53/1.68 eV) 与 PM6 (≈2.0/2.17 eV) 的 0-0/0-1 振动峰，
振幅用最小二乘拟合到文献典型的 PM6:Y6 (1:1.2) 消光系数曲线（k_max≈0.6 @ ~810 nm）。
SPT：Keggin 型多酸在可见/近红外透明，仅在紫外 (<350 nm) 有 O→W 电荷转移吸收，n≈1.65。

如有自己的椭偏实测数据，直接把 nk/PM6Y6.csv、nk/SPT.csv 替换为实测数据即可
（格式：wavelength_nm,n,k），主程序 tmm_sim.py 不需要改。
"""
import os
import numpy as np
from scipy.special import dawsn
from scipy.optimize import nnls

HERE = os.path.dirname(os.path.abspath(__file__))
HC = 1239.84193  # eV·nm


def gauss_eps(E, A, E0, s):
    eps2 = A * (np.exp(-((E - E0) / s) ** 2) - np.exp(-((E + E0) / s) ** 2))
    eps1 = 2 * A / np.sqrt(np.pi) * (dawsn((E + E0) / s) - dawsn((E - E0) / s))
    return eps1 + 1j * eps2


def nk_from_eps(eps):
    nc = np.sqrt(eps)
    return nc.real, np.abs(nc.imag)


wl = np.arange(300, 1201, 2.0)
E = HC / wl

# ---------------- PM6:Y6 ----------------
# 文献典型 k 曲线（PM6:Y6 共混膜，用于拟合振子振幅）
k_target_pts = np.array([
    [300, 0.32], [350, 0.30], [400, 0.21], [450, 0.17], [500, 0.22],
    [550, 0.36], [580, 0.44], [620, 0.44], [660, 0.34], [700, 0.36],
    [740, 0.48], [780, 0.58], [810, 0.62], [840, 0.50], [870, 0.26],
    [900, 0.08], [930, 0.02], [960, 0.004], [1000, 0.0], [1200, 0.0]])
osc = [  # (E0 eV, sigma eV)
    (1.43, 0.06), (1.52, 0.08), (1.64, 0.10), (1.80, 0.14),   # Y6
    (1.99, 0.09), (2.15, 0.12), (2.35, 0.20),        # PM6
    (2.75, 0.35), (3.5, 0.45), (4.3, 0.6),           # 高能吸收
]
eps_inf_pm = 2.9
k_t = np.interp(wl, k_target_pts[:, 0], k_target_pts[:, 1])
n_guess = 1.9
# eps2 ≈ 2 n k  → 对振幅线性 → 非负最小二乘
M = np.array([gauss_eps(E, 1.0, E0, s).imag for E0, s in osc]).T
A_fit = np.zeros(len(osc))
for _ in range(6):  # 迭代用自洽的 n 修正 eps2 目标
    eps2_t = 2 * n_guess * k_t
    A_fit, _ = nnls(M, eps2_t)
    eps = eps_inf_pm + sum(gauss_eps(E, A, E0, s) for A, (E0, s) in zip(A_fit, osc))
    n_guess, _ = nk_from_eps(eps)
n_pm, k_pm = nk_from_eps(eps)

# ---------------- SPT (NaPW12) ----------------
eps_spt = 2.55 + gauss_eps(E, 1.6, 4.9, 0.7) + gauss_eps(E, 0.35, 4.0, 0.35)
n_spt, k_spt = nk_from_eps(eps_spt)

os.makedirs(os.path.join(HERE, "..", "2_材料库_materials"), exist_ok=True)
np.savetxt(os.path.join(HERE, "..", "2_材料库_materials", "PM6Y6.csv"), np.c_[wl, n_pm, k_pm], delimiter=",", fmt="%.5g",
           header="wavelength_nm,n,k  | PM6:Y6 Gaussian-oscillator model (KK-consistent), fitted to typical literature k; "
                  "replace with measured ellipsometry data if available", comments="# ")
np.savetxt(os.path.join(HERE, "..", "2_材料库_materials", "SPT.csv"), np.c_[wl, n_spt, k_spt], delimiter=",", fmt="%.5g",
           header="wavelength_nm,n,k  | SPT (sodium phosphotungstate, NaPW12) Gaussian-oscillator model: transparent in VIS/NIR, "
                  "UV absorption only; replace with measured data if available", comments="# ")

if __name__ == "__main__":
    for w in (350, 450, 550, 620, 700, 800, 850, 900, 1000):
        i = np.argmin(abs(wl - w))
        print(f"{w:5d} nm  PM6:Y6 n={n_pm[i]:.3f} k={k_pm[i]:.3f} (target {k_t[i]:.2f})   SPT n={n_spt[i]:.3f} k={k_spt[i]:.4f}")
    print("PM6:Y6 振子振幅:", np.round(A_fit, 3))
