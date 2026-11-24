# -*- coding: utf-8 -*-
"""
ITO / SPT / PM6:Y6 / ETL / Al 器件光学仿真（传输矩阵法 TMM）

器件：玻璃 / ITO(120) / SPT-HTL / PM6:Y6(100) / ETL / Al(100)
实验组（ETL 变量）：
    w/o ETL   : 无 ETL，Al 直接蒸镀在活性层上
    SPT 5 nm  : ETL = SPT 5 nm
    SPT 8 nm  : ETL = SPT 8 nm

对于平面多层膜、正入射，FDTD 得到的 |E|^2 截面图沿 x 方向是均匀的，
与 TMM 解析解严格一致；因此这里用 TMM 计算（无数值色散/网格误差，计算快）。
截面图按 FDTD 常见的 x–z 形式绘制，便于与参考图对比。

运行：python3 tmm_sim.py      结果输出到 results/
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
NK = os.path.join(HERE, "nk")
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

# ======================= 参数设置（按需修改） =======================
D_ITO = 120.0      # nm
D_HTL = 10.0       # nm，SPT 空穴传输层厚度（题目未给出，默认 10 nm，可修改）
D_ACT = 100.0      # nm，PM6:Y6
D_AL = 100.0       # nm
GROUPS = [          # (组名, ETL 厚度 nm)
    ("w/o ETL", 0.0),
    ("SPT 5 nm", 5.0),
    ("SPT 8 nm", 8.0),
]
MAP_WAVELENGTHS = [500, 600, 800]   # 输出 |E|^2 截面图的波长 (nm)
WL = np.arange(300, 1001, 1.0)      # 光谱计算范围 (nm)
DZ = 0.25                           # 纵向采样步长 (nm)
GLASS_INCOHERENT = True             # 计入空气/玻璃前表面反射及玻璃内非相干多次反射
# ===================================================================

Q = 1.602176634e-19
H = 6.62607015e-34
C = 2.99792458e8


def load_nk(name):
    d = np.loadtxt(os.path.join(NK, name + ".csv"), delimiter=",", comments="#")
    n = np.interp(WL, d[:, 0], d[:, 1])
    k = np.interp(WL, d[:, 0], d[:, 2])
    return n + 1j * k


MAT = {m: load_nk(f) for m, f in
       [("glass", "glass"), ("ITO", "ITO"), ("SPT", "SPT"), ("PM6:Y6", "PM6Y6"), ("Al", "Al")]}
AM = np.loadtxt(os.path.join(NK, "AM1.5G_ASTM_G173.csv"), delimiter=",", comments="#")
IRR = np.interp(WL, AM[:, 0], AM[:, 1])               # W m^-2 nm^-1
PHOTON = IRR * WL * 1e-9 / (H * C)                     # photons s^-1 m^-2 nm^-1


def stack(d_etl):
    """返回 [(材料, 厚度, 显示名)]，不含入射介质玻璃与出射介质空气"""
    s = [("ITO", D_ITO, "ITO"), ("SPT", D_HTL, "SPT"), ("PM6:Y6", D_ACT, "PM6:Y6")]
    if d_etl > 0:
        s.append(("SPT", d_etl, "SPT-ETL"))
    s.append(("Al", D_AL, "Al"))
    return s


def tmm_field(layers, iw):
    """
    Pettersson 传输矩阵法，正入射。
    入射介质：玻璃（半无限），出射介质：空气。
    返回 z(nm)、E(z)（以玻璃中入射电场 E0=1 归一化）、各层索引、R、T
    """
    lam = WL[iw]
    n0 = MAT["glass"][iw]
    ns = [n0] + [MAT[m][iw] for m, _, _ in layers] + [1.0 + 0j]
    ds = [0.0] + [d for _, d, _ in layers] + [0.0]
    N = len(ns)

    def I(a, b):
        r = (ns[a] - ns[b]) / (ns[a] + ns[b])
        t = 2 * ns[a] / (ns[a] + ns[b])
        return np.array([[1, r], [r, 1]], dtype=complex) / t

    def L(j):
        xi = 2 * np.pi * ns[j] / lam
        return np.array([[np.exp(-1j * xi * ds[j]), 0], [0, np.exp(1j * xi * ds[j])]], dtype=complex)

    S = I(0, 1)
    for j in range(1, N - 1):
        S = S @ L(j) @ I(j, j + 1)
    r = S[1, 0] / S[0, 0]
    t = 1 / S[0, 0]
    R = abs(r) ** 2
    T = abs(t) ** 2 * ns[-1].real / n0.real

    zs, Es, idx = [], [], []
    z0 = 0.0
    for j in range(1, N - 1):
        Sp = I(0, 1)
        for k in range(1, j):
            Sp = Sp @ L(k) @ I(k, k + 1)
        Spp = I(j, j + 1)
        for k in range(j + 1, N - 1):
            Spp = Spp @ L(k) @ I(k, k + 1)
        xi = 2 * np.pi * ns[j] / lam
        rpp = Spp[1, 0] / Spp[0, 0]
        tp = 1 / (Sp[0, 0] + Sp[0, 1] * rpp * np.exp(2j * xi * ds[j]))
        x = np.arange(0, ds[j], DZ) + DZ / 2
        E = tp * (np.exp(1j * xi * x) + rpp * np.exp(1j * xi * (2 * ds[j] - x)))
        zs.append(z0 + x)
        Es.append(E)
        idx.append(np.full(x.size, j - 1))
        z0 += ds[j]
    return np.concatenate(zs), np.concatenate(Es), np.concatenate(idx), R, T


def glass_factor(iw, R_stack):
    """空气→玻璃前表面透过，以及玻璃内非相干多次反射对进入膜系光强的修正"""
    if not GLASS_INCOHERENT:
        return 1.0
    ng = MAT["glass"][iw].real
    Rag = ((ng - 1) / (ng + 1)) ** 2
    return (1 - Rag) / (1 - Rag * R_stack)


def simulate(d_etl):
    layers = stack(d_etl)
    nL = len(layers)
    res = dict(layers=layers, A=np.zeros((nL, WL.size)), R=np.zeros(WL.size), T=np.zeros(WL.size))
    z, _, li, _, _ = tmm_field(layers, 0)
    res["z"], res["li"] = z, li
    res["E2"] = np.zeros((WL.size, z.size))
    res["G"] = np.zeros(z.size)       # 产生率 m^-3 s^-1
    for iw, lam in enumerate(WL):
        z, E, li, R, T = tmm_field(layers, iw)
        f = glass_factor(iw, R)
        E2 = abs(E) ** 2 * f
        n0 = MAT["glass"][iw].real
        ncs = np.array([MAT[m][iw] for m, _, _ in layers])[li]
        alpha = 4 * np.pi * ncs.imag / (lam * 1e-9)          # m^-1
        q = alpha * ncs.real / n0 * E2                        # 单位入射光强的吸收功率密度 (m^-1)
        for j in range(nL):
            res["A"][j, iw] = np.sum(q[li == j]) * DZ * 1e-9
        res["G"] += q * PHOTON[iw]                            # ×1 nm 波长步长
        res["E2"][iw] = E2
        res["T"][iw] = T * f
        # 器件总反射（从空气侧看）
        if GLASS_INCOHERENT:
            ng = n0
            Rag = ((ng - 1) / (ng + 1)) ** 2
            res["R"][iw] = Rag + (1 - Rag) ** 2 * R / (1 - Rag * R)
        else:
            res["R"][iw] = R
    act = [j for j, (_, _, nm) in enumerate(layers) if nm == "PM6:Y6"][0]
    res["act"] = act
    res["Jsc"] = Q * np.sum(res["A"][act] * PHOTON) * 1.0 / 10.0   # A/m^2 → mA/cm^2: ×1000/1e4
    return res


def boundaries(layers):
    b = np.cumsum([0] + [d for _, d, _ in layers])
    return b


def label_of(m, d, nm):
    if nm == "PM6:Y6":
        return f"PM6:Y6 ({d:g} nm)"
    if nm == "ITO":
        return f"ITO ({d:g}nm)"
    if nm == "Al":
        return f"Al ({d:g}nm)"
    if nm == "SPT-ETL":
        return f"SPT ({d:g} nm)"
    return f"SPT ({d:g} nm)"


def plot_maps(results, lam, vmax=None, tag=None):
    iw = int(np.argmin(abs(WL - lam)))
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw=dict(wspace=0.32))
    if vmax is None:
        vmax = max(r["E2"][iw].max() for r in results)
    letters = "abc"
    for ax, (name, _), r, L in zip(axes, GROUPS, results, letters):
        z = r["z"]
        x = np.array([-200, 200])
        img = r["E2"][iw][:, None]
        pc = ax.pcolormesh(x, np.r_[z - DZ / 2, z[-1] + DZ / 2], img, cmap="jet", vmin=0, vmax=vmax, shading="flat")
        b = boundaries(r["layers"])
        for zb in b[1:-1]:
            ax.axhline(zb, color="red", ls="--", lw=1.8)
        for (m, d, nm), z0 in zip(r["layers"], b[:-1]):
            zc = z0 + d / 2
            if d < 12:   # 薄层把标签放在层上方，避免遮挡
                zc = z0 + d + 10 if nm == "SPT-ETL" else z0 - 10
            ax.text(-190, zc, label_of(m, d, nm), color="white", fontsize=13, fontweight="bold", va="center",
                    path_effects=None)
        act = r["act"]
        sel = r["li"] == act
        ax.text(190, b[act] + 14, f"{r['E2'][iw][sel].max():.1f}", color="white", fontsize=15,
                fontweight="bold", ha="right", va="center")
        ax.set_xlim(-200, 200)
        ax.set_ylim(0, b[-1])
        ax.set_xticks([-200, -100, 0, 100, 200])
        ax.set_xlabel("x/nm", fontsize=13, fontweight="bold")
        ax.set_ylabel("z/nm", fontsize=13, fontweight="bold")
        ax.set_title(f"({L}) {name}", fontsize=13, fontweight="bold")
        ax.tick_params(labelsize=11)
    cb = fig.colorbar(pc, ax=axes, fraction=0.025, pad=0.02)
    cb.set_label(r"$|E|^2/|E_0|^2$", fontsize=12)
    title = f"λ = {lam:g} nm" if tag is None else tag
    fig.suptitle(f"|E|$^2$ distribution, {title}", fontsize=13, y=1.0)
    fn = os.path.join(OUT, f"E2_map_{lam:g}nm.png" if tag is None else f"E2_map_{tag.split()[0]}.png")
    fig.savefig(fn, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return fn


def main():
    results = [simulate(d) for _, d in GROUPS]
    colors = ["#1f77b4", "#d62728", "#2ca02c"]

    # ---------- 1. |E|^2 截面图（参考图风格） ----------
    for lam in MAP_WAVELENGTHS:
        plot_maps(results, lam)

    # AM1.5G 光子通量加权平均 |E|^2（300–1000 nm）
    w = PHOTON / PHOTON.sum()
    for r in results:
        r["E2avg"] = (w[:, None] * r["E2"]).sum(0)
    saved = []
    for r in results:
        saved.append(r["E2"])
        r["E2"] = np.tile(r["E2avg"], (WL.size, 1))
    plot_maps(results, WL[0], tag="AM1.5G-weighted (300–1000 nm)")
    for r, s in zip(results, saved):
        r["E2"] = s

    # ---------- 2. 沿 z 的 |E|^2 剖面 ----------
    fig, axes = plt.subplots(1, len(MAP_WAVELENGTHS), figsize=(5 * len(MAP_WAVELENGTHS), 4))
    for ax, lam in zip(axes, MAP_WAVELENGTHS):
        iw = int(np.argmin(abs(WL - lam)))
        for (name, _), r, c in zip(GROUPS, results, colors):
            ax.plot(r["z"], r["E2"][iw], color=c, label=name, lw=1.6)
            b = boundaries(r["layers"])
            ax.axvspan(b[r["act"]], b[r["act"] + 1], color="0.92", zorder=0)
        ax.set_title(f"λ = {lam} nm")
        ax.set_xlabel("z (nm)  [0 = glass/ITO]")
        ax.set_ylabel(r"$|E|^2/|E_0|^2$")
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "E2_profile_z.png"), dpi=200)
    plt.close(fig)

    # ---------- 3. 活性层吸收谱 & 反射谱 ----------
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for (name, _), r, c in zip(GROUPS, results, colors):
        ax[0].plot(WL, r["A"][r["act"]], color=c, label=f"{name}  (Jsc,max={r['Jsc']:.2f} mA/cm²)")
        ax[1].plot(WL, r["R"], color=c, label=name)
    ax[0].set_ylabel("Absorptance of PM6:Y6")
    ax[1].set_ylabel("Device reflectance R")
    for a in ax:
        a.set_xlabel("Wavelength (nm)")
        a.set_xlim(WL[0], WL[-1])
        a.legend(frameon=False, fontsize=9)
    ax[0].set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "absorption_reflection.png"), dpi=200)
    plt.close(fig)

    # 两组相对 w/o ETL 的吸收差
    fig, ax = plt.subplots(figsize=(6, 4))
    base = results[0]["A"][results[0]["act"]]
    for (name, _), r, c in list(zip(GROUPS, results, colors))[1:]:
        ax.plot(WL, r["A"][r["act"]] - base, color=c, label=f"{name} − w/o ETL")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("ΔAbsorptance (PM6:Y6)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "delta_absorption.png"), dpi=200)
    plt.close(fig)

    # ---------- 4. 活性层内产生率 G(z) ----------
    fig, ax = plt.subplots(figsize=(6, 4))
    for (name, _), r, c in zip(GROUPS, results, colors):
        sel = r["li"] == r["act"]
        b = boundaries(r["layers"])
        ax.plot(r["z"][sel] - b[r["act"]], r["G"][sel] * 1e-6, color=c, label=name)
    ax.set_xlabel("Position in PM6:Y6 from SPT-HTL side (nm)")
    ax.set_ylabel(r"Generation rate G (cm$^{-3}$ s$^{-1}$)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "generation_rate.png"), dpi=200)
    plt.close(fig)

    # ---------- 5. 汇总表 & 数据导出 ----------
    lines = []
    hdr = ("组别,ETL厚度(nm),Jsc_max(mA/cm2),PM6:Y6平均吸收(400-900nm),"
           "ITO寄生吸收Jloss(mA/cm2),SPT寄生吸收Jloss(mA/cm2),Al寄生吸收Jloss(mA/cm2),"
           "反射损失Jloss(mA/cm2),活性层平均|E|2(AM1.5G加权),能量守恒最大误差")
    lines.append(hdr)
    m49 = (WL >= 400) & (WL <= 900)
    for (name, d), r in zip(GROUPS, results):
        lay = r["layers"]

        def jl(arr):
            return Q * np.sum(arr * PHOTON) / 10.0
        parts = {}
        for j, (m, _, nm) in enumerate(lay):
            key = "SPT" if m == "SPT" else nm
            parts[key] = parts.get(key, 0) + jl(r["A"][j])
        sel = r["li"] == r["act"]
        # 能量守恒：R + T + ΣA_layers + 玻璃吸收(忽略) = 1
        cons = np.max(abs(r["R"] + r["T"] + r["A"].sum(0) - 1))
        lines.append(f"{name},{d:g},{r['Jsc']:.3f},{r['A'][r['act']][m49].mean():.4f},"
                     f"{parts['ITO']:.3f},{parts['SPT']:.3f},{parts['Al']:.3f},{jl(r['R']):.3f},"
                     f"{r['E2avg'][sel].mean():.3f},{cons:.1e}")
    with open(os.path.join(OUT, "summary.csv"), "w", encoding="utf-8-sig") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    head = "wavelength_nm," + ",".join(f"A_PM6Y6[{n}],R[{n}]" for n, _ in GROUPS)
    cols = [WL] + sum([[r["A"][r["act"]], r["R"]] for r in results], [])
    np.savetxt(os.path.join(OUT, "absorption_spectra.csv"), np.c_[tuple(cols)], delimiter=",",
               header=head, comments="", fmt="%.6g")
    for (name, _), r in zip(GROUPS, results):
        tag = name.replace("/", "").replace(" ", "_")
        iws = [int(np.argmin(abs(WL - l))) for l in MAP_WAVELENGTHS]
        layer_names = np.array([r["layers"][j][2] for j in r["li"]])
        data = np.c_[r["z"], *[r["E2"][i] for i in iws], r["E2avg"], r["G"] * 1e-6]
        with open(os.path.join(OUT, f"E2_profile_{tag}.csv"), "w") as fh:
            fh.write("z_nm,layer," + ",".join(f"E2_{l}nm" for l in MAP_WAVELENGTHS) +
                     ",E2_AM15G_weighted,G_cm-3s-1\n")
            for row, ln in zip(data, layer_names):
                fh.write(f"{row[0]:.3f},{ln}," + ",".join(f"{v:.6g}" for v in row[1:]) + "\n")

    # nk 一览图
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for m in ["ITO", "SPT", "PM6:Y6", "Al", "glass"]:
        ax[0].plot(WL, MAT[m].real, label=m)
        ax[1].plot(WL, MAT[m].imag, label=m)
    ax[0].set_ylabel("n"); ax[1].set_ylabel("k"); ax[1].set_yscale("symlog", linthresh=0.01)
    for a in ax:
        a.set_xlabel("Wavelength (nm)"); a.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "nk_materials.png"), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
