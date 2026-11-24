# -*- coding: utf-8 -*-
"""
optsim —— 多层薄膜器件（有机/钙钛矿太阳能电池等）光学仿真工具（传输矩阵法 TMM）

用法：
    python optsim.py configs/example_SPT_ETL.yaml
    python optsim.py 我的配置.yaml --quick        # 快速预览（波长步长 5 nm，不做扫描）

所有器件结构、厚度、材料、实验组、出图参数都写在 YAML 配置文件里，
模板见 configs/template.yaml（每一项都有中文注释）。

计算内容（正入射，光从衬底一侧入射）：
    |E(z)|^2 光场分布、各层吸收谱 A(λ)、器件反射 R / 透射 T、
    活性层 Jsc,max（AM1.5G，IQE 可设）、产生率 G(z)、各层寄生吸收损失、
    可选：单层厚度扫描（Jsc–厚度曲线）、双层厚度二维扫描（Jsc 等高图）
"""
import os
import re
import sys
import shutil
import argparse
import logging
import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

ROOT = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CONFIG = os.path.join(ROOT, "configs", "example_SPT_ETL.yaml")   # 不带参数运行时使用

Q_E = 1.602176634e-19
H_P = 6.62607015e-34
C_L = 2.99792458e8
trapz = getattr(np, "trapezoid", None) or np.trapz

# ------------------------------------------------------------------ 默认参数
DEFAULTS = {
    "output": None,                       # 输出目录，默认 results/<配置文件名>
    "wavelength": {"start": 300, "stop": 1000, "step": 1},
    "substrate": {"material": "glass", "incoherent": True},   # material: null → 光从空气直接入射到膜层
    "exit_medium": {"n": 1.0, "k": 0.0},
    "spectrum": "spectra/AM1.5G_ASTM_G173.csv",
    "iqe": 1.0,
    "dz": 0.25,
    "materials": {},
    "stack": [],
    "groups": None,
    "plot": {
        "map_wavelengths": [500, 600, 800],
        "map_weighted": True,
        "map_x_range": [-200, 200],
        "map_zlim": None,
        "map_vmax": None,
        "map_ncols": 3,
        "colormap": "jet",
        "interface_color": "red",
        "label_color": "white",
        "label_min_thickness": 15,
        "annotate_max": True,
        "font_size": 13,
        "dpi": 300,
        "avg_band": [400, 900],
        "formats": ["png"],
    },
    "sweep": None,
    "sweep2d": None,
}


def deep_update(base, new):
    out = dict(base)
    for k, v in (new or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_update(out[k], v)
        else:
            out[k] = v
    return out


def setup_font():
    """自动寻找中文字体，保证中文层名/组名能正常显示"""
    cands = ["Microsoft YaHei", "SimHei", "PingFang SC", "Heiti SC", "Noto Sans CJK SC",
             "Source Han Sans SC", "WenQuanYi Zen Hei", "WenQuanYi Micro Hei", "Arial Unicode MS"]
    have = {f.name for f in font_manager.fontManager.ttflist}
    for c in cands:
        if c in have:
            plt.rcParams["font.sans-serif"] = [c, "DejaVu Sans"]
            break
    plt.rcParams["axes.unicode_minus"] = False


def safe_name(s):
    return re.sub(r"[^\w\-.]+", "_", str(s)).strip("_") or "group"


# ------------------------------------------------------------------ 文件 & 材料
class Ctx:
    def __init__(self, cfg, cfg_dir):
        self.cfg = cfg
        self.cfg_dir = cfg_dir
        w = cfg["wavelength"]
        self.wl = np.arange(w["start"], w["stop"] + 1e-9, w["step"], dtype=float)
        self.wq = trapz_weights(self.wl)
        self._mat = {}
        spec = read_table(self.path(cfg["spectrum"]))
        irr = np.interp(self.wl, spec[:, 0], spec[:, 1], left=0, right=0)        # W m^-2 nm^-1
        self.photon = irr * self.wl * 1e-9 / (H_P * C_L)                         # s^-1 m^-2 nm^-1
        self.dz = float(cfg["dz"])

    def path(self, p):
        if os.path.isabs(p):
            return p
        for base in (self.cfg_dir, ROOT, os.getcwd()):
            q = os.path.join(base, p)
            if os.path.exists(q):
                return q
        raise FileNotFoundError(f"找不到文件：{p}（已在配置文件目录、{ROOT} 和当前目录中查找）")

    def nk(self, name):
        """返回材料在仿真波长上的复折射率 n + ik"""
        if isinstance(name, dict):                      # 层内直接写 {n:.., k:..}
            return self._spec_nk(name, str(name))
        if name in self._mat:
            return self._mat[name]
        spec = self.cfg["materials"].get(name)
        if spec is None:                                # 未定义 → 自动在 materials/ 下找同名 csv
            for ext in (".csv", ".txt", ".dat", ".nk"):
                p = os.path.join(ROOT, "materials", name + ext)
                if os.path.exists(p):
                    spec = p
                    break
        if spec is None:
            raise KeyError(f"材料 '{name}' 未定义：请在配置 materials: 中给出，或在 materials/ 下放 {name}.csv")
        self._mat[name] = self._spec_nk(spec, name)
        return self._mat[name]

    def _spec_nk(self, spec, name):
        wl = self.wl
        if isinstance(spec, str):
            spec = {"file": spec}
        if "file" in spec:
            d = read_table(self.path(spec["file"]))
            if d.shape[1] < 3:
                raise ValueError(f"材料 {name}: 文件需要 3 列 (波长, n, k)")
            x = d[:, 0] * (1000.0 if d[:, 0].max() < 20 else 1.0)    # µm → nm 自动识别
            o = np.argsort(x)
            x, n, k = x[o], d[o, 1], d[o, 2]
            if wl.min() < x.min() - 1 or wl.max() > x.max() + 1:
                print(f"  [警告] 材料 {name} 的数据范围 {x.min():.0f}–{x.max():.0f} nm 未覆盖仿真范围，超出部分取端点值")
            return np.interp(wl, x, n) + 1j * np.interp(wl, x, k)
        if "cauchy" in spec:                                          # n = A + B/λ² + C/λ⁴, λ 单位 µm
            c = list(spec["cauchy"]) + [0, 0, 0]
            lu = wl / 1000
            n = c[0] + c[1] / lu ** 2 + c[2] / lu ** 4
            return n + 1j * float(spec.get("k", 0.0)) * np.ones_like(wl)
        if "n" in spec:
            return float(spec["n"]) + 1j * float(spec.get("k", 0.0)) * np.ones_like(wl)
        raise ValueError(f"材料 {name} 的定义无法识别：{spec}")


def trapz_weights(x):
    if x.size == 1:
        return np.ones(1)
    w = np.zeros_like(x)
    dx = np.diff(x)
    w[:-1] += dx / 2
    w[1:] += dx / 2
    return w


def read_table(p):
    """读取任意分隔符（逗号/制表符/空格/分号）的数值表，自动跳过表头和注释"""
    rows = []
    with open(p, "r", encoding="utf-8-sig", errors="ignore") as fh:
        for line in fh:
            parts = [s for s in re.split(r"[,\t; ]+", line.strip()) if s]
            try:
                rows.append([float(s) for s in parts])
            except ValueError:
                continue
    rows = [r for r in rows if len(r) >= 2]
    if not rows:
        raise ValueError(f"文件 {p} 中没有读到数值数据")
    m = min(len(r) for r in rows)
    return np.array([r[:m] for r in rows])


# ------------------------------------------------------------------ 器件结构
def build_stack(cfg, group):
    """根据基础 stack 与实验组的修改，生成该组的层列表（厚度为 0 的层自动去掉）"""
    if group.get("stack"):
        base = group["stack"]
    else:
        base = cfg["stack"]
    layers = []
    names = [L.get("name", L["material"] if isinstance(L.get("material"), str) else f"L{i}")
             for i, L in enumerate(base)]
    sets = group.get("set") or {}
    for key in sets:
        if key not in names:
            raise KeyError(f"实验组 '{group['name']}' 的 set 中引用了不存在的层名 '{key}'，可用层名：{names}")
    for L, nm in zip(base, names):
        L = dict(L)
        L["name"] = nm
        if nm in sets:
            v = sets[nm]
            L.update(v if isinstance(v, dict) else {"thickness": v})
        if float(L.get("thickness", 0)) <= 0:
            continue
        L["thickness"] = float(L["thickness"])
        L.setdefault("active", False)
        if "label" not in L:
            mname = L["material"] if isinstance(L["material"], str) else nm
            L["label"] = f"{mname} ({L['thickness']:g} nm)"
        else:
            L["label"] = str(L["label"]).replace("{d}", f"{L['thickness']:g}")
        layers.append(L)
    if not any(L["active"] for L in layers):
        raise ValueError(f"实验组 '{group['name']}' 没有活性层：请给活性层加上 active: true")
    return layers


# ------------------------------------------------------------------ TMM 核心（对波长向量化）
def _I(na, nb):
    r = (na - nb) / (na + nb)
    t = 2 * na / (na + nb)
    M = np.empty(na.shape + (2, 2), complex)
    M[:, 0, 0] = M[:, 1, 1] = 1 / t
    M[:, 0, 1] = M[:, 1, 0] = r / t
    return M


def _L(n, d, wl):
    xi = 2 * np.pi * n / wl
    M = np.zeros(n.shape + (2, 2), complex)
    M[:, 0, 0] = np.exp(-1j * xi * d)
    M[:, 1, 1] = np.exp(1j * xi * d)
    return M


def tmm(nks, ds, n_in, n_out, wl, dz):
    """
    Pettersson 传输矩阵法。nks: 各层复折射率 [(nw,)...]，ds: 厚度 (nm)。
    返回 z 采样点、各点所属层、各点积分权重、E(nw,nz)（入射介质中 E0+=1 归一化）、r、t
    """
    N = len(nks)
    n_all = [n_in] + list(nks) + [n_out]
    P = [None] * (N + 2)
    P[1] = _I(n_all[0], n_all[1])
    for j in range(2, N + 1):
        P[j] = P[j - 1] @ _L(n_all[j - 1], ds[j - 2], wl) @ _I(n_all[j - 1], n_all[j])
    S = P[N] @ _L(n_all[N], ds[N - 1], wl) @ _I(n_all[N], n_all[N + 1])
    Qm = [None] * (N + 2)
    Qm[N] = _I(n_all[N], n_all[N + 1])
    for j in range(N - 1, 0, -1):
        Qm[j] = _I(n_all[j], n_all[j + 1]) @ _L(n_all[j + 1], ds[j], wl) @ Qm[j + 1]
    r = S[:, 1, 0] / S[:, 0, 0]
    t = 1 / S[:, 0, 0]
    zs, Es, li, wz = [], [], [], []
    z0 = 0.0
    for j in range(1, N + 1):
        d = ds[j - 1]
        nx = max(2, int(np.ceil(d / dz)))
        h = d / nx
        x = (np.arange(nx) + 0.5) * h
        xi = (2 * np.pi * n_all[j] / wl)[:, None]
        rpp = (Qm[j][:, 1, 0] / Qm[j][:, 0, 0])[:, None]
        tp = 1 / (P[j][:, 0, 0][:, None] + P[j][:, 0, 1][:, None] * rpp * np.exp(2j * xi * d))
        Es.append(tp * (np.exp(1j * xi * x) + rpp * np.exp(1j * xi * (2 * d - x))))
        zs.append(z0 + x)
        li.append(np.full(nx, j - 1))
        wz.append(np.full(nx, h))
        z0 += d
    return np.concatenate(zs), np.concatenate(li), np.concatenate(wz), np.concatenate(Es, axis=1), r, t


def simulate(ctx, layers, keep_field=True):
    cfg = ctx.cfg
    wl = ctx.wl
    sub = cfg.get("substrate") or {}
    if sub.get("material"):
        n_in = ctx.nk(sub["material"])
        incoh = bool(sub.get("incoherent", True))
    else:
        n_in = np.ones_like(wl, dtype=complex)
        incoh = False
    ex = cfg.get("exit_medium") or {"n": 1.0}
    n_out = ctx.nk(ex)
    nks = [ctx.nk(L["material"]) for L in layers]
    ds = [L["thickness"] for L in layers]
    z, li, wz, E, r, t = tmm(nks, ds, n_in, n_out, wl, ctx.dz)
    R = np.abs(r) ** 2
    T = np.abs(t) ** 2 * n_out.real / n_in.real
    if incoh:   # 厚衬底：空气/衬底前表面反射 + 衬底内非相干多次反射
        ng = n_in.real
        Rag = ((ng - 1) / (ng + 1)) ** 2
        f = (1 - Rag) / (1 - Rag * R)
        Rdev = Rag + (1 - Rag) ** 2 * R / (1 - Rag * R)
    else:
        f = np.ones_like(R)
        Rdev = R
    E2 = np.abs(E) ** 2 * f[:, None]
    ncz = np.array(nks)[li].T                                  # (nw, nz)
    alpha = 4 * np.pi * ncz.imag / (wl[:, None] * 1e-9)        # m^-1
    q = alpha * ncz.real / n_in.real[:, None] * E2             # 单位入射光强下的吸收功率密度 (m^-1)
    nL = len(layers)
    A = np.zeros((nL, wl.size))
    for j in range(nL):
        s = li == j
        A[j] = (q[:, s] * wz[s]).sum(1) * 1e-9
    G = ((q * ctx.photon[:, None]) * ctx.wq[:, None]).sum(0)    # m^-3 s^-1
    act = [j for j, L in enumerate(layers) if L["active"]]
    res = dict(layers=layers, z=z, li=li, wz=wz, A=A, R=Rdev, T=T * f, G=G, act=act,
               A_act=A[act].sum(0))
    res["cons_err"] = float(np.max(np.abs(Rdev + T * f + A.sum(0) - 1)))
    res["Jsc"] = jcur(ctx, res["A_act"]) * float(cfg["iqe"])
    if keep_field:
        res["E2"] = E2
        w = ctx.photon * ctx.wq
        res["E2w"] = (w[:, None] * E2).sum(0) / w.sum()
    return res


def jcur(ctx, A):
    """吸收谱 → 等效电流密度 mA/cm²"""
    return Q_E * np.sum(A * ctx.photon * ctx.wq) / 10.0


def bounds(layers):
    return np.cumsum([0] + [L["thickness"] for L in layers])


# ------------------------------------------------------------------ 出图
def savefig(fig, out, name, cfg):
    for fmt in cfg["plot"]["formats"]:
        fig.savefig(os.path.join(out, f"{name}.{fmt}"), dpi=cfg["plot"]["dpi"], bbox_inches="tight")
    plt.close(fig)


def plot_maps(ctx, groups, results, out, field_key, iw=None, title="", fname="map"):
    cfg = ctx.cfg
    pc_ = cfg["plot"]
    fs = pc_["font_size"]
    n = len(results)
    ncols = min(n, int(pc_["map_ncols"]))
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.6 * nrows), squeeze=False,
                             gridspec_kw=dict(wspace=0.32, hspace=0.35))
    data = [(r[field_key][iw] if iw is not None else r[field_key]) for r in results]
    vmax = pc_["map_vmax"] or max(d.max() for d in data)
    x0, x1 = pc_["map_x_range"]
    pcm = None
    for k, (ax, g, r, d) in enumerate(zip(axes.flat, groups, results, data)):
        z, wz = r["z"], r["wz"]
        edges = np.r_[z - wz / 2, z[-1] + wz[-1] / 2]
        pcm = ax.pcolormesh([x0, x1], edges, d[:, None], cmap=pc_["colormap"], vmin=0, vmax=vmax,
                            shading="flat", rasterized=True)
        b = bounds(r["layers"])
        for zb in b[1:-1]:
            ax.axhline(zb, color=pc_["interface_color"], ls="--", lw=1.8)
        tot = b[-1]
        for j, L in enumerate(r["layers"]):
            if L.get("show_label", True) is False:
                continue
            zc = (b[j] + b[j + 1]) / 2
            if L["thickness"] < pc_["label_min_thickness"]:      # 薄层：标签放到层外侧
                zc = b[j + 1] + 9 if zc > tot / 2 else b[j] - 9
            zc += float(L.get("label_offset", 0))
            ax.text(x0 + 0.025 * (x1 - x0), zc, L["label"], color=pc_["label_color"], fontsize=fs,
                    fontweight="bold", va="center", clip_on=True)
        if pc_["annotate_max"]:
            for j in r["act"]:
                s = r["li"] == j
                ax.text(x1 - 0.025 * (x1 - x0), b[j] + 0.15 * r["layers"][j]["thickness"],
                        f"{d[s].max():.2f}", color=pc_["label_color"], fontsize=fs + 2,
                        fontweight="bold", ha="right", va="center")
        ax.set_xlim(x0, x1)
        zl = pc_["map_zlim"] or [0, tot]
        ax.set_ylim(zl[0], zl[1] if zl[1] is not None else tot)
        ax.set_xlabel("x/nm", fontsize=fs, fontweight="bold")
        ax.set_ylabel("z/nm", fontsize=fs, fontweight="bold")
        ax.set_title(f"({chr(97 + k)}) {g['name']}", fontsize=fs, fontweight="bold")
        ax.tick_params(labelsize=fs - 2)
    for ax in list(axes.flat)[n:]:
        ax.set_visible(False)
    cb = fig.colorbar(pcm, ax=axes, fraction=0.025, pad=0.02)
    cb.set_label(r"$|E|^2/|E_0|^2$", fontsize=fs - 1)
    fig.suptitle(title, fontsize=fs, y=0.88 + 0.65 / fig.get_figheight())
    savefig(fig, out, fname, cfg)


def make_plots(ctx, groups, results, out):
    cfg = ctx.cfg
    pc_ = cfg["plot"]
    wl = ctx.wl
    cols = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    names = [g["name"] for g in groups]
    map_wls = [w for w in pc_["map_wavelengths"] if wl.min() <= w <= wl.max()]
    for w in pc_["map_wavelengths"]:
        if w not in map_wls:
            print(f"  [提示] 截面图波长 {w} nm 超出计算范围，已跳过")

    # 1. |E|² 截面图
    for w in map_wls:
        iw = int(np.argmin(np.abs(wl - w)))
        plot_maps(ctx, groups, results, out, "E2", iw, f"|E|$^2$ distribution, λ = {wl[iw]:g} nm",
                  f"E2_map_{wl[iw]:g}nm")
    if pc_["map_weighted"]:
        plot_maps(ctx, groups, results, out, "E2w", None,
                  f"|E|$^2$ distribution, AM1.5G photon-flux weighted ({wl.min():g}–{wl.max():g} nm)",
                  "E2_map_AM1.5G_weighted")

    # 2. |E|² 沿 z 剖面
    panels = [(w, int(np.argmin(np.abs(wl - w)))) for w in map_wls] + ([("AM1.5G", None)] if pc_["map_weighted"] else [])
    if panels:
        fig, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 4), squeeze=False)
        for ax, (w, iw) in zip(axes[0], panels):
            for k, r in enumerate(results):
                y = r["E2"][iw] if iw is not None else r["E2w"]
                ax.plot(r["z"], y, color=cols[k % 10], label=names[k], lw=1.6)
                b = bounds(r["layers"])
                for zb in b[1:-1]:
                    ax.axvline(zb, color=cols[k % 10], lw=0.6, ls=":")
            ax.set_title(f"λ = {wl[iw]:g} nm" if iw is not None else "AM1.5G weighted")
            ax.set_xlabel("z (nm)   [0 = substrate/first layer]")
            ax.set_ylabel(r"$|E|^2/|E_0|^2$")
            ax.legend(frameon=False, fontsize=9)
        fig.tight_layout()
        savefig(fig, out, "E2_profile_z", cfg)

    # 3. 活性层吸收 & 反射
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for k, r in enumerate(results):
        ax[0].plot(wl, r["A_act"], color=cols[k % 10], label=f"{names[k]}  (Jsc={r['Jsc']:.2f} mA/cm²)")
        ax[1].plot(wl, r["R"], color=cols[k % 10], label=names[k])
    ax[0].set_ylabel("Absorptance of active layer")
    ax[0].set_ylim(0, 1)
    ax[1].set_ylabel("Device reflectance R")
    for a in ax:
        a.set_xlabel("Wavelength (nm)")
        a.set_xlim(wl.min(), wl.max())
        a.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    savefig(fig, out, "absorption_reflection", cfg)

    if len(results) > 1:
        fig, ax = plt.subplots(figsize=(6, 4))
        for k, r in enumerate(results[1:], 1):
            ax.plot(wl, r["A_act"] - results[0]["A_act"], color=cols[k % 10], label=f"{names[k]} − {names[0]}")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("ΔAbsorptance (active layer)")
        ax.legend(frameon=False)
        fig.tight_layout()
        savefig(fig, out, "delta_absorption", cfg)

    # 4. 产生率
    fig, ax = plt.subplots(figsize=(6, 4))
    for k, r in enumerate(results):
        s = np.isin(r["li"], r["act"])
        z0 = bounds(r["layers"])[r["act"][0]]
        ax.plot(r["z"][s] - z0, r["G"][s] * 1e-6, color=cols[k % 10], label=names[k])
    ax.set_xlabel("Position in active layer from substrate side (nm)")
    ax.set_ylabel(r"Generation rate G (cm$^{-3}$ s$^{-1}$)")
    ax.legend(frameon=False)
    fig.tight_layout()
    savefig(fig, out, "generation_rate", cfg)

    # 5. 光子损失分解（堆叠柱状图）
    cats = ["Active (Jsc)"]
    for r in results:
        for L in r["layers"]:
            if not L["active"] and f"{L['name']} abs." not in cats:
                cats.append(f"{L['name']} abs.")
    cats += ["Reflection", "Transmission"]
    vals = np.zeros((len(results), len(cats)))
    for i, r in enumerate(results):
        vals[i, 0] = jcur(ctx, r["A_act"])
        for j, L in enumerate(r["layers"]):
            if not L["active"]:
                vals[i, cats.index(f"{L['name']} abs.")] += jcur(ctx, r["A"][j])
        vals[i, -2] = jcur(ctx, r["R"])
        vals[i, -1] = jcur(ctx, r["T"])
    fig, ax = plt.subplots(figsize=(1.6 * len(results) + 4, 4.5))
    bot = np.zeros(len(results))
    cmap = plt.get_cmap("tab20")
    for c in range(len(cats)):
        ax.bar(names, vals[:, c], bottom=bot, color=cmap(c % 20), label=cats[c], width=0.6)
        for i in range(len(results)):
            if vals[i, c] > 0.6:
                ax.text(i, bot[i] + vals[i, c] / 2, f"{vals[i, c]:.2f}", ha="center", va="center", fontsize=8)
        bot += vals[:, c]
    ax.set_ylabel(f"Photon current (mA/cm²), {wl.min():g}–{wl.max():g} nm")
    ax.set_ylim(0, bot.max() * 1.05)
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    savefig(fig, out, "loss_breakdown", cfg)

    # 6. 材料 n,k
    mats = []
    for r in results:
        for L in r["layers"]:
            m = L["material"]
            if isinstance(m, str) and m not in mats:
                mats.append(m)
    sub = (cfg.get("substrate") or {}).get("material")
    if sub and sub not in mats:
        mats.append(sub)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for m in mats:
        v = ctx.nk(m)
        ax[0].plot(wl, v.real, label=m)
        ax[1].plot(wl, v.imag, label=m)
    ax[0].set_ylabel("n")
    ax[1].set_ylabel("k")
    ax[1].set_yscale("symlog", linthresh=0.01)
    for a in ax:
        a.set_xlabel("Wavelength (nm)")
        a.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    savefig(fig, out, "nk_materials", cfg)
    return cats, vals


# ------------------------------------------------------------------ 数据导出
def write_csv(path, header, rows):
    with open(path, "w", encoding="utf-8-sig") as fh:
        fh.write(",".join(header) + "\n")
        for row in rows:
            fh.write(",".join(f"{v:.6g}" if isinstance(v, (float, np.floating)) else str(v) for v in row) + "\n")


def export(ctx, groups, results, out, cats, vals):
    cfg = ctx.cfg
    wl = ctx.wl
    band = cfg["plot"]["avg_band"]
    mb = (wl >= band[0]) & (wl <= band[1])
    map_wls = [w for w in cfg["plot"]["map_wavelengths"] if wl.min() <= w <= wl.max()]
    iws = [int(np.argmin(np.abs(wl - w))) for w in map_wls]

    hdr = ["group", "structure", "Jsc_max(mA/cm2)", f"active_avg_A({band[0]}-{band[1]}nm)",
           "active_avg_E2(AM1.5G_weighted)"] + [f"{c}(mA/cm2)" for c in cats] + ["energy_conservation_err"]
    rows = []
    for g, r, v in zip(groups, results, vals):
        s = np.isin(r["li"], r["act"])
        struct = " / ".join(f"{L['material'] if isinstance(L['material'], str) else L['name']} {L['thickness']:g}"
                            for L in r["layers"])
        rows.append([g["name"], struct, r["Jsc"], r["A_act"][mb].mean(), np.average(r["E2w"][s], weights=r["wz"][s]),
                     *v, r["cons_err"]])
    write_csv(os.path.join(out, "summary.csv"), hdr, rows)

    write_csv(os.path.join(out, "absorption_all_groups.csv"),
              ["wavelength_nm"] + [f"A_active[{g['name']}]" for g in groups] + [f"R[{g['name']}]" for g in groups],
              np.c_[wl, *[r["A_act"] for r in results], *[r["R"] for r in results]].tolist())

    for g, r in zip(groups, results):
        tag = safe_name(g["name"])
        write_csv(os.path.join(out, f"spectra_{tag}.csv"),
                  ["wavelength_nm"] + [f"A[{L['name']}]" for L in r["layers"]] + ["A_active_total", "R", "T"],
                  np.c_[wl, *r["A"], r["A_act"], r["R"], r["T"]].tolist())
        lname = [r["layers"][j]["name"] for j in r["li"]]
        rows = [[r["z"][i], lname[i], *[r["E2"][iw][i] for iw in iws], r["E2w"][i], r["G"][i] * 1e-6]
                for i in range(r["z"].size)]
        write_csv(os.path.join(out, f"field_{tag}.csv"),
                  ["z_nm", "layer"] + [f"E2_{wl[iw]:g}nm" for iw in iws] + ["E2_AM1.5G_weighted", "G_cm-3s-1"], rows)
    print("\n汇总结果 (summary.csv)：")
    print(f"  {'组别':<16}{'Jsc,max (mA/cm²)':>18}{'能量守恒误差':>14}")
    for g, r in zip(groups, results):
        print(f"  {g['name']:<16}{r['Jsc']:>18.3f}{r['cons_err']:>14.1e}")


# ------------------------------------------------------------------ 厚度扫描
def sweep_values(spec):
    if "values" in spec:
        return [float(v) for v in spec["values"]]
    a, b, s = spec["range"]
    return list(np.arange(a, b + 1e-9, s))


def run_sweep(ctx, groups, out):
    sw = ctx.cfg["sweep"]
    vals = sweep_values(sw)
    target = sw["layer"]
    print(f"\n厚度扫描：层 '{target}'，{len(vals)} 个厚度 × {len(groups)} 组 ...")
    table = np.zeros((len(vals), len(groups)))
    for k, g in enumerate(groups):
        for i, v in enumerate(vals):
            g2 = dict(g)
            g2["set"] = dict(g.get("set") or {})
            old = g2["set"].get(target)
            g2["set"][target] = dict(old, thickness=v) if isinstance(old, dict) else v
            table[i, k] = simulate(ctx, build_stack(ctx.cfg, g2), keep_field=False)["Jsc"]
    tag = safe_name(target)
    write_csv(os.path.join(out, f"sweep_{tag}.csv"), [f"{target}_thickness_nm"] + [f"Jsc[{g['name']}]" for g in groups],
              np.c_[vals, table].tolist())
    fig, ax = plt.subplots(figsize=(6, 4))
    for k, g in enumerate(groups):
        ax.plot(vals, table[:, k], "o-", ms=4, label=g["name"])
        i = int(np.argmax(table[:, k]))
        ax.plot(vals[i], table[i, k], "*", ms=12, color=ax.lines[-1].get_color())
    ax.set_xlabel(f"{target} thickness (nm)")
    ax.set_ylabel("Jsc,max (mA/cm²)")
    ax.legend(frameon=False)
    fig.tight_layout()
    savefig(fig, out, f"sweep_{tag}", ctx.cfg)


def run_sweep2d(ctx, groups, out):
    sw = ctx.cfg["sweep2d"]
    gname = sw.get("group", groups[0]["name"])
    g = next((x for x in groups if x["name"] == gname), None)
    if g is None:
        raise KeyError(f"sweep2d.group '{gname}' 不在实验组中")
    lx, ly = sw["x"]["layer"], sw["y"]["layer"]
    xv, yv = sweep_values(sw["x"]), sweep_values(sw["y"])
    print(f"\n二维扫描：{lx} × {ly}（基于组 '{gname}'），{len(xv) * len(yv)} 个结构 ...")
    J = np.zeros((len(yv), len(xv)))
    for i, y in enumerate(yv):
        for k, x in enumerate(xv):
            g2 = dict(g)
            st = dict(g.get("set") or {})
            for lay, v in ((lx, x), (ly, y)):
                old = st.get(lay)
                st[lay] = dict(old, thickness=v) if isinstance(old, dict) else v
            g2["set"] = st
            J[i, k] = simulate(ctx, build_stack(ctx.cfg, g2), keep_field=False)["Jsc"]
    tag = f"{safe_name(lx)}_vs_{safe_name(ly)}"
    write_csv(os.path.join(out, f"sweep2d_{tag}.csv"), [f"{ly}\\{lx}"] + [f"{x:g}" for x in xv],
              [[y, *J[i]] for i, y in enumerate(yv)])
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    cs = ax.contourf(xv, yv, J, levels=20, cmap="viridis")
    ax.contour(xv, yv, J, levels=10, colors="w", linewidths=0.4)
    i, k = np.unravel_index(np.argmax(J), J.shape)
    ax.plot(xv[k], yv[i], "r*", ms=14)
    ax.set_title(f"max Jsc = {J[i, k]:.2f} mA/cm² @ {lx}={xv[k]:g}, {ly}={yv[i]:g} nm", fontsize=10)
    ax.set_xlabel(f"{lx} thickness (nm)")
    ax.set_ylabel(f"{ly} thickness (nm)")
    fig.colorbar(cs, label="Jsc,max (mA/cm²)")
    fig.tight_layout()
    savefig(fig, out, f"sweep2d_{tag}", ctx.cfg)


# ------------------------------------------------------------------ 主程序
def main(argv=None):
    ap = argparse.ArgumentParser(description="多层薄膜器件光学仿真（TMM）")
    ap.add_argument("config", nargs="?", default=DEFAULT_CONFIG, help="YAML 配置文件")
    ap.add_argument("--quick", action="store_true", help="快速预览：波长步长 5 nm，跳过厚度扫描")
    ap.add_argument("-o", "--output", help="输出目录（覆盖配置中的 output）")
    a = ap.parse_args(argv)

    with open(a.config, "r", encoding="utf-8") as fh:
        user = yaml.safe_load(fh) or {}
    cfg = deep_update(DEFAULTS, user)
    if a.quick:
        cfg["wavelength"]["step"] = max(5, cfg["wavelength"]["step"])
        cfg["sweep"] = cfg["sweep2d"] = None
    cfg_dir = os.path.dirname(os.path.abspath(a.config))
    if not cfg["stack"]:
        raise ValueError("配置中没有 stack（器件结构）")
    groups = cfg["groups"] or [{"name": "device"}]
    for i, g in enumerate(groups):
        g.setdefault("name", f"group{i + 1}")

    out = a.output or cfg["output"] or os.path.join("results", os.path.splitext(os.path.basename(a.config))[0])
    if not os.path.isabs(out):
        out = os.path.join(ROOT, out)
    os.makedirs(out, exist_ok=True)
    shutil.copy(a.config, os.path.join(out, "used_config.yaml"))

    setup_font()
    plt.rcParams["font.size"] = 11
    ctx = Ctx(cfg, cfg_dir)
    print(f"配置：{a.config}\n输出：{out}\n波长：{ctx.wl.min():g}–{ctx.wl.max():g} nm，步长 {cfg['wavelength']['step']} nm")
    results = []
    for g in groups:
        layers = build_stack(cfg, g)
        print(f"  [{g['name']}]  " + " / ".join(L["label"] for L in layers))
        results.append(simulate(ctx, layers))
    cats, vals = make_plots(ctx, groups, results, out)
    export(ctx, groups, results, out, cats, vals)
    if cfg.get("sweep"):
        run_sweep(ctx, groups, out)
    if cfg.get("sweep2d"):
        run_sweep2d(ctx, groups, out)
    print(f"\n完成，图和数据已保存到：{out}")


if __name__ == "__main__":
    main()
