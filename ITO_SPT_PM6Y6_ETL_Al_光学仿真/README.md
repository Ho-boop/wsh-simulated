# ITO / SPT / PM6:Y6 / ETL / Al 光学仿真

## 器件结构

```
光 → 玻璃 / ITO (120 nm) / SPT 空穴传输层 (10 nm) / PM6:Y6 (100 nm) / ETL / Al (100 nm)
```

SPT = 磷钨酸钠（NaPW12，参考图 `参考图/4.7.tif` 中的 NaPW12）。

| 组别 | ETL |
|---|---|
| (a) w/o ETL | 无，Al 直接蒸镀在活性层上 |
| (b) SPT 5 nm | SPT 5 nm |
| (c) SPT 8 nm | SPT 8 nm |

> SPT 空穴传输层厚度没有给出，默认按 10 nm 计算，可以在 `tmm_sim.py` 顶部改 `D_HTL`。

## 方法

- 用传输矩阵法（TMM，Pettersson 1999）计算正入射下的光场 |E(z)|²/|E₀|²、各层吸收、反射率，再用 AM1.5G（ASTM G173）积分得到 Jsc,max（假设 IQE = 100%）和活性层产生率 G(z)。
- 平面多层膜在正入射下，FDTD 的 x–z 截面图沿 x 方向是均匀的，结果与 TMM 解析解严格一致。截面图按参考图的风格画（jet 色标、红色虚线标界面、标出活性层最大 |E|²）。
- 玻璃按非相干处理：考虑了空气/玻璃前表面反射以及玻璃内的多次反射。|E|² 以玻璃中的入射场归一化。
- 能量守恒自检：R + T + ΣA = 1，误差 < 2×10⁻⁵。

## 材料光学常数（`nk/`）

| 材料 | 来源 |
|---|---|
| Al | McPeak et al., *ACS Photonics* 2, 326 (2015)，蒸镀 Al（refractiveindex.info 数据库） |
| ITO | Minenkov et al., *ACS AMI* 16, 9517 (2024)，玻璃上 110 nm ITO（refractiveindex.info） |
| 玻璃 | Rubin, *Sol. Energy Mater.* 12, 275 (1985)，钠钙玻璃（refractiveindex.info） |
| PM6:Y6 | **模型数据**：KK 自洽的高斯振子模型，拟合到文献中典型的 PM6:Y6 消光系数（k_max ≈ 0.63 @ 800 nm） |
| SPT (NaPW12) | **模型数据**：可见/近红外透明（n ≈ 1.70），只在 < 350 nm 有紫外吸收 |

**注意：** PM6:Y6 和 SPT 没有公开的椭偏数据库，用的是 `build_model_nk.py` 生成的近似数据。如果有自己测的椭偏 n,k，直接替换 `nk/PM6Y6.csv`、`nk/SPT.csv`（格式 `wavelength_nm,n,k`），重跑 `tmm_sim.py` 就行。替换后 Jsc 的绝对值会变，三组之间的相对趋势通常比较稳。

## 结果（`results/`）

| 组别 | Jsc,max (mA/cm²) | PM6:Y6 平均吸收 (400–900 nm) | Al 寄生吸收损失 (mA/cm²) | ITO 寄生吸收损失 (mA/cm²) | 反射损失 (mA/cm²) | 活性层平均 \|E\|²（AM1.5G 加权） |
|---|---|---|---|---|---|---|
| w/o ETL | 25.74 | 0.750 | 3.59 | 1.69 | 7.06 | 1.043 |
| SPT 5 nm | 25.89 | 0.755 | 3.41 | 1.60 | 7.18 | 1.057 |
| SPT 8 nm | 25.90 | 0.756 | 3.31 | 1.56 | 7.31 | 1.061 |

主要结论：
1. 加入 SPT ETL 后，Jsc,max 从 25.74 增加到 25.89 / 25.90 mA/cm²（+0.6%）。光学增益很小，说明实验中 ETL 带来的 Jsc 提升主要应来自电学因素（界面偶极、功函数、电子提取/复合），而不是光学效应。
2. 增益来自 SPT 光学间隔层把活性层与 Al 隔开：Al 的寄生吸收从 3.59 降到 3.31 mA/cm²，650–800 nm 的吸收略有增加（700 nm 处从 0.761 增加到 0.790）。代价是 < 420 nm 的吸收略有下降（380 nm 处从 0.746 降到 0.703）。
3. 产生率分布（`generation_rate.png`）：加了 ETL 之后，峰值从距 HTL 约 29 nm 后移到约 37 nm，靠近 ETL/阴极一侧的产生率增加，分布更均匀。

文件清单：
- `E2_map_500nm.png`、`E2_map_600nm.png`、`E2_map_800nm.png`：三组器件的 |E|² 截面图（参考图风格）
- `E2_map_AM1.5G-weighted.png`：按 AM1.5G 光子通量加权（300–1000 nm）的平均 |E|²
- `E2_profile_z.png`：沿 z 方向的 |E|² 剖面对比
- `absorption_reflection.png`、`delta_absorption.png`：活性层吸收谱、器件反射谱，以及相对 w/o ETL 的吸收差
- `generation_rate.png`：活性层内的产生率 G(z)
- `nk_materials.png`：所用材料的 n,k
- `summary.csv`、`absorption_spectra.csv`、`E2_profile_*.csv`：原始数据（可导入 Origin 作图）

## 运行

```bash
pip install numpy scipy matplotlib
python3 build_model_nk.py   # 可选：重新生成 PM6:Y6 / SPT 的模型 n,k
python3 tmm_sim.py          # 输出到 results/
```

可以在 `tmm_sim.py` 顶部修改的参数：各层厚度、`GROUPS`（实验组）、`MAP_WAVELENGTHS`（截面图波长）、`GLASS_INCOHERENT`。
