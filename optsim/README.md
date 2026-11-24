# optsim —— 多层薄膜器件光学仿真工具（TMM）

适用于有机/钙钛矿太阳能电池等平面多层器件。器件结构、厚度、材料、实验组和出图参数都写在一个 YAML 配置文件里，运行一条命令就能得到全部图和数据。

## 快速开始

```bash
pip install -r requirements.txt
python optsim.py configs/example_SPT_ETL.yaml
```

新建自己的器件时，复制 `configs/template.yaml`（每项都有中文注释），改好后运行：

```bash
python optsim.py configs/我的器件.yaml
python optsim.py configs/我的器件.yaml --quick     # 快速预览：波长步长 5 nm，跳过厚度扫描
python optsim.py configs/我的器件.yaml -o D:/out   # 指定输出目录
```

在 Spyder/PyCharm 里直接点运行也可以：修改 `optsim.py` 顶部的 `DEFAULT_CONFIG`。

## 目录

```
optsim/
├─ optsim.py              主程序（计算 + 出图 + 导出数据）
├─ configs/
│   ├─ template.yaml      全注释模板（所有可设参数）
│   └─ example_SPT_ETL.yaml  例子：ITO/SPT/PM6:Y6/ETL/Al，ETL = 无 / SPT 5 nm / SPT 8 nm
├─ materials/             n,k 材料库（每种材料一个 csv：波长, n, k）
├─ spectra/               AM1.5G 光谱（ASTM G173）
├─ tools/build_model_nk.py   PM6:Y6、SPT 模型 n,k 的生成脚本
└─ results/<配置名>/       输出
```

## 最简配置

```yaml
stack:                                   # 从衬底（入光侧）往上写，单位 nm
  - {name: ITO,    material: ITO,   thickness: 120}
  - {name: HTL,    material: SPT,   thickness: 10}
  - {name: Active, material: PM6Y6, thickness: 100, active: true}
  - {name: ETL,    material: SPT,   thickness: 0}
  - {name: Al,     material: Al,    thickness: 100}
groups:                                  # 实验组：在 stack 基础上改某几层
  - {name: "w/o ETL",  set: {ETL: 0}}    # 厚度 0 = 去掉这一层
  - {name: "SPT 5 nm", set: {ETL: 5}}
  - {name: "PDINN 8 nm", set: {ETL: {material: PDINN, thickness: 8}}}   # 换材料
```

## 全部参数及默认值

| 参数 | 默认值 | 说明 |
|---|---|---|
| `output` | `results/<配置名>` | 输出目录 |
| `wavelength` | `{start: 300, stop: 1000, step: 1}` | 计算波长 (nm) |
| `substrate.material` | `glass` | 入光侧衬底；`null` 表示光直接从空气入射到第一层 |
| `substrate.incoherent` | `true` | 厚衬底按非相干处理，计入空气/玻璃前表面反射和衬底内多次反射 |
| `exit_medium` | `{n: 1, k: 0}` | 最后一层后面的介质 |
| `spectrum` | `spectra/AM1.5G_ASTM_G173.csv` | 太阳光谱 |
| `iqe` | `1.0` | 内量子效率，Jsc = IQE × 吸收的光子电流 |
| `dz` | `0.25` | 纵向采样步长 (nm) |
| `materials` | 空 | 额外材料定义（见下） |
| `stack` | 必填 | 层列表：`name` / `material` / `thickness` / `active` / `label` / `show_label` / `label_offset` |
| `groups` | 单组 | `name` + `set`（改层），或 `name` + `stack`（完全自定义） |
| `plot.map_wavelengths` | `[500, 600, 800]` | 截面图波长 |
| `plot.map_weighted` | `true` | 额外输出 AM1.5G 加权的平均 \|E\|² 截面图 |
| `plot.map_x_range` / `map_zlim` / `map_vmax` | `[-200,200]` / 整个器件 / 自动 | 截面图坐标范围和色标上限 |
| `plot.colormap` / `interface_color` / `label_color` | `jet` / `red` / `white` | 配色 |
| `plot.label_min_thickness` | `15` | 比这个薄的层，标签放到层外侧 |
| `plot.annotate_max` | `true` | 在活性层标出 \|E\|² 最大值 |
| `plot.font_size` / `dpi` / `formats` | `13` / `300` / `[png]` | 可输出 `pdf`、`svg`、`tif` |
| `plot.avg_band` | `[400, 900]` | summary 中“平均吸收”的统计波段 |
| `sweep` | 无 | 一维厚度扫描：`{layer: Active, range: [60, 200, 10]}` 或 `values: [...]` |
| `sweep2d` | 无 | 二维扫描：`{group: .., x: {layer, range}, y: {layer, range}}` |

### 材料定义（三种方式）

```yaml
materials:
  MyPolymer: my_nk/xxx.txt              # 数据文件：波长, n, k 三列
  PDINN: {n: 1.70, k: 0}                # 常数折射率
  ZnO: {cauchy: [1.90, 0.02, 0.001], k: 0}   # n = A + B/λ² + C/λ⁴（λ 单位 µm）
```

- 把 `XXX.csv` 放进 `materials/`，配置里直接写 `material: XXX` 就能用，不需要在 `materials:` 里登记。
- 数据文件的分隔符可以是逗号、Tab 或空格；表头会自动跳过；波长单位是 nm 还是 µm 会自动识别。

## 输出内容

| 文件 | 内容 |
|---|---|
| `E2_map_<λ>nm.png` | 各组 \|E\|² 截面图（x–z 图，红色虚线标界面，标出活性层最大值） |
| `E2_map_AM1.5G_weighted.png` | 按 AM1.5G 光子通量加权的平均 \|E\|² 截面图 |
| `E2_profile_z.png` | 沿 z 的 \|E\|² 剖面对比 |
| `absorption_reflection.png` | 活性层吸收谱（图例中给出 Jsc）和器件反射谱 |
| `delta_absorption.png` | 各组相对第一组的吸收差 |
| `generation_rate.png` | 活性层内产生率 G(z) |
| `loss_breakdown.png` | 光子去向分解：活性层 / 各层寄生吸收 / 反射 / 透射 (mA/cm²) |
| `nk_materials.png` | 所用材料的 n,k |
| `sweep_<层>.png/csv` | Jsc 随厚度变化 |
| `sweep2d_<层>_vs_<层>.png/csv` | 二维厚度扫描的 Jsc 等高图 |
| `summary.csv` | 每组的 Jsc、平均吸收、平均 \|E\|²、各项损失、能量守恒误差 |
| `absorption_all_groups.csv` | 各组的活性层吸收和反射 |
| `spectra_<组>.csv` | 该组每一层的吸收谱 A(λ)，以及 R、T |
| `field_<组>.csv` | 该组的 z、层名、各波长 \|E\|²、加权 \|E\|²、G(z) |
| `used_config.yaml` | 本次运行所用配置的备份（便于复现） |

CSV 文件用 UTF-8-BOM 编码，Excel 和 Origin 可以直接打开。

## 计算方法与假设

- 传输矩阵法（Pettersson, *J. Appl. Phys.* 86, 487 (1999)），光正入射。平面多层膜在正入射下的结果与 FDTD 一致。截面图沿 x 方向是均匀的，x 坐标只用于显示。
- |E|²/|E₀|² 以衬底中的入射场归一化，已计入空气/玻璃前表面的损失。
- 吸收功率密度 Q(z) = α·n·|E|²/n₀，α = 4πk/λ。Jsc,max = q∫A(λ)Φ(λ)dλ × IQE。
- 每次运行都检查能量守恒（R + T + ΣA = 1），误差列在 summary 中，一般 < 1e-4。
- **适用范围：** 只适用于平面膜层。陷光结构、纳米颗粒、粗糙界面、斜入射都不适用，这些情况需要 FDTD。

## 材料数据来源

| 材料 | 来源 |
|---|---|
| Al | McPeak et al., *ACS Photonics* 2, 326 (2015)，refractiveindex.info |
| ITO | Minenkov et al., *ACS AMI* 16, 9517 (2024)，玻璃上 110 nm ITO，refractiveindex.info |
| glass | Rubin, *Sol. Energy Mater.* 12, 275 (1985)，钠钙玻璃，refractiveindex.info |
| PM6Y6 | **模型数据**：KK 自洽的高斯振子模型，拟合文献典型 k 曲线（`tools/build_model_nk.py`） |
| SPT | **模型数据**：磷钨酸钠 NaPW12，可见/近红外透明，n ≈ 1.70，只有紫外吸收 |

有椭偏实测数据时，直接覆盖 `materials/PM6Y6.csv`、`materials/SPT.csv` 即可。
