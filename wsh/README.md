# wsh · 有机太阳能电池光学仿真

```
wsh/
├─ 01_SPT-ETL光学仿真结果/     ← 📊 结果：ITO/SPT/PM6:Y6/ETL/Al，ETL = 无 / SPT 5 nm / SPT 8 nm
│   ├─ README_结果说明.md        结论、参数、先看哪张图
│   ├─ 00_配置备份/
│   ├─ 01_电场分布图_E2/
│   ├─ 02_吸收与反射谱/
│   ├─ 03_产生率与损失分析/
│   ├─ 04_厚度扫描/
│   ├─ 05_数据表格_csv/
│   └─ 06_材料nk/
│
├─ 02_光学仿真程序_optsim/     ← 🛠 可复用程序：自己定义器件结构，一键出图和数据
│   ├─ 使用说明.md               所有参数说明
│   ├─ 运行仿真.bat              Windows 双击运行，或把 yaml 拖上来运行
│   ├─ optsim.py                 主程序
│   ├─ 1_配置文件_configs/       ★ 在这里写器件（复制 template.yaml）
│   ├─ 2_材料库_materials/       n,k 数据
│   ├─ 3_光谱_spectra/           AM1.5G 光谱
│   ├─ 4_工具_tools/             模型 n,k 生成脚本
│   └─ 5_输出结果_results/       新器件的结果默认输出到这里
│
└─ 03_参考图_原始/             ← 🖼 原始参考图（FDTD |E|² 图）
    ├─ 4.7.tif
    └─ 4.7_预览.png
```

## 我想……

| 需求 | 去哪里 |
|---|---|
| 看 SPT-ETL 三组器件的结论 | `01_SPT-ETL光学仿真结果/README_结果说明.md` |
| 拿图写论文 | `01_.../01_电场分布图_E2/`、`02_吸收与反射谱/` |
| 拿数据自己用 Origin 作图 | `01_.../05_数据表格_csv/` |
| 仿真新器件 | 复制 `02_.../1_配置文件_configs/template.yaml` → 修改 → 运行 `python optsim.py 1_配置文件_configs/你的文件.yaml` |
| 换成自己测的 n,k | 把 csv 放进 `02_.../2_材料库_materials/`（三列：波长, n, k） |

## 安装（第一次使用）

```bash
pip install -r 02_光学仿真程序_optsim/requirements.txt
```
