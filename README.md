# 智理污水 · 数字孪生仿真系统

<p align="center">
  <img src="https://img.shields.io/badge/Python-3.8+-3776AB?style=flat-square&logo=python&logoColor=white"/>
  <img src="https://img.shields.io/badge/matplotlib-3.5+-11557C?style=flat-square&logo=python&logoColor=white"/>
  <img src="https://img.shields.io/badge/scipy-1.7+-8CAAE6?style=flat-square"/>
  <img src="https://img.shields.io/badge/License-MIT-green?style=flat-square"/>
  <img src="https://img.shields.io/badge/标准-GB_39731--2020-blue?style=flat-square"/>
</p>

<p align="center">
  <b>制药废水处理数字孪生平台 · 电芬顿氧化 vs 臭氧氧化 · 实时仿真对比</b>
</p>

---

## 📖 项目简介

本项目是**智能化新型污水处理塔**的数字孪生仿真模块，针对制药废水中的重金属有机络合物（如 Cu(II)-EDTA）去除问题，基于 ODE 动力学方程，实时对比两种核心处理工艺的去除效率：

| 系统 | 工艺 | 核心机制 |
|------|------|----------|
| **主系统** | 电芬顿氧化（Electro-Fenton） | 电催化生成 •OH 自由基，链式氧化降解 |
| **对比系统** | 臭氧氧化（Ozone Oxidation） | 臭氧直接氧化 + 间接羟基自由基路径 |

仿真结果严格对照 **GB 39731-2020** 现行排放标准，并通过自适应控制逻辑模拟真实工况下的 PID 调节过程。

> 🏫 **项目来源**：华中科技大学能源与动力工程学院 · 源梦小组  
> 👥 **设计者**：刘松奎、王赞宇、郑崇勋、郭子鑫、阮天佑  
> 📋 **指导老师**：罗光前教授（华中科技大学能源与动力工程学院）

---

## ✨ 功能特性

- **🔬 双工艺实时对比**：电芬顿与臭氧氧化并行仿真，参数独立演化
- **📊 六联动可视化面板**：趋势曲线 / 综合性能雷达图 / 处理效率柱状图 / 传感器数据面板
- **🤖 自适应控制仿真**：模拟传感器反馈驱动的电压、臭氧浓度、MVR 转速调节
- **⚡ 随机工况扰动**：注入高斯噪声模拟电网波动和气源压力变化
- **📏 国标对标**：实时标注 GB 39731-2020 排放限值参考线
- **🎨 专业暗色主题 UI**：GitHub 风格深色配色，适合工程演示和答辩展示

---

## 🏗️ 系统架构

```
制药废水
    │
    ▼
┌─────────────────────────────────────┐
│          预处理层                     │
│  还原破络（多羟基亚铁）+ 电芬顿氧化    │
└──────────────┬──────────────────────┘
               │
               ▼
┌─────────────────────────────────────┐
│          中处理层                     │
│  多级离子沉淀（Cu / Ni / Pb / Cr）   │
│  分子筛负压过滤                       │
└──────────────┬──────────────────────┘
               │
               ▼
┌─────────────────────────────────────┐
│          后处理层                     │
│  MVR 蒸发结晶（磁浮式压缩机）         │
│  NaCl / Na₂SO₄ 提纯回收              │
└─────────────────────────────────────┘
               │
               ▼
         达标排放 / 回流
```

数字孪生模块对上述完整流程进行建模，重点仿真**预处理层**的双工艺动力学。

---

## 📦 安装与运行

### 环境要求

- Python ≥ 3.8
- pip（或 conda）

### 快速开始

```bash
# 1. 克隆仓库
git clone https://github.com/your-username/wastewater-digital-twin.git
cd wastewater-digital-twin

# 2. 安装依赖
pip install -r requirements.txt

# 3. 运行仿真
python digital_twin.py
```

### Jupyter Notebook 中运行

```python
%matplotlib notebook
from digital_twin import WastewaterDigitalTwin, ProcessConfig
from matplotlib.animation import FuncAnimation

cfg = ProcessConfig()
twin = WastewaterDigitalTwin(cfg)
ani = FuncAnimation(twin.fig, twin.update, frames=cfg.TIME_STEPS,
                    interval=cfg.ANIMATION_INTERVAL_MS, blit=False)
```

---

## ⚙️ 参数配置

所有可调参数集中在 `ProcessConfig` 类中，无需修改核心逻辑：

```python
class ProcessConfig:
    # 初始水质 (mg/L)
    INIT_ORGANIC          = 100.0   # 有机物（COD）
    INIT_HEAVY_METAL      = 50.0    # 重金属浓度
    INIT_SUSPENDED_SOLIDS = 80.0    # 悬浮固体（SS）
    INIT_SALINITY         = 5.0     # 盐度 (g/L)

    # 电芬顿参数
    EF_VOLTAGE_INIT       = 12.0    # 初始电压 (V)
    EF_VOLTAGE_MIN        = 10.0    # 最低电压
    EF_VOLTAGE_MAX        = 15.0    # 最高电压

    # 臭氧参数
    OZONE_CONC_INIT       = 20.0    # 初始浓度 (mg/L)

    # MVR 压缩机
    MVR_RPM_INIT          = 3000    # 初始转速 (rpm)

    # 仿真控制
    TIME_STEPS            = 120     # 总时间步数
    ANIMATION_INTERVAL_MS = 180     # 帧间隔（越小越快）
```

---

## 📐 动力学模型

### 电芬顿氧化（主系统）

$$\frac{d[\text{Organic}]}{dt} = -0.15 \cdot [\text{Organic}] \cdot \frac{U}{10}$$

$$\frac{d[\text{Metal}]}{dt} = -0.08 \cdot [\text{Metal}] \cdot \frac{U}{10}$$

### 臭氧氧化（对比系统）

$$\frac{d[\text{Organic}]}{dt} = -0.10 \cdot [\text{Organic}] \cdot \frac{C_{O_3}}{20}$$

$$\frac{d[\text{Metal}]}{dt} = -0.04 \cdot [\text{Metal}] \cdot \frac{C_{O_3}}{20}$$

### 中后处理（共用）

$$\frac{d[\text{SS}]}{dt} = 0.05 \cdot [\text{Organic}] - 0.10 \cdot [\text{SS}]$$

$$\frac{d[\text{Salinity}]}{dt} = -0.02 \cdot [\text{Salinity}] \cdot \frac{n_{MVR}}{3000}$$

> 方程组通过 `scipy.integrate.odeint` 求解，每步使用 LSODA 自适应步长积分器。

---

## 📊 可视化说明

| 面板 | 内容 |
|------|------|
| 有机物浓度变化 | 双系统 COD 趋势 + GB 39731-2020 国标限值 |
| 重金属浓度变化 | 双系统重金属趋势 + 国标限值 |
| 悬浮固体变化 | SS 浓度趋势（中处理层代理指标） |
| 传感器实时数据 | 各水质参数 + 设备控制参数 + 达标状态指示 |
| 综合性能雷达图 | 五维度对比（有机物 / 重金属 / SS / 盐度 / 稳定性） |
| 实时效率柱状图 | 有机物 / 重金属去除率实时对比 |

---

## 💾 保存动画

取消 `main()` 中的注释，并确保已安装 `ffmpeg`：

```bash
# 安装 ffmpeg（macOS）
brew install ffmpeg

# 安装 ffmpeg（Ubuntu/Debian）
sudo apt install ffmpeg
```

```python
ani.save("wastewater_twin.mp4", writer="ffmpeg", fps=6, dpi=120)
```

---

## 📁 项目结构

```
wastewater-digital-twin/
├── digital_twin.py          # 主仿真程序（数字孪生核心）
├── requirements.txt         # Python 依赖
├── README.md                # 项目说明（本文件）
├── docs/
│   ├── design_report.pdf    # 智能化新型污水处理塔设计说明书
│   └── algorithm_diagram.vsdx  # 数字孪生算法架构图
└── assets/
    └── demo.gif             # 仿真效果演示
```

---

## 🔬 技术依据

| 创新点 | 参考文献 |
|--------|----------|
| 多羟基亚铁还原破络动力学 | 邰伟等, 能源环境保护, 2024 |
| 电芬顿法处理 Ni-EDTA | 邵天元等, 环境科学学报, 2015(3):745-749 |
| MVR 离心压缩机匹配特性 | 周东等, 工程设计学报, 2022, 29(5):595-606 |

---

## 📄 许可证

本项目采用 [MIT License](LICENSE) 开源，欢迎学术引用与二次开发。

---

## 🤝 联系方式

- **联系人**：刘松奎（大一本科生）
- **单位**：华中科技大学网络空间安全学院
- **邮箱**：3575365448@qq.com
