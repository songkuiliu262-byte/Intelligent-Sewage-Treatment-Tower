"""
智理污水 · 数字孪生仿真系统
Zhili Wastewater · Digital Twin Simulation System

基于电芬顿氧化 vs 臭氧氧化的制药废水处理对比模拟
Pharmaceutical wastewater treatment comparison: Electro-Fenton vs Ozone Oxidation

Authors: 刘松奎, 王赞宇, 郑崇勋, 郭子鑫, 阮天佑
Institution: 华中科技大学能源与动力工程学院
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.animation import FuncAnimation
from scipy.integrate import odeint
from collections import deque
import warnings

warnings.filterwarnings("ignore")

# ── 字体配置 ──────────────────────────────────────────────────────────────────
plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.facecolor"] = "#0d1117"
plt.rcParams["axes.facecolor"] = "#161b22"
plt.rcParams["axes.edgecolor"] = "#30363d"
plt.rcParams["axes.labelcolor"] = "#c9d1d9"
plt.rcParams["xtick.color"] = "#8b949e"
plt.rcParams["ytick.color"] = "#8b949e"
plt.rcParams["text.color"] = "#c9d1d9"
plt.rcParams["grid.color"] = "#21262d"
plt.rcParams["grid.linestyle"] = "--"
plt.rcParams["grid.alpha"] = 0.6

# ── 颜色主题 ──────────────────────────────────────────────────────────────────
COLOR = {
    "electro_primary": "#58a6ff",   # 电芬顿 - 蓝色
    "electro_second":  "#1f6feb",
    "ozone_primary":   "#3fb950",   # 臭氧 - 绿色
    "ozone_second":    "#238636",
    "warning":         "#f85149",
    "accent":          "#d2a8ff",
    "text_muted":      "#8b949e",
    "panel_bg":        "#161b22",
}

# ── 国标排放限值 (GB 39731-2020, mg/L) ───────────────────────────────────────
DISCHARGE_STANDARD = {
    "organic":  30.0,   # COD 限值（近似）
    "heavy_metal": 0.5, # 总铜限值
    "suspended_solids": 70.0,
}


class ProcessConfig:
    """工艺参数配置类"""

    # 初始水质参数 (mg/L)
    INIT_ORGANIC = 100.0
    INIT_HEAVY_METAL = 50.0
    INIT_SUSPENDED_SOLIDS = 80.0
    INIT_SALINITY = 5.0       # g/L
    INIT_PH = 7.0

    # 电芬顿工艺参数
    EF_VOLTAGE_INIT = 12.0    # V
    EF_VOLTAGE_MIN = 10.0
    EF_VOLTAGE_MAX = 15.0

    # 臭氧氧化参数
    OZONE_CONC_INIT = 20.0    # mg/L
    OZONE_CONC_MIN = 10.0
    OZONE_CONC_MAX = 30.0

    # MVR 压缩机参数
    MVR_RPM_INIT = 3000       # rpm
    MVR_RPM_MIN = 2500
    MVR_RPM_MAX = 3500

    # 模拟参数
    TIME_STEPS = 120
    HISTORY_SIZE = 60
    ANIMATION_INTERVAL_MS = 180  # 帧间隔


class WastewaterDigitalTwin:
    """
    制药废水处理数字孪生仿真器

    模拟两种处理工艺的动态响应：
      - 主系统：电芬顿氧化（Electro-Fenton, EF）
      - 对比系统：臭氧氧化（Ozone Oxidation）
    """

    def __init__(self, cfg: ProcessConfig = None):
        self.cfg = cfg or ProcessConfig()
        self._init_state()
        self._init_history()
        self._build_figure()

    # ── 初始化 ────────────────────────────────────────────────────────────────

    def _init_state(self):
        """初始化系统状态变量"""
        c = self.cfg
        # 主系统（电芬顿）
        self.organic = c.INIT_ORGANIC
        self.heavy_metal = c.INIT_HEAVY_METAL
        self.suspended_solids = c.INIT_SUSPENDED_SOLIDS
        self.salinity = c.INIT_SALINITY
        self.ph = c.INIT_PH

        # 对比系统（臭氧）
        self.oz_organic = c.INIT_ORGANIC
        self.oz_heavy_metal = c.INIT_HEAVY_METAL

        # 设备参数
        self.ef_voltage = c.EF_VOLTAGE_INIT
        self.ozone_conc = c.OZONE_CONC_INIT
        self.mvr_rpm = c.MVR_RPM_INIT

        self.current_step = 0

    def _init_history(self):
        """初始化历史数据队列"""
        sz = self.cfg.HISTORY_SIZE
        self.t_hist = deque(maxlen=sz)
        # 电芬顿历史
        self.ef_organic_hist = deque(maxlen=sz)
        self.ef_metal_hist = deque(maxlen=sz)
        self.ef_solids_hist = deque(maxlen=sz)
        self.ef_salinity_hist = deque(maxlen=sz)
        # 臭氧历史
        self.oz_organic_hist = deque(maxlen=sz)
        self.oz_metal_hist = deque(maxlen=sz)

    def _build_figure(self):
        """构建可视化界面"""
        self.fig = plt.figure(figsize=(14, 9))
        self.fig.patch.set_facecolor("#0d1117")

        gs = gridspec.GridSpec(
            2, 3,
            figure=self.fig,
            hspace=0.45,
            wspace=0.35,
            left=0.07, right=0.97,
            top=0.90, bottom=0.08,
        )

        self.ax_organic = self.fig.add_subplot(gs[0, 0])
        self.ax_metal   = self.fig.add_subplot(gs[0, 1])
        self.ax_solids  = self.fig.add_subplot(gs[0, 2])
        self.ax_params  = self.fig.add_subplot(gs[1, 0])
        self.ax_radar   = self.fig.add_subplot(gs[1, 1], polar=True)
        self.ax_eff     = self.fig.add_subplot(gs[1, 2])

        self.fig.suptitle(
            "智理污水  ·  数字孪生仿真平台\n"
            "Pharmaceutical Wastewater Digital Twin  —  Electro-Fenton  vs  Ozone",
            fontsize=13,
            color="#c9d1d9",
            y=0.97,
        )

    # ── 动力学方程 ────────────────────────────────────────────────────────────

    def _dynamics(self, y, t, method: str = "electro"):
        """
        ODE 系统动力学

        Args:
            y: [organic, heavy_metal, suspended_solids, salinity]
            t: 时间
            method: "electro" | "ozone"
        """
        organic, metal, solids, salinity = y

        if method == "electro":
            # 电芬顿：电压驱动的高级氧化（羟基自由基链式反应）
            k_ef = self.ef_voltage / 10.0
            d_organic = -0.15 * organic * k_ef
            d_metal   = -0.08 * metal   * k_ef
        else:
            # 臭氧氧化：臭氧浓度驱动的直接氧化
            k_oz = self.ozone_conc / 20.0
            d_organic = -0.10 * organic * k_oz
            d_metal   = -0.04 * metal   * k_oz

        # 中后处理共用动力学
        d_solids   = 0.05 * organic - 0.10 * solids
        d_salinity = -0.02 * salinity * (self.mvr_rpm / 3000.0)

        return [d_organic, d_metal, d_solids, d_salinity]

    # ── 系统更新 ──────────────────────────────────────────────────────────────

    def _step_ozone(self):
        """推进臭氧对比系统一个时间步"""
        t = np.linspace(0, 1, 2)
        y0 = [self.oz_organic, self.oz_heavy_metal,
              self.suspended_solids, self.salinity]
        y = odeint(self._dynamics, y0, t, args=("ozone",))
        self.oz_organic    = max(0.0, y[-1, 0])
        self.oz_heavy_metal = max(0.0, y[-1, 1])
        # 臭氧系统工况扰动（模拟气源压力波动）
        if np.random.random() < 0.08:
            self.oz_organic     = max(0.0, self.oz_organic     + np.random.normal(0, 4))
            self.oz_heavy_metal = max(0.0, self.oz_heavy_metal + np.random.normal(0, 1.5))

    def _step_electro(self):
        """推进电芬顿主系统一个时间步"""
        t = np.linspace(0, 1, 2)
        y0 = [self.organic, self.heavy_metal,
              self.suspended_solids, self.salinity]
        y = odeint(self._dynamics, y0, t, args=("electro",))
        self.organic          = max(0.0, y[-1, 0])
        self.heavy_metal      = max(0.0, y[-1, 1])
        self.suspended_solids = max(0.0, y[-1, 2])
        self.salinity         = max(0.1, y[-1, 3])
        self.ph = 6.5 + 0.5 * np.sin(self.current_step * 0.12)
        # 电气扰动（模拟电网波动）
        if np.random.random() < 0.08:
            self.organic     = max(0.0, self.organic     + np.random.normal(0, 4))
            self.heavy_metal = max(0.0, self.heavy_metal + np.random.normal(0, 1.5))

    def _adaptive_control(self):
        """基于传感器反馈的自适应控制逻辑（PID 简化版）"""
        c = self.cfg
        # 电压控制：有机物浓度高则升压
        if self.organic > 60:
            self.ef_voltage  = min(c.EF_VOLTAGE_MAX,  self.ef_voltage  + 0.12)
            self.ozone_conc  = min(c.OZONE_CONC_MAX,  self.ozone_conc  + 0.5)
        else:
            self.ef_voltage  = max(c.EF_VOLTAGE_MIN,  self.ef_voltage  - 0.06)
            self.ozone_conc  = max(c.OZONE_CONC_MIN,  self.ozone_conc  - 0.25)
        # MVR 转速控制：盐度高则加速蒸发
        if self.salinity > 3.0:
            self.mvr_rpm = min(c.MVR_RPM_MAX, self.mvr_rpm + 12)
        else:
            self.mvr_rpm = max(c.MVR_RPM_MIN, self.mvr_rpm - 6)

    def _record_history(self):
        """记录当前步历史数据"""
        self.t_hist.append(self.current_step)
        self.ef_organic_hist.append(self.organic)
        self.ef_metal_hist.append(self.heavy_metal)
        self.ef_solids_hist.append(self.suspended_solids)
        self.ef_salinity_hist.append(self.salinity)
        self.oz_organic_hist.append(self.oz_organic)
        self.oz_metal_hist.append(self.oz_heavy_metal)

    def step(self, frame: int):
        """执行一个完整仿真步"""
        self.current_step = frame
        self._step_ozone()
        self._step_electro()
        self._adaptive_control()
        self._record_history()

    # ── 绘图辅助 ──────────────────────────────────────────────────────────────

    @staticmethod
    def _style_ax(ax, title: str, xlabel: str, ylabel: str):
        ax.set_title(title, fontsize=10, color="#c9d1d9", pad=6)
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.grid(True)
        ax.legend(fontsize=8, loc="upper right",
                  facecolor="#21262d", edgecolor="#30363d", labelcolor="#c9d1d9")

    def _draw_trend(self, ax, title: str, ylabel: str,
                    ef_data, oz_data, std_limit=None):
        """绘制双系统趋势对比折线"""
        t = list(self.t_hist)
        ax.clear()
        ax.plot(t, list(ef_data), color=COLOR["electro_primary"],
                linewidth=1.8, label="电芬顿 EF")
        ax.plot(t, list(oz_data), color=COLOR["ozone_primary"],
                linewidth=1.8, linestyle="--", label="臭氧 O₃")
        if std_limit is not None:
            ax.axhline(y=std_limit, color=COLOR["warning"],
                       linewidth=1.2, linestyle=":", label=f"国标限值 {std_limit}")
        self._style_ax(ax, title, "时间步 / step", ylabel)

    # ── 效率雷达图 ────────────────────────────────────────────────────────────

    def _draw_radar(self, ax):
        ax.clear()
        if len(self.ef_organic_hist) < 5:
            return

        categories = ["有机物\n去除率", "重金属\n去除率", "悬浮物\n去除率", "盐度\n削减率", "系统\n稳定性"]
        n = len(categories)

        def removal(init, current):
            return max(0.0, (init - current) / init * 100) if init > 0 else 0

        cfg = self.cfg
        ef_vals = [
            removal(cfg.INIT_ORGANIC,          self.ef_organic_hist[-1]),
            removal(cfg.INIT_HEAVY_METAL,       self.ef_metal_hist[-1]),
            removal(cfg.INIT_SUSPENDED_SOLIDS,  self.ef_solids_hist[-1]),
            removal(cfg.INIT_SALINITY,          self.ef_salinity_hist[-1]),
            85.0,  # 稳定性（固定评分）
        ]
        oz_vals = [
            removal(cfg.INIT_ORGANIC,      self.oz_organic_hist[-1]),
            removal(cfg.INIT_HEAVY_METAL,  self.oz_metal_hist[-1]),
            removal(cfg.INIT_SUSPENDED_SOLIDS, self.ef_solids_hist[-1]) * 0.85,
            removal(cfg.INIT_SALINITY,     self.ef_salinity_hist[-1]) * 0.9,
            72.0,
        ]

        angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist()
        angles += angles[:1]
        ef_vals  += ef_vals[:1]
        oz_vals  += oz_vals[:1]

        ax.set_facecolor("#161b22")
        ax.plot(angles, ef_vals, color=COLOR["electro_primary"], linewidth=2)
        ax.fill(angles, ef_vals, color=COLOR["electro_primary"], alpha=0.25)
        ax.plot(angles, oz_vals, color=COLOR["ozone_primary"], linewidth=2, linestyle="--")
        ax.fill(angles, oz_vals, color=COLOR["ozone_primary"], alpha=0.15)

        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories, size=8, color="#c9d1d9")
        ax.set_ylim(0, 100)
        ax.set_yticks([25, 50, 75, 100])
        ax.set_yticklabels(["25%", "50%", "75%", "100%"], size=6, color="#8b949e")
        ax.tick_params(colors="#8b949e")
        ax.set_title("综合性能雷达图", fontsize=10, color="#c9d1d9", pad=15)
        ax.spines["polar"].set_color("#30363d")
        ax.grid(color="#30363d")

    # ── 效率柱状对比图 ────────────────────────────────────────────────────────

    def _draw_efficiency(self, ax):
        ax.clear()
        if len(self.ef_organic_hist) < 5:
            return

        cfg = self.cfg

        def eff(init, hist):
            return max(0.0, (init - hist[-1]) / init * 100) if init > 0 else 0

        ef_eff = [
            eff(cfg.INIT_ORGANIC,     self.ef_organic_hist),
            eff(cfg.INIT_HEAVY_METAL, self.ef_metal_hist),
        ]
        oz_eff = [
            eff(cfg.INIT_ORGANIC,     self.oz_organic_hist),
            eff(cfg.INIT_HEAVY_METAL, self.oz_metal_hist),
        ]

        x = np.arange(2)
        w = 0.32
        bars1 = ax.bar(x - w / 2, ef_eff, w,
                       color=COLOR["electro_primary"], label="电芬顿 EF", alpha=0.85)
        bars2 = ax.bar(x + w / 2, oz_eff, w,
                       color=COLOR["ozone_primary"],   label="臭氧 O₃",  alpha=0.85)

        for bar in list(bars1) + list(bars2):
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.8,
                    f"{h:.1f}%", ha="center", va="bottom", fontsize=8, color="#c9d1d9")

        ax.set_xticks(x)
        ax.set_xticklabels(["有机物去除率", "重金属去除率"], fontsize=9)
        ax.set_ylim(0, 110)
        ax.set_ylabel("去除效率 (%)", fontsize=8)
        ax.set_title("实时处理效率对比", fontsize=10, color="#c9d1d9", pad=6)
        ax.legend(fontsize=8, facecolor="#21262d", edgecolor="#30363d", labelcolor="#c9d1d9")
        ax.axhline(y=90, color=COLOR["warning"], linewidth=1, linestyle=":",
                   label="达标线 90%")
        ax.grid(True, axis="y")

    # ── 参数面板 ──────────────────────────────────────────────────────────────

    def _draw_params(self, ax):
        ax.clear()
        ax.set_facecolor("#0d1117")
        ax.axis("off")

        def status_icon(val, threshold):
            return "●" if val <= threshold else "⚠"

        lines = [
            ("── 传感器实时数据 ──", COLOR["accent"],  12),
            (f"有机物浓度      {self.organic:.1f} mg/L  {status_icon(self.organic, 30)}",
             COLOR["electro_primary"] if self.organic < 60 else COLOR["warning"], 10),
            (f"重金属浓度      {self.heavy_metal:.2f} mg/L  {status_icon(self.heavy_metal, 0.5)}",
             COLOR["electro_primary"] if self.heavy_metal < 25 else COLOR["warning"], 10),
            (f"悬浮固体        {self.suspended_solids:.1f} mg/L", COLOR["text_muted"], 10),
            (f"盐度            {self.salinity:.2f} g/L",           COLOR["text_muted"], 10),
            (f"pH 值           {self.ph:.2f}",                      COLOR["text_muted"], 10),
            ("", "#ffffff", 8),
            ("── 设备控制参数 ──", COLOR["accent"],  12),
            (f"电芬顿电压      {self.ef_voltage:.1f} V",  COLOR["electro_primary"], 10),
            (f"臭氧浓度        {self.ozone_conc:.1f} mg/L", COLOR["ozone_primary"],   10),
            (f"MVR 转速        {self.mvr_rpm:.0f} rpm",    COLOR["text_muted"], 10),
            ("", "#ffffff", 8),
            (f"仿真进度  {self.current_step + 1}/{self.cfg.TIME_STEPS} 步",
             COLOR["text_muted"], 9),
        ]

        y = 0.96
        for text, color, size in lines:
            ax.text(0.05, y, text, transform=ax.transAxes,
                    fontsize=size, color=color, verticalalignment="top",
                    fontfamily="monospace")
            y -= 0.075

    # ── 主更新函数 ────────────────────────────────────────────────────────────

    def update(self, frame: int):
        """动画帧更新回调"""
        self.step(frame)

        if len(self.t_hist) < 2:
            return

        # 趋势图
        self._draw_trend(
            self.ax_organic, "有机物浓度变化 (COD)", "浓度 (mg/L)",
            self.ef_organic_hist, self.oz_organic_hist,
            std_limit=DISCHARGE_STANDARD["organic"],
        )
        self._draw_trend(
            self.ax_metal, "重金属浓度变化", "浓度 (mg/L)",
            self.ef_metal_hist, self.oz_metal_hist,
            std_limit=DISCHARGE_STANDARD["heavy_metal"],
        )
        self._draw_trend(
            self.ax_solids, "悬浮固体变化 (SS)", "浓度 (mg/L)",
            self.ef_solids_hist, self.ef_solids_hist,  # SS 仅主系统
        )

        self._draw_params(self.ax_params)
        self._draw_radar(self.ax_radar)
        self._draw_efficiency(self.ax_eff)

        return []


# ── 入口 ──────────────────────────────────────────────────────────────────────

def main():
    cfg = ProcessConfig()
    twin = WastewaterDigitalTwin(cfg)

    ani = FuncAnimation(
        twin.fig,
        twin.update,
        frames=cfg.TIME_STEPS,
        interval=cfg.ANIMATION_INTERVAL_MS,
        blit=False,
        repeat=False,
    )

    plt.show()

    # 保存动画（需要 ffmpeg）：
    # ani.save("wastewater_twin.mp4", writer="ffmpeg", fps=6, dpi=120)


if __name__ == "__main__":
    main()
