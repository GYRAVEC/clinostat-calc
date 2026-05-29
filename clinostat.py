#!/usr/bin/env python3
"""
2軸クリノスタット 回転設定計算ツール

  ハードウェア: 内側フレーム r = 8.5 cm / 外側フレーム r = 14.5 cm

  サンプル配置とモードの組み合わせは3通り:

  Tilt モード (既定, サンプル=2軸交点 r≈0, 0 ≤ g_target ≤ 1g):
    外軸を鉛直から ψ = arccos(g_target) 傾け、ω_out 高速 + ω_in 低速で
    細胞は |g| = g0·cos ψ を短時間平均として知覚し、長時間平均で蓄積ゼロ.
    懸濁状態を保ちたい培養に推奨.

  Centrifuge モード (--mode centrifuge, サンプル=4角 r=8.5cm, 0 ≤ g_target):
    内側回転 ω_in² × r で遠心力を g_target として負荷.
    外側はクリノローテーション (低速で地球重力を時間平均で消す).
    遠心方向に沈殿が起きる (Mars/Moon の実環境と同じ).
    Mars analog 実験など、沈殿が"再現すべき現象"の場合に使う.

  Switching モード (--mode switching, サンプル=2軸交点 r≈0, 0 ≤ g_target ≤ 0.5g):
    傾斜なし運用. Mode A (1:1+位相ドリフト 0.5g) と Mode B (黄金比 0g) を
    時間配分で混合する従来方式. 沈殿なしだが上限 0.5g.
"""

import argparse
import math
import sys

# 🐱 物理定数とフレーム諸元
G_EARTH = 9.81
R_INNER_FRAME = 0.085
R_OUTER_FRAME = 0.145
PHI = (1 + math.sqrt(5)) / 2

# 🐱 動作パラメータ既定値
DEFAULT_BASE_RPM = 5.0
DEFAULT_DRIFT_MIN = 20.0
DEFAULT_CYCLE_MIN = 60.0

# 🐱 細菌培養想定 (沈降診断)
V_SED_BACTERIA = 1.0e-6        # 自由沈降速度 [m/s]
R_CHAMBER_DEFAULT = 0.005      # チャンバ半径 [m]
T_BIO_BACTERIA_S = 600.0       # 生物応答時間 [s]

G_MAX_SWITCHING = 0.5          # スイッチングモードの達成上限
G_MAX_TILT = 1.0               # 傾斜モードの達成上限

PRESETS = {
    "earth":  ("地球",     1.000),
    "mars":   ("火星",     0.378),
    "moon":   ("月",       0.166),
    "ceres":  ("ケレス",   0.029),
    "europa": ("エウロパ", 0.134),
    "titan":  ("タイタン", 0.138),
    "iss":    ("ISS",      1.0e-4),
    "ug":     ("微小重力", 1.0e-3),
}


# ---------------- Tilt モード -----------------

def compute_tilt(g_target: float,
                 base_rpm: float = DEFAULT_BASE_RPM,
                 drift_min: float = DEFAULT_DRIFT_MIN,
                 chamber_r: float = R_CHAMBER_DEFAULT) -> dict:
    """外軸傾斜モード. ψ = arccos(g_target) で 0..1g 全域カバー."""
    if g_target < 0:
        raise ValueError("g_target は非負")
    if g_target > G_MAX_TILT + 1e-9:
        raise ValueError(f"g_target {g_target:g}g > 1g は地球重力超過のため不可")

    psi_rad = math.acos(min(g_target, 1.0))
    psi_deg = math.degrees(psi_rad)

    omega_out = base_rpm * 2 * math.pi / 60.0          # rad/s
    T_drift = drift_min * 60.0                          # s
    omega_in = 2 * math.pi / T_drift if T_drift > 0 else 0.0

    # 🐱 沈降診断: 短期平均 g_target でドリフト方向が回るので、振幅 = v_sed*g*T/π
    drift_amp = V_SED_BACTERIA * g_target * T_drift / math.pi
    chamber_ratio = drift_amp / chamber_r if chamber_r > 0 else float("inf")

    warnings = []
    if T_drift < T_BIO_BACTERIA_S:
        warnings.append(
            f"T_drift ({drift_min:.1f}min) < T_bio (~{T_BIO_BACTERIA_S/60:.0f}min): "
            "細胞は g_target を知覚せず ~0g として応答する可能性"
        )
    if drift_amp > chamber_r:
        warnings.append(
            f"ドリフト振幅 {drift_amp*1000:.2f} mm > チャンバ半径 {chamber_r*1000:.1f} mm: "
            "T_drift を短く / base-rpm を上げる必要"
        )
    if g_target > 0.9999 and base_rpm > 0:
        warnings.append("g_target ≈ 1g: ψ = 0 のため外軸回転は意味なし (傾斜なし静置でよい)")

    return {
        "mode": "tilt",
        "g_target": g_target,
        "g_effective": g_target,
        "psi_deg": psi_deg,
        "psi_rad": psi_rad,
        "base_rpm": base_rpm,
        "drift_min": drift_min,
        "chamber_r_mm": chamber_r * 1000.0,
        "omega_out_rpm": omega_out * 60.0 / (2 * math.pi),
        "omega_in_rpm": omega_in * 60.0 / (2 * math.pi),
        "omega_out_rad_s": omega_out,
        "omega_in_rad_s": omega_in,
        "diagnostics": {
            "drift_amplitude_mm": drift_amp * 1000.0,
            "chamber_safety_ratio": chamber_ratio,
        },
        "warnings": warnings,
    }


# ---------------- Switching モード (Mode A + Mode B) -----------------

def compute_switching(g_target: float,
                      base_rpm: float = DEFAULT_BASE_RPM,
                      drift_min: float = DEFAULT_DRIFT_MIN,
                      cycle_min: float = DEFAULT_CYCLE_MIN,
                      chamber_r: float = R_CHAMBER_DEFAULT) -> dict:
    """Mode A (1:1+位相ドリフト 0.5g) と Mode B (黄金比 0g) の時間配分."""
    if g_target < 0:
        raise ValueError("g_target は非負")
    if g_target > G_MAX_SWITCHING + 1e-9:
        raise ValueError(
            f"g_target {g_target:g}g > {G_MAX_SWITCHING}g は switching モードでは不可。"
            " --mode tilt を使うか base-rpm/drift-min を調整"
        )

    omega_base = base_rpm * 2 * math.pi / 60.0
    T_drift = drift_min * 60.0
    T_cycle = cycle_min * 60.0

    eps = 2 * math.pi / T_drift
    omega_out_A = omega_base
    omega_in_A = omega_base + eps
    omega_out_B = omega_base
    omega_in_B = omega_base / PHI

    duty_A = min(2.0 * g_target, 1.0)
    duty_B = 1.0 - duty_A

    drift_amp = V_SED_BACTERIA * 0.5 * T_drift / math.pi
    chamber_ratio = drift_amp / chamber_r if chamber_r > 0 else float("inf")

    warnings = []
    if T_drift < T_BIO_BACTERIA_S:
        warnings.append(
            f"T_drift ({drift_min:.1f}min) < T_bio (~{T_BIO_BACTERIA_S/60:.0f}min): "
            "細胞は 0.5g を知覚できず ~0g として応答する可能性"
        )
    if drift_amp > chamber_r:
        warnings.append(
            f"ドリフト振幅 {drift_amp*1000:.2f} mm > チャンバ半径 {chamber_r*1000:.1f} mm"
        )

    return {
        "mode": "switching",
        "g_target": g_target,
        "g_effective": 0.5 * duty_A,
        "base_rpm": base_rpm,
        "drift_min": drift_min,
        "cycle_min": cycle_min,
        "chamber_r_mm": chamber_r * 1000.0,
        "mode_A": {
            "name": "1:1+位相ドリフト",
            "duty": duty_A,
            "duration_min": duty_A * cycle_min,
            "omega_out_rpm": omega_out_A * 60.0 / (2 * math.pi),
            "omega_in_rpm": omega_in_A * 60.0 / (2 * math.pi),
            "omega_out_rad_s": omega_out_A,
            "omega_in_rad_s": omega_in_A,
            "epsilon_rad_s": eps,
            "g_perceived": 0.5,
        },
        "mode_B": {
            "name": "黄金比(微小重力)",
            "duty": duty_B,
            "duration_min": duty_B * cycle_min,
            "omega_out_rpm": omega_out_B * 60.0 / (2 * math.pi),
            "omega_in_rpm": omega_in_B * 60.0 / (2 * math.pi),
            "omega_out_rad_s": omega_out_B,
            "omega_in_rad_s": omega_in_B,
            "g_perceived": 0.0,
        },
        "diagnostics": {
            "drift_amplitude_mm": drift_amp * 1000.0,
            "chamber_safety_ratio": chamber_ratio,
        },
        "warnings": warnings,
    }


# ---------------- Centrifuge モード (サンプル=内側4角) -----------------

def compute_centrifuge(g_target: float,
                       sample_r: float = R_INNER_FRAME,
                       outer_rpm: float = DEFAULT_BASE_RPM,
                       chamber_r: float = R_CHAMBER_DEFAULT,
                       culture_hours: float = 24.0) -> dict:
    """サンプルを内軸から sample_r (既定 8.5 cm = 4角配置) に置き、
    内側遠心力 ω_in² × sample_r で g_target を負荷する.
    外側はクリノローテーション (低速) で地球重力を時間平均で消す.
    沈殿は遠心方向に必ず発生する (Mars/Moon の実環境と同じ性質).
    """
    if g_target < 0:
        raise ValueError("g_target は非負")
    if sample_r <= 0:
        raise ValueError(f"sample_r {sample_r} は正の値. r=0 の場合は tilt モード推奨")

    # 🐱 内側遠心で g_target 負荷: ω_in² × sample_r = g_target × g0
    omega_in = math.sqrt(g_target * G_EARTH / sample_r)
    rpm_in = omega_in * 60.0 / (2 * math.pi)

    # 🐱 外側は低速クリノローテーション (黄金比で非共鳴)
    omega_out = outer_rpm * 2 * math.pi / 60.0

    # 🐱 外側位置 (フレーム半径 R_OUTER_FRAME) での外乱遠心
    a_outer_disturb = omega_out ** 2 * R_OUTER_FRAME
    g_outer_disturb = a_outer_disturb / G_EARTH

    # 🐱 沈降診断: チャンバ底到達時間 = R_ch / (v_sed × g_target)
    sed_velocity = V_SED_BACTERIA * g_target            # [m/s]
    if sed_velocity > 0:
        t_reach_wall_s = chamber_r / sed_velocity
        t_reach_wall_h = t_reach_wall_s / 3600.0
    else:
        t_reach_wall_s = float("inf")
        t_reach_wall_h = float("inf")

    # 🐱 培養時間中、何時間ぶん細胞が壁に張り付いているか
    if t_reach_wall_h < culture_hours:
        wall_dwell_h = culture_hours - t_reach_wall_h
        wall_dwell_ratio = wall_dwell_h / culture_hours
    else:
        wall_dwell_h = 0.0
        wall_dwell_ratio = 0.0

    warnings = []
    if g_outer_disturb > g_target * 0.1 and g_target > 1e-4:
        warnings.append(
            f"外側遠心外乱 {g_outer_disturb:.4f}g が目標の10%超: "
            f"--outer-rpm を下げる"
        )
    if t_reach_wall_h < 1.0 and g_target > 1e-3:
        warnings.append(
            f"細胞は {t_reach_wall_h*60:.1f} 分で壁に到達: "
            "短時間培養しかできない or 大きいチャンバが必要"
        )
    if g_target > 5.0:
        warnings.append(
            f"g_target {g_target:g} は ω_in = {rpm_in:.0f} RPM 必要、装置の機械的限界に注意"
        )

    return {
        "mode": "centrifuge",
        "g_target": g_target,
        "g_effective": g_target,
        "sample_r_cm": sample_r * 100.0,
        "outer_rpm": outer_rpm,
        "chamber_r_mm": chamber_r * 1000.0,
        "culture_hours": culture_hours,
        "omega_in_rpm": rpm_in,
        "omega_in_rad_s": omega_in,
        "omega_out_rpm": outer_rpm,
        "omega_out_rad_s": omega_out,
        "diagnostics": {
            "outer_disturb_g": g_outer_disturb,
            "sed_velocity_um_s": sed_velocity * 1e6,
            "t_reach_wall_h": t_reach_wall_h,
            "wall_dwell_h": wall_dwell_h,
            "wall_dwell_ratio": wall_dwell_ratio,
        },
        "warnings": warnings,
    }


# ---------------- Hybrid モード (任意のサンプル半径) -----------------

def compute_hybrid(g_target: float,
                   sample_r: float = 0.02,
                   base_rpm: float = DEFAULT_BASE_RPM,
                   drift_min: float = DEFAULT_DRIFT_MIN,
                   chamber_r: float = R_CHAMBER_DEFAULT,
                   culture_hours: float = 24.0) -> dict:
    """サンプルが軸交点から r 離れた任意の位置にある場合のハイブリッド計算.

    戦略: Tilt モードベース (ω_in 低速ドリフト) + 遠心寄与の補正
    Body系での重力 2成分:
      g_tilt        = g0 · cos ψ        (ゆっくりドリフト方向)
      g_centrifugal = ω_in² × r          (body 固定方向)
    これらは直交するので合成は二乗和:
      |⟨g⟩| = √(g_tilt² + g_cent²) = g_target × g0
    """
    if g_target < 0 or g_target > G_MAX_TILT + 1e-9:
        raise ValueError(f"g_target {g_target}g は [0, {G_MAX_TILT}] の範囲")
    if sample_r < 0:
        raise ValueError(f"sample_r {sample_r} は非負")

    omega_out = base_rpm * 2 * math.pi / 60.0
    T_drift = drift_min * 60.0
    omega_in = 2 * math.pi / T_drift if T_drift > 0 else 0.0

    # 🐱 サンプル位置での遠心寄与 (body-fixed)
    g_cent_abs = omega_in ** 2 * sample_r            # [m/s²]
    g_cent_norm = g_cent_abs / G_EARTH               # [g]

    # 🐱 tilt 寄与を残りに割り当て (直交合成)
    g_tilt_normsq = g_target ** 2 - g_cent_norm ** 2
    if g_tilt_normsq < -1e-12:
        raise ValueError(
            f"遠心寄与 {g_cent_norm:.4e}g (= ω_in²·r/g0) > 目標 {g_target}g: "
            f"sample_r ({sample_r*100:.2f}cm) を下げるか drift_min を増やす"
        )
    g_tilt_norm = math.sqrt(max(g_tilt_normsq, 0.0))
    if g_tilt_norm > 1.0:
        raise ValueError(f"必要な tilt 寄与 {g_tilt_norm}g > 1g: 物理的に不可")
    psi_rad = math.acos(min(g_tilt_norm, 1.0))
    psi_deg = math.degrees(psi_rad)

    # 🐱 沈殿診断
    # tilt 寄与: 閉軌道 (1ドリフト周期で正味ゼロ), 振幅 v_sed*g_tilt*T/π
    drift_amp = V_SED_BACTERIA * g_tilt_norm * T_drift / math.pi      # m
    # 遠心寄与: body-fixed で時間に比例して累積
    v_cent = V_SED_BACTERIA * g_cent_norm                              # m/s
    cent_disp_total = v_cent * culture_hours * 3600.0                  # m
    total_disp = drift_amp + cent_disp_total
    chamber_ratio = total_disp / chamber_r if chamber_r > 0 else float("inf")

    warnings = []
    if g_cent_norm > 0.05 * g_target and g_target > 1e-3:
        warnings.append(
            f"遠心寄与 {g_cent_norm:.4e}g (= ω_in²·r/g0) が目標の 5% 超: "
            "sample_r を下げる or drift_min を増やす or centrifuge モードに切替"
        )
    if cent_disp_total > chamber_r:
        warnings.append(
            f"遠心由来累積変位 {cent_disp_total*1000:.2f}mm > チャンバ半径 {chamber_r*1000:.1f}mm: "
            f"培養 {culture_hours}h 内に壁到達"
        )
    if T_drift < T_BIO_BACTERIA_S:
        warnings.append(
            f"T_drift ({drift_min:.1f}min) < T_bio (~{T_BIO_BACTERIA_S/60:.0f}min): "
            "細胞は g_target を知覚できない可能性"
        )
    if sample_r == 0.0:
        warnings.append("sample_r = 0 は tilt モードと等価. --mode tilt 推奨")

    return {
        "mode": "hybrid",
        "g_target": g_target,
        "g_effective": g_target,
        "g_tilt_component": g_tilt_norm,
        "g_centrifugal_component": g_cent_norm,
        "sample_r_cm": sample_r * 100.0,
        "psi_deg": psi_deg,
        "psi_rad": psi_rad,
        "base_rpm": base_rpm,
        "drift_min": drift_min,
        "chamber_r_mm": chamber_r * 1000.0,
        "culture_hours": culture_hours,
        "omega_out_rpm": omega_out * 60.0 / (2 * math.pi),
        "omega_in_rpm": omega_in * 60.0 / (2 * math.pi),
        "omega_out_rad_s": omega_out,
        "omega_in_rad_s": omega_in,
        "diagnostics": {
            "tilt_drift_amp_mm": drift_amp * 1000.0,
            "centrifugal_displacement_mm": cent_disp_total * 1000.0,
            "total_displacement_mm": total_disp * 1000.0,
            "chamber_safety_ratio": chamber_ratio,
        },
        "warnings": warnings,
    }


# ---------------- 表示系 -----------------

ANSI = {
    "reset":   "\033[0m",
    "bold":    "\033[1m",
    "dim":     "\033[2m",
    "cyan":    "\033[36m",
    "yellow":  "\033[33m",
    "green":   "\033[32m",
    "magenta": "\033[35m",
    "red":     "\033[31m",
    "blue":    "\033[34m",
}


def c(text: str, *styles: str, use_color: bool = True) -> str:
    if not use_color:
        return text
    return "".join(ANSI[s] for s in styles) + text + ANSI["reset"]


def print_header(use_color: bool = True) -> None:
    line = "=" * 72
    print(c(line, "cyan", use_color=use_color))
    print(c("  2軸クリノスタット 回転設定計算ツール", "bold", "cyan", use_color=use_color))
    print(c(f"  内側フレーム r = {R_INNER_FRAME*100:.1f} cm  /  外側フレーム r = {R_OUTER_FRAME*100:.1f} cm",
            "dim", use_color=use_color))
    print(c(f"  g0 = {G_EARTH} m/s²  /  φ = {PHI:.6f}  /  "
            "サンプル位置はモードにより異なる (tilt/switching=軸交点, centrifuge=4角)",
            "dim", use_color=use_color))
    print(c(line, "cyan", use_color=use_color))


def render_tilt(r: dict, label: str = "", use_color: bool = True) -> None:
    head = f"目標重力: {r['g_target']:.4g} g  ({r['g_target'] * G_EARTH:.4f} m/s²)"
    if label:
        head = f"[{label}] " + head
    print(c(head + "   [Tilt モード]", "bold", "yellow", use_color=use_color))
    print(f"  実効時間平均 : {r['g_effective']:.4f} g  (= g0·cos ψ)")
    print(f"  位相ドリフト周期: {r['drift_min']:.1f} min")

    print(c("\n  [ハードウェア設定]", "bold", "blue", use_color=use_color))
    psi_str = c(f"{r['psi_deg']:8.4f}°", "bold", use_color=use_color)
    print(f"    外軸傾斜角 ψ : {psi_str}  ({r['psi_rad']:.6f} rad)  [鉛直から]")

    print(c("\n  [回転設定]", "bold", "green", use_color=use_color))
    out_str = c(f"{r['omega_out_rpm']:9.5f} RPM", "bold", use_color=use_color)
    in_str = c(f"{r['omega_in_rpm']:9.5f} RPM", "bold", use_color=use_color)
    print(f"    外側 ω_out : {out_str}  ({r['omega_out_rad_s']:.6f} rad/s)  [傾斜軸周り高速]")
    print(f"    内側 ω_in  : {in_str}  ({r['omega_in_rad_s']:.6e} rad/s)  [方向ドリフト用低速]")

    print(c(f"\n  [診断] 細菌想定 v_sed={V_SED_BACTERIA*1e6:.1f} μm/s, "
            f"チャンバ半径 {r['chamber_r_mm']:.1f} mm",
            "dim", use_color=use_color))
    d = r["diagnostics"]
    safety = d["chamber_safety_ratio"]
    print(f"    ドリフト振幅: {d['drift_amplitude_mm']:.3f} mm  ({safety*100:.1f}% of chamber)  "
          + (c("[!] チャンバ越え", "red", use_color=use_color) if safety > 1
             else c("[OK]", "dim", use_color=use_color)))

    for w in r["warnings"]:
        print("  " + c("[!] " + w, "red", use_color=use_color))


def render_switching(r: dict, label: str = "", use_color: bool = True) -> None:
    head = f"目標重力: {r['g_target']:.4g} g  ({r['g_target'] * G_EARTH:.4f} m/s²)"
    if label:
        head = f"[{label}] " + head
    print(c(head + "   [Switching モード]", "bold", "yellow", use_color=use_color))
    print(f"  実効時間平均 : {r['g_effective']:.4f} g  (= 0.5 × τ_A, τ_A = {r['mode_A']['duty']:.4f})")
    print(f"  基準速度: {r['base_rpm']:.3f} RPM  /  切替周期: {r['cycle_min']:.1f} min  "
          f"/  位相ドリフト周期: {r['drift_min']:.1f} min")

    for key, color in (("mode_A", "green"), ("mode_B", "magenta")):
        m = r[key]
        tag = "A" if key == "mode_A" else "B"
        print(c(f"\n  [Mode {tag}] {m['name']}  →  知覚重力 {m['g_perceived']:.1f} g",
                "bold", color, use_color=use_color))
        out_str = c(f"{m['omega_out_rpm']:9.5f} RPM", "bold", use_color=use_color)
        in_str = c(f"{m['omega_in_rpm']:9.5f} RPM", "bold", use_color=use_color)
        print(f"    時間配分  : {m['duty']*100:6.2f}%  ({m['duration_min']:.3f} min / cycle)")
        print(f"    ω_out     : {out_str}  ({m['omega_out_rad_s']:.6f} rad/s)")
        print(f"    ω_in      : {in_str}  ({m['omega_in_rad_s']:.6f} rad/s)")
        if key == "mode_A":
            print(f"    ε         : {m['epsilon_rad_s']:.4e} rad/s")

    print(c(f"\n  [診断] 細菌想定 v_sed={V_SED_BACTERIA*1e6:.1f} μm/s, "
            f"チャンバ半径 {r['chamber_r_mm']:.1f} mm",
            "dim", use_color=use_color))
    d = r["diagnostics"]
    safety = d["chamber_safety_ratio"]
    print(f"    ドリフト振幅: {d['drift_amplitude_mm']:.3f} mm  ({safety*100:.1f}% of chamber)  "
          + (c("[!] チャンバ越え", "red", use_color=use_color) if safety > 1
             else c("[OK]", "dim", use_color=use_color)))

    for w in r["warnings"]:
        print("  " + c("[!] " + w, "red", use_color=use_color))


def render_centrifuge(r: dict, label: str = "", use_color: bool = True) -> None:
    head = f"目標重力: {r['g_target']:.4g} g  ({r['g_target'] * G_EARTH:.4f} m/s²)"
    if label:
        head = f"[{label}] " + head
    print(c(head + "   [Centrifuge モード]", "bold", "yellow", use_color=use_color))
    print(f"  実効時間平均 : {r['g_effective']:.4f} g  (= ω_in² × r_sample / g0)")
    print(f"  サンプル位置 : 内軸から r = {r['sample_r_cm']:.1f} cm (例: 内側フレーム4角)")

    print(c("\n  [回転設定]", "bold", "green", use_color=use_color))
    in_str = c(f"{r['omega_in_rpm']:9.4f} RPM", "bold", use_color=use_color)
    out_str = c(f"{r['omega_out_rpm']:9.4f} RPM", "bold", use_color=use_color)
    print(f"    内側 ω_in  : {in_str}  ({r['omega_in_rad_s']:.4f} rad/s)  [遠心力で g_target 負荷]")
    print(f"    外側 ω_out : {out_str}  ({r['omega_out_rad_s']:.4f} rad/s)  [低速クリノローテーション]")

    d = r["diagnostics"]
    print(c(f"\n  [診断] 細菌想定 v_sed={V_SED_BACTERIA*1e6:.1f} μm/s, "
            f"チャンバ半径 {r['chamber_r_mm']:.1f} mm, 培養時間 {r['culture_hours']:.0f} h",
            "dim", use_color=use_color))
    print(f"    外側位置 (r={R_OUTER_FRAME*100:.1f}cm) での遠心外乱: "
          f"{d['outer_disturb_g']:.4f} g  "
          + (c("[!] 大", "red", use_color=use_color) if d['outer_disturb_g'] > r['g_target']*0.1
             else c("[OK]", "dim", use_color=use_color)))
    print(f"    細胞沈降速度    : {d['sed_velocity_um_s']:.3f} μm/s "
          f"(= v_sed × g_target)")

    twr = d["t_reach_wall_h"]
    if twr == float("inf"):
        wall_str = "壁到達なし (g_target=0)"
    elif twr < 1:
        wall_str = f"{twr*60:.1f} 分後"
    else:
        wall_str = f"{twr:.2f} 時間後"
    print(f"    細胞が壁に到達  : {wall_str}  "
          + (c("[!] 速い", "red", use_color=use_color) if twr < 1
             else c("", "dim", use_color=use_color)))

    if d["wall_dwell_ratio"] > 0:
        dwell = c(f"{d['wall_dwell_ratio']*100:.1f}%", "red" if d['wall_dwell_ratio'] > 0.5 else "yellow",
                  use_color=use_color)
        print(f"    培養時間中の壁滞在: {d['wall_dwell_h']:.2f} h ({dwell})  "
              "← 火星/月でも同じ現象、Mars analog 実験では正常な挙動")

    for w in r["warnings"]:
        print("  " + c("[!] " + w, "red", use_color=use_color))


def render_hybrid(r: dict, label: str = "", use_color: bool = True) -> None:
    head = f"目標重力: {r['g_target']:.4g} g  ({r['g_target'] * G_EARTH:.4f} m/s²)"
    if label:
        head = f"[{label}] " + head
    print(c(head + "   [Hybrid モード]", "bold", "yellow", use_color=use_color))
    print(f"  実効時間平均 : {r['g_effective']:.4f} g  "
          f"(= √(g_tilt² + g_cent²))")
    print(f"  サンプル位置 : 軸交点から r = {r['sample_r_cm']:.2f} cm")
    print(f"  位相ドリフト周期: {r['drift_min']:.1f} min")
    print(f"  内訳         : g_tilt = {r['g_tilt_component']:.4g} g  "
          f"+ g_cent = {r['g_centrifugal_component']:.4e} g (直交合成)")

    print(c("\n  [ハードウェア設定]", "bold", "blue", use_color=use_color))
    psi_str = c(f"{r['psi_deg']:8.4f}°", "bold", use_color=use_color)
    print(f"    外軸傾斜角 ψ : {psi_str}  ({r['psi_rad']:.6f} rad)  [鉛直から]")

    print(c("\n  [回転設定]", "bold", "green", use_color=use_color))
    out_str = c(f"{r['omega_out_rpm']:9.5f} RPM", "bold", use_color=use_color)
    in_str = c(f"{r['omega_in_rpm']:9.5f} RPM", "bold", use_color=use_color)
    print(f"    外側 ω_out : {out_str}  ({r['omega_out_rad_s']:.6f} rad/s)  [傾斜軸周り高速]")
    print(f"    内側 ω_in  : {in_str}  ({r['omega_in_rad_s']:.6e} rad/s)  [方向ドリフト用低速]")

    d = r["diagnostics"]
    print(c(f"\n  [診断] 細菌想定 v_sed={V_SED_BACTERIA*1e6:.1f} μm/s, "
            f"チャンバ半径 {r['chamber_r_mm']:.1f} mm, 培養時間 {r['culture_hours']:.0f} h",
            "dim", use_color=use_color))
    print(f"    tilt ドリフト振幅 (閉軌道, 振動)  : {d['tilt_drift_amp_mm']:.4f} mm")
    print(f"    遠心由来累積変位 (body fixed, 蓄積): {d['centrifugal_displacement_mm']:.4f} mm")
    print(f"    合計変位                          : {d['total_displacement_mm']:.4f} mm  "
          f"({d['chamber_safety_ratio']*100:.2f}% of chamber)  "
          + (c("[!] チャンバ越え", "red", use_color=use_color)
             if d['chamber_safety_ratio'] > 1
             else c("[OK]", "dim", use_color=use_color)))

    for w in r["warnings"]:
        print("  " + c("[!] " + w, "red", use_color=use_color))


def render(result: dict, label: str = "", use_color: bool = True) -> None:
    mode = result["mode"]
    if mode == "tilt":
        render_tilt(result, label=label, use_color=use_color)
    elif mode == "centrifuge":
        render_centrifuge(result, label=label, use_color=use_color)
    elif mode == "hybrid":
        render_hybrid(result, label=label, use_color=use_color)
    else:
        render_switching(result, label=label, use_color=use_color)


def compute(g_target, mode, **kw):
    if mode == "tilt":
        return compute_tilt(g_target, base_rpm=kw["base_rpm"],
                            drift_min=kw["drift_min"], chamber_r=kw["chamber_r"])
    if mode == "switching":
        return compute_switching(g_target, base_rpm=kw["base_rpm"],
                                 drift_min=kw["drift_min"], cycle_min=kw["cycle_min"],
                                 chamber_r=kw["chamber_r"])
    if mode == "centrifuge":
        return compute_centrifuge(g_target, sample_r=kw["sample_r"],
                                  outer_rpm=kw["base_rpm"], chamber_r=kw["chamber_r"],
                                  culture_hours=kw["culture_hours"])
    if mode == "hybrid":
        return compute_hybrid(g_target, sample_r=kw["sample_r"],
                              base_rpm=kw["base_rpm"], drift_min=kw["drift_min"],
                              chamber_r=kw["chamber_r"], culture_hours=kw["culture_hours"])
    raise ValueError(f"unknown mode: {mode}")


def parse_input(text: str) -> float:
    s = text.strip().lower()
    if s in PRESETS:
        return PRESETS[s][1]
    for suffix in ("m/s^2", "m/s2", "ms2", "ms-2"):
        if s.endswith(suffix):
            return float(s[: -len(suffix)].strip()) / G_EARTH
    return float(s)


# ---------------- 対話モード -----------------

def repl(mode: str, use_color: bool = True, **kw) -> None:
    print_header(use_color=use_color)
    print(f"モード: {mode}    数値 (g単位) または プリセット名 を入力.")
    print("  例:  0.378   |  mars   |  3.72 m/s2")
    print("  プリセット: " + ", ".join(PRESETS.keys()))
    print(c("  終了するには  quit / exit / Ctrl-D", "dim", use_color=use_color))
    print()
    while True:
        try:
            raw = input(c("g> ", "bold", "cyan", use_color=use_color))
        except (EOFError, KeyboardInterrupt):
            print()
            break
        s = raw.strip()
        if not s:
            continue
        if s.lower() in ("q", "quit", "exit"):
            break
        if s.lower() in ("h", "help", "?"):
            print("  プリセット: " + ", ".join(PRESETS.keys()))
            continue
        if s.lower() in ("list", "presets"):
            for k, (jp, g) in PRESETS.items():
                print(f"  {k:<7} ({jp:<8}) : {g:g} g")
            continue
        try:
            g = parse_input(s)
            res = compute(g, mode, **kw)
        except ValueError as e:
            print(c(f"  エラー: {e}", "red", use_color=use_color))
            continue
        label = PRESETS[s.lower()][0] if s.lower() in PRESETS else ""
        render(res, label=label, use_color=use_color)
        print()


# ---------------- CLI -----------------

def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="2軸クリノスタット 回転設定計算 (内側r=8.5cm/外側r=14.5cm)",
    )
    p.add_argument("gravity", nargs="?",
                   help="目標重力(g単位 or プリセット or 'X m/s2'). 省略時は対話モード.")
    p.add_argument("--mode",
                   choices=["tilt", "centrifuge", "hybrid", "switching", "all"],
                   default="tilt",
                   help="計算モード (既定: tilt). "
                        "tilt=サンプル=軸交点, 外軸傾斜で 0..1g 沈殿なし; "
                        "centrifuge=サンプル=4角 (r=8.5cm), 内側遠心で g_target、沈殿あり; "
                        "hybrid=任意の sample_r で tilt+遠心の合成、小オフセットを許容; "
                        "switching=軸交点+1:1/黄金比切替で 0..0.5g 傾斜なし; "
                        "all=全モード並べて表示")
    p.add_argument("--all", action="store_true", help="代表プリセットを一括表示")
    p.add_argument("--base-rpm", type=float, default=DEFAULT_BASE_RPM,
                   help=f"基準回転速度 [RPM] (既定: {DEFAULT_BASE_RPM}; "
                        "centrifuge モードでは外側速度として使用)")
    p.add_argument("--drift-min", type=float, default=DEFAULT_DRIFT_MIN,
                   help=f"位相ドリフト周期 [min] (既定: {DEFAULT_DRIFT_MIN}, "
                        "tilt/switching モードのみ)")
    p.add_argument("--cycle-min", type=float, default=DEFAULT_CYCLE_MIN,
                   help=f"switching モード時の切替周期 [min] (既定: {DEFAULT_CYCLE_MIN})")
    p.add_argument("--chamber-mm", type=float, default=R_CHAMBER_DEFAULT * 1000,
                   help=f"チャンバ半径 [mm] (既定: {R_CHAMBER_DEFAULT*1000:.1f})")
    p.add_argument("--sample-r-cm", type=float, default=R_INNER_FRAME * 100,
                   help=f"centrifuge/hybrid モード: サンプル位置の内軸からの半径 [cm] "
                        f"(centrifuge 既定: {R_INNER_FRAME*100:.1f} = 4角, "
                        f"hybrid では 0〜10 cm 程度を想定)")
    p.add_argument("--culture-hours", type=float, default=24.0,
                   help="centrifuge モードの沈降診断で使う培養時間 [h] (既定: 24)")
    p.add_argument("--no-color", action="store_true", help="ANSIカラー無効")
    args = p.parse_args(argv)

    use_color = not args.no_color and sys.stdout.isatty()
    kw = dict(base_rpm=args.base_rpm, drift_min=args.drift_min,
              cycle_min=args.cycle_min, chamber_r=args.chamber_mm / 1000.0,
              sample_r=args.sample_r_cm / 100.0, culture_hours=args.culture_hours)

    def emit(g, label=""):
        if args.mode == "all":
            for m in ("tilt", "centrifuge", "hybrid", "switching"):
                try:
                    render(compute(g, m, **kw), label=label, use_color=use_color)
                    print()
                except ValueError as e:
                    print(c(f"[{m}] {e}", "red", use_color=use_color))
                    print()
            return
        render(compute(g, args.mode, **kw), label=label, use_color=use_color)
        print()

    if args.all:
        print_header(use_color=use_color)
        for k, (jp, g) in PRESETS.items():
            try:
                emit(g, label=f"{k} ({jp})")
            except ValueError as e:
                print(c(f"[{k}] {e}", "red", use_color=use_color))
        return 0

    if args.gravity is None:
        if args.mode == "all":
            print(c("対話モードでは --mode tilt / centrifuge / switching のいずれかを選択してください",
                    "red", use_color=use_color), file=sys.stderr)
            return 1
        repl(mode=args.mode, use_color=use_color, **kw)
        return 0

    try:
        g = parse_input(args.gravity)
    except ValueError as e:
        print(c(f"エラー: {e}", "red", use_color=use_color), file=sys.stderr)
        return 1

    print_header(use_color=use_color)
    label = PRESETS[args.gravity.lower()][0] if args.gravity.lower() in PRESETS else ""
    try:
        emit(g, label=label)
    except ValueError as e:
        print(c(f"エラー: {e}", "red", use_color=use_color), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
