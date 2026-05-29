#!/usr/bin/env python3
"""
2軸クリノスタット 回転設定計算ツール

  ハードウェア: 内側フレーム r = 8.5 cm / 外側フレーム r = 14.5 cm
  配置        : サンプルは 2軸の交点 (r→0) に置き、遠心力は使わない

  Tilt モード (既定, 0 ≤ g_target ≤ 1g):
    外軸を鉛直から ψ = arccos(g_target) 傾け、ω_out 高速 + ω_in 低速で
    細胞は |g| = g0·cos ψ を短時間平均として知覚し、長時間平均で蓄積ゼロ.

  Switching モード (--mode switching, 0 ≤ g_target ≤ 0.5g):
    傾斜なし運用. Mode A (1:1+位相ドリフト 0.5g) と Mode B (黄金比 0g) を
    時間配分で混合する従来方式.
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
    print(c(f"  内側 r = {R_INNER_FRAME*100:.1f} cm  /  外側 r = {R_OUTER_FRAME*100:.1f} cm  "
            f"/  サンプル位置 = 2軸交点", "dim", use_color=use_color))
    print(c(f"  g0 = {G_EARTH} m/s²  /  φ = {PHI:.6f}",
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


def render(result: dict, label: str = "", use_color: bool = True) -> None:
    if result["mode"] == "tilt":
        render_tilt(result, label=label, use_color=use_color)
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
    p.add_argument("--mode", choices=["tilt", "switching", "both"], default="tilt",
                   help="計算モード (既定: tilt). "
                        "tilt=外軸傾斜で 0..1g 全域, "
                        "switching=Mode A+B 混合で 0..0.5g 傾斜なし, "
                        "both=両方を並べて表示")
    p.add_argument("--all", action="store_true", help="代表プリセットを一括表示")
    p.add_argument("--base-rpm", type=float, default=DEFAULT_BASE_RPM,
                   help=f"基準回転速度 [RPM] (既定: {DEFAULT_BASE_RPM})")
    p.add_argument("--drift-min", type=float, default=DEFAULT_DRIFT_MIN,
                   help=f"位相ドリフト周期 [min] (既定: {DEFAULT_DRIFT_MIN})")
    p.add_argument("--cycle-min", type=float, default=DEFAULT_CYCLE_MIN,
                   help=f"switching モード時の切替周期 [min] (既定: {DEFAULT_CYCLE_MIN})")
    p.add_argument("--chamber-mm", type=float, default=R_CHAMBER_DEFAULT * 1000,
                   help=f"チャンバ半径 [mm] (既定: {R_CHAMBER_DEFAULT*1000:.1f})")
    p.add_argument("--no-color", action="store_true", help="ANSIカラー無効")
    args = p.parse_args(argv)

    use_color = not args.no_color and sys.stdout.isatty()
    kw = dict(base_rpm=args.base_rpm, drift_min=args.drift_min,
              cycle_min=args.cycle_min, chamber_r=args.chamber_mm / 1000.0)

    def emit(g, label=""):
        if args.mode == "both":
            try:
                render(compute(g, "tilt", **kw), label=label, use_color=use_color)
                print()
            except ValueError as e:
                print(c(f"[tilt] {e}", "red", use_color=use_color))
            try:
                render(compute(g, "switching", **kw), label=label, use_color=use_color)
                print()
            except ValueError as e:
                print(c(f"[switching] {e}", "red", use_color=use_color))
                print()
        else:
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
        if args.mode == "both":
            print(c("対話モードでは --mode tilt または switching を選択してください",
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
