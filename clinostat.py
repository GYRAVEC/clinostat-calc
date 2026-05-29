#!/usr/bin/env python3
"""
2軸クリノスタット 回転設定計算ツール (Mode A + Mode B 混合方式)

  ハードウェア: 内側フレーム r = 8.5 cm / 外側フレーム r = 14.5 cm
  配置        : サンプルは 2軸の交点 (r→0) に置き、遠心力は使わない
  原理        : 1:1+位相ドリフト(Mode A; 0.5g) と 黄金比(Mode B; 0g) の
                時間配分で 0 ≤ g_target ≤ 0.5g を実現
"""

import argparse
import math
import sys

# 🐱 物理定数とフレーム諸元 (フレーム半径は機械的諸元として保持)
G_EARTH = 9.81                 # 地球重力加速度 [m/s²]
R_INNER_FRAME = 0.085          # 内側フレーム半径 [m] (8.5 cm)
R_OUTER_FRAME = 0.145          # 外側フレーム半径 [m] (14.5 cm)
PHI = (1 + math.sqrt(5)) / 2   # 黄金比

# 🐱 動作パラメータ既定値
DEFAULT_BASE_RPM = 5.0         # 基準角速度 [RPM]
DEFAULT_DRIFT_MIN = 20.0       # Mode A 位相ドリフト周期 [min]
DEFAULT_CYCLE_MIN = 60.0       # モード切替サイクル [min]

# 🐱 細菌培養想定の参照値（沈降診断に使用）
V_SED_BACTERIA = 1.0e-6        # 自由沈降速度 [m/s] (~1 μm/s)
R_CHAMBER_DEFAULT = 0.005      # チャンバ半径仮定 [m] (5 mm)
T_BIO_BACTERIA_S = 600.0       # 生物応答時間 [s] (~10 min)

# 🐱 達成可能上限 (1:1+ドリフトモードの時間平均値)
G_MAX_2AXIS = 0.5

# 🐱 プリセット (1g/earth は等速2軸では出せないので除外)
PRESETS = {
    "mars":   ("火星",     0.378),
    "moon":   ("月",       0.166),
    "ceres":  ("ケレス",   0.029),
    "europa": ("エウロパ", 0.134),
    "titan":  ("タイタン", 0.138),
    "iss":    ("ISS",      1.0e-4),
    "ug":     ("微小重力", 1.0e-3),
}


# ---------------- 物理計算 -----------------

def compute_profile(g_target: float,
                    base_rpm: float = DEFAULT_BASE_RPM,
                    drift_min: float = DEFAULT_DRIFT_MIN,
                    cycle_min: float = DEFAULT_CYCLE_MIN,
                    chamber_r: float = R_CHAMBER_DEFAULT) -> dict:
    """目標重力 → 2軸モーター設定値を返す.

    Mode A: ω_in = ω_out + ε  (1:1+位相ドリフト) → 知覚重力 0.5g
    Mode B: ω_in = ω_out / φ  (黄金比)            → 知覚重力 0g
    τ_A = 2·g_target で混合し、長時間平均 g_target を達成する.
    """
    if g_target < 0:
        raise ValueError("g_target は非負")
    if g_target > G_MAX_2AXIS + 1e-9:
        raise ValueError(
            f"g_target {g_target:g}g > {G_MAX_2AXIS}g は等速2軸+サンプル中心では不可。"
            " 上限は 0.5g (1:1モード時の時間平均)。"
        )

    omega_base = base_rpm * 2 * math.pi / 60.0           # rad/s
    T_drift = drift_min * 60.0                            # s
    T_cycle = cycle_min * 60.0                            # s

    # 🐱 Mode A: 1:1 + 微小オフセット ε で位相が T_drift で1周する
    eps = 2 * math.pi / T_drift                           # rad/s
    omega_out_A = omega_base
    omega_in_A = omega_base + eps

    # 🐱 Mode B: 黄金比で時間平均 0g
    omega_out_B = omega_base
    omega_in_B = omega_base / PHI

    # 🐱 時間配分
    duty_A = min(2.0 * g_target, 1.0)
    duty_B = 1.0 - duty_A
    t_A = duty_A * T_cycle
    t_B = duty_B * T_cycle

    # 🐱 細菌沈降診断
    #   Mode A 中の位相ドリフト中、細胞は v_sed * 0.5 で振幅運動
    #   1ドリフト周期中の最大振幅 = v_sed * 0.5 * T_drift / π (sin/cos の積分)
    drift_amp = V_SED_BACTERIA * 0.5 * T_drift / math.pi   # m
    chamber_ratio = drift_amp / chamber_r if chamber_r > 0 else float("inf")

    # 🐱 妥当性チェック
    warnings = []
    if T_drift < T_BIO_BACTERIA_S:
        warnings.append(
            f"T_drift ({drift_min:.1f}min) < T_bio (~{T_BIO_BACTERIA_S/60:.0f}min): "
            "細胞は 0.5g を知覚できず ~0g として応答する可能性"
        )
    if drift_amp > chamber_r:
        warnings.append(
            f"ドリフト振幅 {drift_amp*1000:.2f} mm > チャンバ半径 {chamber_r*1000:.1f} mm: "
            "T_drift を短く / 基準速度を上げる必要"
        )
    if omega_base <= 0:
        warnings.append("base_rpm が 0 以下")

    return {
        "g_target": g_target,
        "g_effective": 0.5 * duty_A,
        "base_rpm": base_rpm,
        "drift_min": drift_min,
        "cycle_min": cycle_min,
        "chamber_r_mm": chamber_r * 1000.0,
        "mode_A": {
            "name": "1:1+位相ドリフト",
            "duty": duty_A,
            "duration_min": t_A / 60.0,
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
            "duration_min": t_B / 60.0,
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
}


def c(text: str, *styles: str, use_color: bool = True) -> str:
    if not use_color:
        return text
    return "".join(ANSI[s] for s in styles) + text + ANSI["reset"]


def print_header(use_color: bool = True) -> None:
    line = "=" * 68
    print(c(line, "cyan", use_color=use_color))
    print(c("  2軸クリノスタット 部分重力プロファイル計算 (Mode A+B 方式)",
            "bold", "cyan", use_color=use_color))
    print(c(f"  内側 r = {R_INNER_FRAME*100:.1f} cm  /  外側 r = {R_OUTER_FRAME*100:.1f} cm  "
            f"/  サンプル位置 = 2軸交点", "dim", use_color=use_color))
    print(c(f"  g0 = {G_EARTH} m/s²  /  φ = {PHI:.6f}  /  達成可能 g_target ≤ {G_MAX_2AXIS}",
            "dim", use_color=use_color))
    print(c(line, "cyan", use_color=use_color))


def render_profile(r: dict, label: str = "", use_color: bool = True) -> None:
    head = f"目標重力: {r['g_target']:.4g} g  ({r['g_target'] * G_EARTH:.4f} m/s²)"
    if label:
        head = f"[{label}] " + head
    print(c(head, "bold", "yellow", use_color=use_color))
    print(f"  実効時間平均 : {r['g_effective']:.4f} g   "
          f"(= 0.5 × τ_A, τ_A = {r['mode_A']['duty']:.4f})")
    print(f"  基準速度     : {r['base_rpm']:.3f} RPM   "
          f"切替周期: {r['cycle_min']:.1f} min   位相ドリフト周期: {r['drift_min']:.1f} min")

    for key, color in (("mode_A", "green"), ("mode_B", "magenta")):
        m = r[key]
        tag = "A" if key == "mode_A" else "B"
        print(c(f"\n  [Mode {tag}] {m['name']}  →  知覚重力 {m['g_perceived']:.1f} g",
                "bold", color, use_color=use_color))
        out_rpm = m["omega_out_rpm"]
        in_rpm = m["omega_in_rpm"]
        out_rpm_str = c(f"{out_rpm:9.5f} RPM", "bold", use_color=use_color)
        in_rpm_str = c(f"{in_rpm:9.5f} RPM", "bold", use_color=use_color)
        print(f"    時間配分  : {m['duty']*100:6.2f}%  "
              f"({m['duration_min']:.3f} min / cycle)")
        print(f"    ω_out     : {out_rpm_str}  ({m['omega_out_rad_s']:.6f} rad/s)")
        print(f"    ω_in      : {in_rpm_str}  ({m['omega_in_rad_s']:.6f} rad/s)")
        if key == "mode_A":
            print(f"    ε (= ω_in - ω_out): {m['epsilon_rad_s']:.4e} rad/s")

    d = r["diagnostics"]
    print(c(f"\n  [診断] 細菌想定 v_sed={V_SED_BACTERIA*1e6:.1f} μm/s, "
            f"チャンバ半径 {r['chamber_r_mm']:.1f} mm",
            "dim", use_color=use_color))
    safety = d["chamber_safety_ratio"]
    safety_str = f"  ドリフト振幅: {d['drift_amplitude_mm']:.3f} mm  ({safety*100:.1f}% of chamber)"
    print(safety_str + ("  " + c("[!] チャンバ越え", "red", use_color=use_color)
                        if safety > 1 else "  " + c("[OK]", "dim", use_color=use_color)))

    for w in r["warnings"]:
        print("  " + c("[!] " + w, "red", use_color=use_color))


def parse_input(text: str) -> float:
    s = text.strip().lower()
    if s in PRESETS:
        return PRESETS[s][1]
    for suffix in ("m/s^2", "m/s2", "ms2", "ms-2"):
        if s.endswith(suffix):
            return float(s[: -len(suffix)].strip()) / G_EARTH
    return float(s)


# ---------------- 対話モード -----------------

def repl(use_color: bool = True, **kw) -> None:
    print_header(use_color=use_color)
    print("数値 (g単位) または プリセット名 を入力してください.")
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
            res = compute_profile(g, **kw)
        except ValueError as e:
            print(c(f"  エラー: {e}", "red", use_color=use_color))
            continue
        label = PRESETS[s.lower()][0] if s.lower() in PRESETS else ""
        render_profile(res, label=label, use_color=use_color)
        print()


# ---------------- CLI -----------------

def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="2軸クリノスタット 部分重力プロファイル計算 "
                    "(内側r=8.5cm/外側r=14.5cm, サンプル中心配置)",
    )
    p.add_argument("gravity", nargs="?",
                   help="目標重力(g単位 or プリセット名 or 'X m/s2'). 省略時は対話モード.")
    p.add_argument("--all", action="store_true", help="代表プリセットを一括表示")
    p.add_argument("--base-rpm", type=float, default=DEFAULT_BASE_RPM,
                   help=f"基準回転速度 [RPM] (既定: {DEFAULT_BASE_RPM})")
    p.add_argument("--drift-min", type=float, default=DEFAULT_DRIFT_MIN,
                   help=f"Mode A 位相ドリフト周期 [min] (既定: {DEFAULT_DRIFT_MIN})")
    p.add_argument("--cycle-min", type=float, default=DEFAULT_CYCLE_MIN,
                   help=f"モード切替サイクル [min] (既定: {DEFAULT_CYCLE_MIN})")
    p.add_argument("--chamber-mm", type=float, default=R_CHAMBER_DEFAULT * 1000,
                   help=f"チャンバ半径 [mm] (既定: {R_CHAMBER_DEFAULT*1000:.1f}, 沈降診断用)")
    p.add_argument("--no-color", action="store_true", help="ANSIカラー無効")
    args = p.parse_args(argv)

    use_color = not args.no_color and sys.stdout.isatty()
    kw = dict(base_rpm=args.base_rpm, drift_min=args.drift_min,
              cycle_min=args.cycle_min, chamber_r=args.chamber_mm / 1000.0)

    if args.all:
        print_header(use_color=use_color)
        for k, (jp, g) in PRESETS.items():
            try:
                render_profile(compute_profile(g, **kw), label=f"{k} ({jp})", use_color=use_color)
                print()
            except ValueError as e:
                print(c(f"[{k}] {e}", "red", use_color=use_color))
        return 0

    if args.gravity is None:
        repl(use_color=use_color, **kw)
        return 0

    try:
        g = parse_input(args.gravity)
        res = compute_profile(g, **kw)
    except ValueError as e:
        print(c(f"エラー: {e}", "red", use_color=use_color), file=sys.stderr)
        return 1

    print_header(use_color=use_color)
    label = PRESETS[args.gravity.lower()][0] if args.gravity.lower() in PRESETS else ""
    render_profile(res, label=label, use_color=use_color)
    return 0


if __name__ == "__main__":
    sys.exit(main())
