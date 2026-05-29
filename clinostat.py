#!/usr/bin/env python3
"""
2軸クリノスタット 回転速度計算ツール
内側フレーム半径: 8.5 cm / 外側フレーム半径: 14.5 cm
"""

import argparse
import math
import sys

# 🐱 物理定数とフレーム諸元
G_EARTH = 9.81                 # 地球重力加速度 [m/s²]
R_INNER = 0.085                # 内側フレーム回転半径 [m]
R_OUTER = 0.145                # 外側フレーム回転半径 [m]
PHI = (1 + math.sqrt(5)) / 2   # 黄金比（非共鳴比）

# 🐱 代表的な天体重力プリセット (g0 単位)
PRESETS = {
    "earth":   ("地球",       1.000),
    "mars":    ("火星",       0.378),
    "moon":    ("月",         0.166),
    "ceres":   ("ケレス",     0.029),
    "europa":  ("エウロパ",   0.134),
    "titan":   ("タイタン",   0.138),
    "iss":     ("ISS (微小)", 1e-4),
    "ug":      ("微小重力",   1e-3),
}


def clinostat_speeds(g_target: float) -> dict:
    """目標重力を入力に内側・外側の回転速度を計算する."""
    if g_target <= 0:
        raise ValueError("g_target は正の値を指定してください")

    # 🐱 内側: 遠心力で目標g生成  ω_in = sqrt(g_target * g0 / r_in)
    omega_in = math.sqrt(g_target * G_EARTH / R_INNER)
    rpm_in = omega_in * 60.0 / (2 * math.pi)

    # 🐱 外側: 黄金比で非共鳴クリノローテーション
    omega_out = omega_in / PHI
    rpm_out = omega_out * 60.0 / (2 * math.pi)

    a_in = omega_in ** 2 * R_INNER
    a_out = omega_out ** 2 * R_OUTER

    return {
        "g_target": g_target,
        "inner": {"rpm": rpm_in, "omega": omega_in, "a_ms2": a_in, "a_g": a_in / G_EARTH},
        "outer": {"rpm": rpm_out, "omega": omega_out, "a_ms2": a_out, "a_g": a_out / G_EARTH},
    }


# ---------------- 表示系 -----------------

ANSI = {
    "reset":  "\033[0m",
    "bold":   "\033[1m",
    "dim":    "\033[2m",
    "cyan":   "\033[36m",
    "yellow": "\033[33m",
    "green":  "\033[32m",
    "magenta":"\033[35m",
    "red":    "\033[31m",
}


def c(text: str, *styles: str, use_color: bool = True) -> str:
    if not use_color:
        return text
    return "".join(ANSI[s] for s in styles) + text + ANSI["reset"]


def print_header(use_color: bool = True) -> None:
    line = "=" * 62
    print(c(line, "cyan", use_color=use_color))
    print(c("  2軸クリノスタット 回転速度計算ツール", "bold", "cyan", use_color=use_color))
    print(c(f"  内側 r = {R_INNER*100:.1f} cm  /  外側 r = {R_OUTER*100:.1f} cm  /  g0 = {G_EARTH} m/s²",
            "dim", use_color=use_color))
    print(c(line, "cyan", use_color=use_color))


def render(result: dict, label: str = "", use_color: bool = True) -> None:
    g = result["g_target"]
    head = f"目標重力: {g:.4g} g  ({g * G_EARTH:.4f} m/s²)"
    if label:
        head = f"[{label}] " + head
    print(c(head, "bold", "yellow", use_color=use_color))

    inner = result["inner"]
    outer = result["outer"]

    inner_rpm_str = c(f"{inner['rpm']:9.4f} RPM", "bold", use_color=use_color)
    outer_rpm_str = c(f"{outer['rpm']:9.4f} RPM", "bold", use_color=use_color)

    print(c(f"  内側フレーム (r = {R_INNER*100:.1f} cm)", "bold", "green", use_color=use_color))
    print(f"    回転数      : {inner_rpm_str}")
    print(f"    角速度      : {inner['omega']:.6f} rad/s")
    print(f"    遠心加速度  : {inner['a_ms2']:.4f} m/s²  ({inner['a_g']:.4f} g)")

    print(c(f"  外側フレーム (r = {R_OUTER*100:.1f} cm)", "bold", "magenta", use_color=use_color))
    print(f"    回転数      : {outer_rpm_str}")
    print(f"    角速度      : {outer['omega']:.6f} rad/s")
    print(f"    遠心加速度  : {outer['a_ms2']:.4f} m/s²  ({outer['a_g']:.4f} g)  ", end="")
    # 🐱 外側の遠心成分が目標を超えていたら警告
    if outer["a_g"] > g * 0.5 and g > 0.01:
        print(c("[!] 外乱大", "red", use_color=use_color))
    else:
        print(c("[OK]", "dim", use_color=use_color))


def parse_input(text: str) -> float:
    """ユーザー入力を g 単位に変換する.

    対応:
      - 'mars', 'moon' などのプリセット
      - '0.38'  (g 単位)
      - '3.72m/s2' / '3.72 ms2' / '3.72 m/s^2'
    """
    s = text.strip().lower()
    if s in PRESETS:
        return PRESETS[s][1]
    # m/s² 単位
    for suffix in ("m/s^2", "m/s2", "ms2", "ms-2"):
        if s.endswith(suffix):
            val = float(s[: -len(suffix)].strip())
            return val / G_EARTH
    return float(s)


# ---------------- 対話モード -----------------

def repl(use_color: bool = True) -> None:
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
                print(f"  {k:<8} ({jp:<10}) : {g:g} g")
            continue

        try:
            g_target = parse_input(s)
        except ValueError:
            print(c(f"  入力を解釈できませんでした: {s!r}", "red", use_color=use_color))
            continue

        try:
            res = clinostat_speeds(g_target)
        except ValueError as e:
            print(c(f"  エラー: {e}", "red", use_color=use_color))
            continue

        # 🐱 プリセット名なら表示に反映
        label = ""
        if s.lower() in PRESETS:
            label = PRESETS[s.lower()][0]
        render(res, label=label, use_color=use_color)
        print()


# ---------------- CLI エントリ -----------------

def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="2軸クリノスタット 回転速度計算 (内側r=8.5cm / 外側r=14.5cm)",
    )
    p.add_argument("gravity", nargs="?",
                   help="目標重力(g単位). 例: 0.378  /  mars  /  '3.72 m/s2'. 省略時は対話モード.")
    p.add_argument("--all", action="store_true", help="代表的プリセットを全て表示")
    p.add_argument("--no-color", action="store_true", help="ANSIカラーを無効化")
    args = p.parse_args(argv)
    use_color = not args.no_color and sys.stdout.isatty()

    if args.all:
        print_header(use_color=use_color)
        for k, (jp, g) in PRESETS.items():
            render(clinostat_speeds(g), label=f"{k} ({jp})", use_color=use_color)
            print()
        return 0

    if args.gravity is None:
        repl(use_color=use_color)
        return 0

    try:
        g_target = parse_input(args.gravity)
        res = clinostat_speeds(g_target)
    except ValueError as e:
        print(c(f"エラー: {e}", "red", use_color=use_color), file=sys.stderr)
        return 1

    print_header(use_color=use_color)
    label = PRESETS[args.gravity.lower()][0] if args.gravity.lower() in PRESETS else ""
    render(res, label=label, use_color=use_color)
    return 0


if __name__ == "__main__":
    sys.exit(main())
