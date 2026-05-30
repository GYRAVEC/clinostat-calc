# 時系列重力シミュレーション開発ガイド

サンプル位置・フレームサイズ・回転設定をもとに、**サンプル系で時刻 t に細胞が感じる重力ベクトル g_sample(t)** を時間離散化でシミュレーションする方法。地球重力は無視せず正確に変換に組み込む。本ドキュメントは別PCで実装するときの設計図。

---

## 1. ゴール

入力:
- モード (tilt / centrifuge / hybrid / switching)
- 目標重力 g_target (例: 0.378g)
- サンプル位置 (内軸からの半径 r [m]、初期姿勢)
- フレーム諸元 (r_inner = 8.5 cm, r_outer = 14.5 cm)
- 回転設定 (ψ, ω_in, ω_out, 必要なら時系列)
- シミュレーション時間と分解能

出力:
- 時系列 CSV: `t [s], gx_lab, gy_lab, gz_lab, gx_sample, gy_sample, gz_sample, |g_sample|`
- プロット (4枚): 3成分時系列 / |g| 時系列 / 単位球上トレース / 短時間平均ヒストグラム
- 標準出力: 解析的に求めた理論平均値との突合

---

## 2. 必要なライブラリ

```python
import numpy as np                  # ベクトル・行列演算
import matplotlib.pyplot as plt     # プロット
from scipy.spatial.transform import Rotation as R   # 回転行列の生成
import argparse, json, csv
```

それだけ。numpy + matplotlib + scipy.spatial.transform で十分。Cython や Numba は1秒分解能 × 数千秒なら不要。

---

## 3. 物理モデル

### 3.1 座標系と回転の合成

3つの座標系:
- 実験室系 (lab): 床固定、重力 = `-g₀ ẑ_lab`
- 外フレーム系 (outer): 外軸まわりに `α(t) = ω_out × t` で回転
- サンプル系 (body): 外フレーム上の内軸まわりに `β(t) = ω_in × t` で回転

サンプル系 → 実験室系 への変換 (アクティブ回転):

```python
def rotation_lab_from_body(alpha, beta, psi, outer_axis="x"):
    """
    psi: 外軸を鉛直 (lab z) から倒した角度 [rad]
    外軸の lab 系での向き: n̂_out = (sin ψ, 0, cos ψ)
    内軸: 外フレーム内で n̂_out に直交した固定軸
    """
    # 外軸をベクトルで表現
    n_out = np.array([np.sin(psi), 0.0, np.cos(psi)])
    # 内軸 (外フレーム静止時の lab x 軸 を ψ 回した方向、外軸と直交)
    n_in_outer_frame = np.array([np.cos(psi), 0.0, -np.sin(psi)])

    # まず外フレーム空間で内軸まわりに β 回転
    R_inner = R.from_rotvec(beta * n_in_outer_frame)
    # 次に外軸まわりに α 回転
    R_outer = R.from_rotvec(alpha * n_out)
    return R_outer * R_inner   # 合成 (順序重要: 外が後で適用)
```

### 3.2 サンプル系での重力ベクトル

実験室系での重力は単純:

```python
g_lab = np.array([0.0, 0.0, -G_EARTH])   # (0, 0, -9.81)
```

これをサンプル系に変換するには、上の合成回転の **逆** を適用 (lab 系の固定ベクトルをサンプル系で見たもの):

```python
R_total = rotation_lab_from_body(alpha, beta, psi)
g_sample_from_earth = R_total.inv().apply(g_lab)
```

### 3.3 遠心力 (サンプルが r ≠ 0 にいる場合)

サンプル位置を外フレーム系で `p_outer = (0, r, 0)` のように指定 (内軸が外フレームの x 軸とすると、内軸に直交方向に r 離れた位置)。

実験室系での位置:

```python
p_lab_at_t = R_outer.apply(R_inner.apply(p_outer))
```

そこに働く遠心力は **2つ** (内外それぞれの軸まわりから):

**内側軸まわりの遠心力 (サンプル系で固定方向)**:
```python
# サンプル系での position は p_inner = R_inner.inv().apply(p_outer) ... 
# でも内軸はサンプルが回らないとそのまま p_outer なので simpler:
# r_inner_axis_to_sample (in body frame) = p_outer (rotated by R_inner)
# 簡単化: サンプル位置を body 系で (0, 0, r) とすれば
g_cent_inner_body = omega_in**2 * r * np.array([0, 0, 1])  # body 固定方向
g_cent_inner_lab = R_total.apply(g_cent_inner_body)
```

**外側軸まわりの遠心力 (時間変動)**:
```python
# サンプルの外軸からの垂直距離ベクトル (lab 系)
n_out_lab = np.array([np.sin(psi), 0, np.cos(psi)])
p_perp = p_lab_at_t - np.dot(p_lab_at_t, n_out_lab) * n_out_lab
R_perp = np.linalg.norm(p_perp)
g_cent_outer_lab = omega_out**2 * R_perp * (p_perp / R_perp)  # 外向き
```

### 3.4 合計の実効重力

```python
g_eff_lab = g_lab + g_cent_inner_lab + g_cent_outer_lab
g_eff_sample = R_total.inv().apply(g_eff_lab)
```

Coriolis 力は細胞の自発速度が分からないと計算できない (Stokes 終端速度仮定なら無視可)。生物応答だけ見たいなら省略でOK。

---

## 4. メインループ

```python
def simulate(mode, g_target, sample_r, psi, omega_in, omega_out,
             duration_s, dt_s, p_body_init):
    """
    mode: 'tilt' | 'centrifuge' | 'hybrid' | 'switching'
    duration_s: シミュレーション総時間 [秒]
    dt_s: 時間分解能 [秒]  (推奨: 1/(10*max(omega)) より小)
    p_body_init: サンプル位置 (body系) [m] 例: np.array([0, 0, 0.02])
    """
    n_steps = int(duration_s / dt_s)
    ts = np.arange(n_steps) * dt_s

    records = np.zeros((n_steps, 7))  # t, g_x_sample, g_y_sample, g_z_sample, |g|, alpha, beta

    for i, t in enumerate(ts):
        alpha = omega_out * t
        beta = omega_in * t

        # 必要なら switching モードでここで ω を時刻ごとに変える
        if mode == "switching":
            cycle_pos = (t / 3600.0) % 1.0   # 60min サイクル想定 (要パラメータ化)
            if cycle_pos < duty_A:
                # Mode A: 1:1+ドリフト
                omega_in_t = omega_out_base + 2*np.pi/T_drift
            else:
                # Mode B: 黄金比
                omega_in_t = omega_out_base / PHI
            beta = ... # 時刻積分が必要、別管理

        R_total = rotation_lab_from_body(alpha, beta, psi)
        g_lab = np.array([0, 0, -G_EARTH])

        # 遠心 (sample_r > 0 のとき)
        if sample_r > 1e-9:
            g_cent_inner_body = omega_in**2 * sample_r * np.array([0, 0, 1])
            g_cent_inner_lab = R_total.apply(g_cent_inner_body)
            p_lab = R_total.apply(p_body_init)
            n_out_lab = np.array([np.sin(psi), 0, np.cos(psi)])
            p_perp = p_lab - np.dot(p_lab, n_out_lab) * n_out_lab
            R_perp = np.linalg.norm(p_perp)
            if R_perp > 1e-12:
                g_cent_outer_lab = omega_out**2 * R_perp * (p_perp / R_perp)
            else:
                g_cent_outer_lab = np.zeros(3)
        else:
            g_cent_inner_lab = np.zeros(3)
            g_cent_outer_lab = np.zeros(3)

        g_total_lab = g_lab + g_cent_inner_lab + g_cent_outer_lab
        g_sample = R_total.inv().apply(g_total_lab)

        records[i] = [t, g_sample[0], g_sample[1], g_sample[2],
                      np.linalg.norm(g_sample), alpha, beta]

    return ts, records
```

時間分解能の目安: tilt モードで ω_out = 5 RPM (周期 12秒) なら `dt = 0.1秒` で十分。centrifuge モードで ω_in = 63 RPM (周期 0.95秒) なら `dt = 0.01秒`。

---

## 5. プロット (matplotlib 4枚)

```python
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# (1) 3成分の時系列
axes[0,0].plot(ts, records[:,1], label='g_x')
axes[0,0].plot(ts, records[:,2], label='g_y')
axes[0,0].plot(ts, records[:,3], label='g_z')
axes[0,0].axhline(0, color='gray', lw=0.5)
axes[0,0].set_xlabel('time [s]'); axes[0,0].set_ylabel('g_sample component [m/s²]')
axes[0,0].legend(); axes[0,0].set_title('サンプル系での重力 3成分')

# (2) |g_sample| の時系列 + 理論値ライン
axes[0,1].plot(ts, records[:,4])
axes[0,1].axhline(g_target * G_EARTH, color='red', linestyle='--', label=f'g_target = {g_target}g')
axes[0,1].set_xlabel('time [s]'); axes[0,1].set_ylabel('|g_sample| [m/s²]')
axes[0,1].set_title('重力ベクトルの大きさ'); axes[0,1].legend()

# (3) 単位球上のトレース (3D)
ax3d = fig.add_subplot(2, 2, 3, projection='3d')
g_hat = records[:,1:4] / records[:,4:5]  # 単位ベクトル
ax3d.plot(g_hat[:,0], g_hat[:,1], g_hat[:,2], lw=0.3, alpha=0.5)
ax3d.set_title('単位球上の g 方向トレース')

# (4) 短時間平均 (移動平均) の時系列
window = int(600 / dt_s)   # T_bio = 600s 窓
moving_avg = np.zeros((n_steps, 3))
for i in range(3):
    moving_avg[:,i] = np.convolve(records[:,i+1], np.ones(window)/window, mode='same')
moving_avg_mag = np.linalg.norm(moving_avg, axis=1)
axes[1,1].plot(ts, moving_avg_mag, label='|⟨g⟩_T_bio|')
axes[1,1].axhline(g_target * G_EARTH, color='red', linestyle='--', label=f'g_target')
axes[1,1].set_xlabel('time [s]'); axes[1,1].set_ylabel('|⟨g⟩| [m/s²]')
axes[1,1].set_title(f'短時間平均 (T_bio = 600s 窓)'); axes[1,1].legend()

plt.tight_layout()
plt.savefig(f'simulation_{mode}_{g_target}.png', dpi=150)
```

---

## 6. 検算

各モードでシミュレーション結果が理論値と合うかを確認する診断を入れる:

```python
# 短時間平均 (T_bio = 600s 窓) を取って g_target と比較
T_bio_window_steps = int(600 / dt_s)
late_window = records[-T_bio_window_steps:, 1:4]  # 最後の 10分
mean_g_sample = np.mean(late_window, axis=0)
mag_mean = np.linalg.norm(mean_g_sample)

print(f"短時間平均 |⟨g⟩|: {mag_mean:.4f} m/s² ({mag_mean/G_EARTH:.4f} g)")
print(f"理論値 g_target × g0: {g_target * G_EARTH:.4f} m/s² ({g_target:.4f} g)")
print(f"誤差: {abs(mag_mean - g_target*G_EARTH)/G_EARTH:.4e} g")

# 長時間平均 (全期間)
long_avg = np.mean(records[:, 1:4], axis=0)
print(f"長時間平均 ⟨g⟩: ({long_avg[0]:.4e}, {long_avg[1]:.4e}, {long_avg[2]:.4e})")
print(f"  ↑ tilt モードなら全成分が長時間で 0 に近づくはず")
```

期待される結果:
| モード | 短時間平均 \|⟨g⟩\| | 長時間平均 ⟨g⟩ |
|--------|--------------------|------------------|
| tilt   | g_target × g0      | 0 (各成分)        |
| centrifuge | g_target × g0  | g_target × g0 (body固定方向) |
| hybrid | g_target × g0      | 微小 (≈ g_cent × time) |
| switching | g_target × g0 (時間平均) | 0 |

---

## 7. 細胞軌道シミュレーション (オプション、レベル4)

g_sample(t) が出たら、細胞の運動方程式 (低 Re 仮定、Stokes drag balance) で位置を積分:

```python
# 細胞は重力の方向に Stokes 終端速度 v_sed で動く
def cell_trajectory(records, ts, v_sed=1e-6, n_cells=100, chamber_r=0.005):
    # 初期位置: チャンバ中心付近に N 個ランダム配置
    positions = np.random.uniform(-chamber_r*0.1, chamber_r*0.1, (n_cells, 3))
    history = np.zeros((len(ts), n_cells, 3))

    for i, t in enumerate(ts[:-1]):
        dt = ts[i+1] - ts[i]
        g_vec = records[i, 1:4]
        g_unit = g_vec / np.linalg.norm(g_vec)
        # 沈降速度 = v_sed × (|g|/g0) の方向に
        v = v_sed * np.linalg.norm(g_vec) / G_EARTH * g_unit
        positions += v * dt
        # チャンバ壁で止める
        radial = np.linalg.norm(positions, axis=1, keepdims=True)
        mask = (radial > chamber_r).flatten()
        positions[mask] = positions[mask] / radial[mask] * chamber_r
        history[i] = positions

    return history

# ヒストグラム表示
final_positions = history[-1]
plt.hist(np.linalg.norm(final_positions, axis=1), bins=30)
plt.xlabel('center からの距離 [m]'); plt.ylabel('細胞数')
```

ブラウン運動を入れたければ毎ステップ `+ np.random.normal(0, np.sqrt(2*D*dt), 3)` を足す。

---

## 8. CLI 構造 (案)

```python
parser = argparse.ArgumentParser()
parser.add_argument('--mode', choices=['tilt', 'centrifuge', 'hybrid', 'switching'])
parser.add_argument('--g-target', type=float, required=True)
parser.add_argument('--sample-r-cm', type=float, default=0.0)
parser.add_argument('--duration', type=float, default=600.0, help='[秒]')
parser.add_argument('--dt', type=float, default=0.1, help='[秒]')
parser.add_argument('--base-rpm', type=float, default=5.0)
parser.add_argument('--drift-min', type=float, default=20.0)
parser.add_argument('--output-prefix', default='simulation')
parser.add_argument('--cells', type=int, default=0, help='細胞数 (>0 で軌道計算)')
parser.add_argument('--plot', action='store_true')
parser.add_argument('--csv', action='store_true')
args = parser.parse_args()

# モードごとに psi, omega_in, omega_out を計算 (clinostat-calc の関数を import)
from clinostat import compute_tilt, compute_centrifuge, compute_hybrid, compute_switching
...
```

---

## 9. 想定する実装手順 (1日コース)

1. **30分**: numpy + scipy.spatial.transform.Rotation を `pip install` し、3.1 の `rotation_lab_from_body` をテスト
   - 検算: psi=0, alpha=π/2 で `g_lab = (0,0,-9.81)` がサンプル系でも同じになるか
2. **1時間**: 3.4 の合計重力計算と、短いシミュレーション (300秒分) で記録
3. **30分**: 4枚プロット書く
4. **30分**: 検算出力 (短期/長期平均) と理論値の突合
5. **1時間**: CLI 化と各モードの分岐対応
6. **オプション 2時間**: 細胞軌道シミュレーションを追加

合計 3-4時間で動くものができる。

---

## 10. 検証用テストケース

実装したら最初にこれを通す:

| ケース | 入力 | 期待出力 |
|--------|------|---------|
| 静置 | mode=tilt, psi=0, omega=0 | g_sample = (0,0,-9.81) 一定 |
| 横倒し | mode=tilt, psi=π/2, omega=0 | g_sample = (g₀ sin α, 0, -g₀ cos α) で振動 |
| 微小重力 | mode=tilt, psi=π/2, omega_out=5RPM, omega_in=3.09RPM (golden ratio) | 長時間平均で ⟨g⟩ → 0 |
| 火星 tilt | mode=tilt, g_target=0.378 | 短期 \|⟨g⟩\| ≈ 0.378×g₀, 長期 → 0 |
| 火星 centrifuge | mode=centrifuge, r=8.5cm, ω_in=63RPM | \|g_sample\| が定数 ω²r + 振動成分 |

---

## 11. 参考: 既存実装

このディレクトリの [clinostat.py](clinostat.py) には、各モードの ψ, ω_in, ω_out を計算する `compute_tilt()` などが既にある。simulate.py からはそれを import して使うのが速い:

```python
from clinostat import (
    compute_tilt, compute_centrifuge, compute_hybrid, compute_switching,
    G_EARTH, R_INNER_FRAME, R_OUTER_FRAME, PHI, V_SED_BACTERIA
)
```

理論的妥当性 (短時間/長時間平均がそれぞれ何になるか) は [THEORY.md](THEORY.md) §6-§8 (tilt), §12 (centrifuge), §13 (hybrid) に書いてある。シミュレーション結果との突合の根拠として使える。

---

## 12. 余裕があれば

- **3D アニメーション** (matplotlib.animation.FuncAnimation): 装置の枠 + サンプル + 重力ベクトルを毎フレーム描画して GIF/MP4 出力
- **plotly** で対話的可視化: ブラウザで重力トレースをぐりぐり回せる
- **streamlit** で web GUI 化: パラメータをスライダーで動かしながら結果を見る
- **numba.jit** でメインループ高速化: 24時間 × 10ms 分解能 = 8640万ステップを数秒で
- **Coriolis 力** を加える: 細胞の運動方向に依存するので軌道計算とセットで
