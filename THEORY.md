# Tilt モード理論的妥当性検証レポート

**対象**: clinostat-calc (2軸クリノスタット, 内側 r = 8.5 cm / 外側 r = 14.5 cm, サンプル位置 = 2軸交点)
**モード**: Tilt モード (外軸を鉛直から ψ 傾けて回転)
**結論**: サンプル系での **短時間平均重力** が `|⟨g⟩| = g₀·cos ψ` で目標値に一致し、**長時間平均はゼロ** で蓄積方向もないことを示す。よって、入力 g_target に対し ψ = arccos(g_target) と運用パラメータを与えれば、生物的に "g_target を知覚 + 細胞は蓄積しない" 状態を理論的に実現する。

---

## 1. 設計仕様と前提

| 項目 | 記号 | 値・条件 |
|------|------|---------|
| 地球重力加速度 | g₀ | 9.81 m/s² |
| 外軸の鉛直からの傾斜 | ψ | 0 ≤ ψ ≤ π/2 |
| 外軸まわり回転角 | α(t) = ω_out·t | ω_out: 高速 |
| 内軸まわり回転角 | β(t) = ω_in·t | ω_in: 超低速 |
| サンプル位置 | r⃗_sample = 0⃗ | 2軸の交点 (遠心力ゼロ) |
| 細胞自由沈降速度 | v_sed | 〜 1 μm/s (バクテリア) |
| 細胞生物応答時間 | T_bio | 〜 10 min |
| チャンバ半径 | R_ch | 〜 5 mm |

サンプルが 2軸の交点にあるため、**遠心力(ω²r) はゼロ**。働く唯一の慣性力は地球重力で、これを2軸回転で時間平均することによって "実効重力" を作る。

---

## 2. 座標系と回転の定義

### 慣性 (Lab) 系

x̂_lab, ŷ_lab, ẑ_lab の右手直交基底。地球重力は

```
g⃗_lab = -g₀ · ẑ_lab
```

### 回転軸の取り方

外軸方向 (Lab 系で固定):

```
n̂_out = sin ψ · x̂_lab + cos ψ · ẑ_lab
```

ψ = 0 で鉛直 (n̂_out = ẑ_lab)、ψ = π/2 で水平 (n̂_out = x̂_lab)。

内軸は外軸に垂直な、外フレーム内で固定された方向。Lab 系基底で表すなら t=0 時点で:

```
n̂_in = cos ψ · x̂_lab - sin ψ · ẑ_lab    ( ≡ p̂ )
```

(n̂_out · p̂ = sin ψ cos ψ - sin ψ cos ψ = 0 で直交 ✓)

3本目の基底:

```
q̂ = n̂_out × p̂ = ŷ_lab
```

{p̂, q̂, n̂_out} は右手直交基底 (これを以後 "傾斜基底" と呼ぶ)。

### 回転作用素

外回転 (Lab → 外フレーム): n̂_out のまわり -α の受動回転 R_out(α)
内回転 (外フレーム → サンプル): p̂ のまわり -β の受動回転 R_in(β)

サンプル系で見える重力ベクトル:

```
g⃗_sample(t) = R_in(-β) · R_out(-α) · g⃗_lab
```

---

## 3. サンプル系での重力ベクトル導出

### 3.1 Lab → 外フレーム

ẑ_lab を傾斜基底で分解:

```
ẑ_lab = cos ψ · n̂_out − sin ψ · p̂
```

R_out(-α) は n̂_out を不変、(p̂, q̂) 平面内で角度 α 回転:

```
R_out(-α) · p̂ = cos α · p̂ − sin α · q̂
```

よって外フレーム内での重力:

```
g⃗_outer = -g₀ · R_out(-α) · ẑ_lab
        = -g₀ cos ψ · n̂_out + g₀ sin ψ cos α · p̂ − g₀ sin ψ sin α · q̂
```

### 3.2 外フレーム → サンプル

R_in(-β) は p̂ を不変、(q̂, n̂_out) 平面内で角度 β 回転:

```
R_in(-β) · q̂      = cos β · q̂ − sin β · n̂_out
R_in(-β) · n̂_out  = sin β · q̂ + cos β · n̂_out
```

これを g⃗_outer に作用:

```
g⃗_sample(α, β) = -g₀ cos ψ (sin β q̂ + cos β n̂_out)
                + g₀ sin ψ cos α · p̂
                - g₀ sin ψ sin α (cos β q̂ - sin β n̂_out)
```

### 3.3 成分

傾斜基底 {p̂, q̂, n̂_out} で:

```
g_p   =  g₀ sin ψ cos α
g_q   = -g₀ cos ψ sin β - g₀ sin ψ sin α cos β
g_n   = -g₀ cos ψ cos β + g₀ sin ψ sin α sin β
```

これは時刻依存の閉じた式。以下ではこれの時間平均を取る。

---

## 4. 時間平均と知覚重力

3つのタイムスケールを定義:

```
T_out = 2π / ω_out   ≪   T_bio   ≪   T_drift = 2π / ω_in   ≪   T_exp
```

例: 既定値 (ω_out = 5 RPM, ω_in = 0.05 RPM) なら
T_out = 12 s, T_bio ≈ 600 s, T_drift = 1200 s, T_exp = 数時間。

### 4.1 短時間平均 (T_bio スケール, 細胞が知覚する量)

T_bio の間 α は多数回転 (ω_out T_bio = π·5·600/30 ≈ 314 rad → 50回転以上) するので

```
⟨cos α⟩_T_bio ≈ 0,   ⟨sin α⟩_T_bio ≈ 0
```

一方 β は ω_in T_bio = 0.05·2π/60·600 ≈ 0.31 rad で **ほぼ定数** (= β₀)。

成分代入:

```
⟨g_p⟩_T_bio    ≈  g₀ sin ψ · 0 = 0
⟨g_q⟩_T_bio    ≈ -g₀ cos ψ · sin β₀  - g₀ sin ψ · 0 · cos β₀ = -g₀ cos ψ sin β₀
⟨g_n⟩_T_bio    ≈ -g₀ cos ψ · cos β₀  + g₀ sin ψ · 0 · sin β₀ = -g₀ cos ψ cos β₀
```

よって:

```
⟨g⃗⟩_T_bio = -g₀ cos ψ · ( sin β₀ q̂ + cos β₀ n̂_out )
```

**ベクトルの大きさ**:

```
|⟨g⃗⟩_T_bio| = g₀ cos ψ · √(sin²β₀ + cos²β₀) = g₀ cos ψ
```

→ **β₀ の値によらず常に g₀·cos ψ**。これが g_target そのもの。

**方向**: (q̂, n̂_out) 平面内、初期値 -n̂_out 方向から角度 β₀ 回った位置。β₀ は ω_in t で変化するので、**ゆっくり (T_drift = 2π/ω_in 周期で) 1周する**。

### 4.2 長時間平均 (T_exp >> T_drift)

α と β の両方が独立に多数回転するエルゴード平均:

```
⟨cos α⟩ = ⟨sin α⟩ = ⟨cos β⟩ = ⟨sin β⟩ = 0
⟨cross terms⟩ = ⟨..⟩ · ⟨..⟩ = 0
```

各成分:

```
⟨g_p⟩∞ = g₀ sin ψ · 0 = 0
⟨g_q⟩∞ = -g₀ cos ψ · 0 - g₀ sin ψ · 0 = 0
⟨g_n⟩∞ = -g₀ cos ψ · 0 + g₀ sin ψ · 0 = 0
```

→ **⟨g⃗⟩∞ = 0⃗** (長時間平均はゼロ)。

### 4.3 物理的解釈

| タイムスケール | ⟨g⃗⟩         | 物理的意味                                |
|--------------|---------------|------------------------------------------|
| T_out (秒)   | 振動する瞬間値 | 外回転で水平成分が急速に均される          |
| T_bio (分)   | g₀ cos ψ 一定 | **細胞はこの大きさの重力を知覚**          |
| T_drift (時) | 0 へ向かう    | 方向が1周することで蓄積が消える            |

これが「**短時間で部分重力、長時間でゼロ → 部分重力を感じながら蓄積しない**」の理論的構造。

---

## 5. 沈殿挙動の解析

### 5.1 細胞の即時速度

ストークス則の終端速度仮定で、細胞は重力方向に v_sed × |g⃗|/g₀ で運動。サンプル系での細胞の運動方程式:

```
dx⃗/dt = (v_sed / g₀) · g⃗_sample(t)
```

### 5.2 T_bio スケールの変位 (細胞応答時間内)

⟨g⃗⟩_T_bio = -g₀ cos ψ · (sin β₀ q̂ + cos β₀ n̂_out) を使うと:

```
Δx⃗(T_bio) ≈ -v_sed cos ψ · T_bio · (sin β₀ q̂ + cos β₀ n̂_out)
|Δx⃗(T_bio)| = v_sed · g_target · T_bio       (g_target = cos ψ 単位なし)
```

数値 (火星 g_target=0.378, v_sed=1 μm/s, T_bio=600s):
|Δx⃗| = 1 × 0.378 × 600 = **227 μm** ← 細胞が T_bio で漂う距離 (生物的信号)

### 5.3 T_drift スケールの正味変位 (蓄積判定)

β = ω_in t が 0 → 2π を一周する間の積分:

```
Δx⃗(T_drift) = ∫₀^(T_drift) (v_sed/g₀) · g⃗_sample dt
            ≈ (v_sed/g₀) · ∫₀^(T_drift) ⟨g⃗⟩_T_bio (β₀=ω_in t) dt
            = -v_sed cos ψ · ∫₀^(T_drift) (sin(ω_in t) q̂ + cos(ω_in t) n̂_out) dt
            = -v_sed cos ψ · ( [-cos(ω_in t)/ω_in]₀^(T_drift) · q̂ + [sin(ω_in t)/ω_in]₀^(T_drift) · n̂_out )
            = -v_sed cos ψ · ( 0 · q̂ + 0 · n̂_out )
            = 0⃗
```

→ **1ドリフト周期で正味変位ゼロ**、細胞は出発点に戻る。

### 5.4 周期内の最大振幅 (チャンバ越え判定)

β が 0 → 2π を一周する間の `|Δx⃗|` の最大値:

```
Δx⃗(t) = -v_sed cos ψ · ( [-cos(ω_in t) + 1]/ω_in · q̂ + sin(ω_in t)/ω_in · n̂_out )
|Δx⃗(t)|² = (v_sed cos ψ / ω_in)² · ((1 - cos(ω_in t))² + sin²(ω_in t))
        = (v_sed cos ψ / ω_in)² · (2 - 2 cos(ω_in t))
最大は ω_in t = π で、(2 - 2(-1)) = 4 ⇒ |Δx⃗|_max = 2 v_sed cos ψ / ω_in
```

ω_in = 2π/T_drift を代入:

```
|Δx⃗|_max = v_sed · g_target · T_drift / π
```

数値 (火星, T_drift = 1200s):
|Δx⃗|_max = 1 × 0.378 × 1200 / π = **144 μm** ← clinostat-calc の "ドリフト振幅" 出力と一致 ✓

チャンバ半径 5 mm に対して 144 μm = 2.9%、十分小さい。

### 5.5 沈殿が成立しない条件

|Δx⃗|_max < R_ch を満たせば、細胞はチャンバ内で円形軌道を描いて壁に到達しない:

```
v_sed · g_target · T_drift / π  <  R_ch
T_drift  <  π · R_ch / (v_sed · g_target)
```

バクテリア・5 mm チャンバ・火星: T_drift < π · 5000 / (1 · 0.378) = 41,500 s ≈ 11.5 時間。
T_drift = 20 min はこれよりはるかに小さい → 蓄積なし ✓

---

## 6. 数値検証 (clinostat-calc 出力との突合)

ツールの `compute_tilt()` が返す値と、本レポートの解析式を比較:

| 量                         | 解析式                              | ツール出力 (火星) | 一致 |
|----------------------------|------------------------------------|------------------|------|
| ψ                          | arccos(g_target)                    | 67.7901°         | ✓   |
| 実効時間平均 \|⟨g⟩\|_T_bio | g₀·cos ψ                            | 0.3780 g         | ✓   |
| ω_out                      | 設定値                              | 5.00000 RPM      | ✓   |
| ω_in                       | 2π / T_drift                        | 0.05000 RPM      | ✓   |
| ドリフト振幅 \|Δx\|_max    | v_sed · g_target · T_drift / π     | 0.144 mm         | ✓   |

`compute_tilt()` の実装 ([clinostat.py](clinostat.py)) は:

```python
psi_rad = math.acos(min(g_target, 1.0))             # ψ = arccos(g_target)
omega_out = base_rpm * 2 * math.pi / 60.0           # ユーザ指定
omega_in  = 2 * math.pi / T_drift                   # T_drift から逆算
drift_amp = V_SED_BACTERIA * g_target * T_drift / math.pi   # |Δx|_max
```

→ 解析式と完全一致。実装は理論を忠実に反映している。

### 各 g_target に対する妥当性

| g_target | ψ      | 実効重力 = g₀·cos ψ | 沈殿振幅 (T_drift=20min) | 適性     |
|----------|--------|---------------------|--------------------------|----------|
| 1.0g     | 0°    | 9.810 m/s²          | 0.382 mm                 | 静置と等価 (回転不要) |
| 0.8g     | 36.87° | 7.848 m/s²          | 0.306 mm                 | OK       |
| 0.5g     | 60°   | 4.905 m/s²          | 0.191 mm                 | OK       |
| 0.378g (Mars) | 67.79° | 3.708 m/s²    | 0.144 mm                 | OK       |
| 0.166g (Moon) | 80.44° | 1.629 m/s²    | 0.063 mm                 | OK       |
| 0.01g    | 89.43° | 0.098 m/s²          | 0.004 mm                 | OK       |

全範囲で **|⟨g⟩| = g_target × g₀ 厳密に成立**, 沈殿振幅は線形に g_target に比例しチャンバ内に収まる。

---

## 7. 仮定と適用限界

### 7.1 物理的仮定

1. **サンプル位置 = 2軸交点** ─ 偏心していると遠心力 ω²r が混入し、g_target がずれる
2. **剛体回転** ─ 流体の慣性 (角加速度時の流体スロッシング) は無視
3. **Coriolis 力無視** ─ 細胞速度 v_sed ≪ 回転速度 ωR で十分成立
4. **線形抵抗 (Stokes)** ─ Re が大きい大型粒子では破綻
5. **均一回転** ─ ω_out, ω_in 一定。脈動が大きいと短時間平均が時間変動する

### 7.2 時間スケール条件

理論が成立するには:

```
ω_out · T_bio  ≫  2π   ⇔   T_bio  ≫  T_out
ω_in  · T_bio  ≪  2π   ⇔   T_bio  ≪  T_drift
T_drift · v_sed · g_target  <  π · R_ch
```

既定値はこれを満たしている。ω_out を下げすぎたり ω_in を上げすぎると破綻 → スクリプトは警告を出す。

### 7.3 適用限界

- **g_target > 1.0g 不可** (ψ が虚数になる; 地球重力を超えるには遠心モード等の追加機構が必要)
- **g_target ≈ 1.0g** で ψ = 0 となり外回転は意味を持たない (傾斜なし静置でよい)
- **真空中・低粘度流体** では v_sed が増大し T_drift 上限が下がる
- **大粒子サンプル** (細胞凝集塊, 種子等) では Stokes 則が破綻し別解析が必要

---

## 8. 結論

| 項目 | 結論 |
|------|------|
| **大きさ** | 短時間平均 \|⟨g⃗⟩\|_T_bio が **厳密に g₀·cos ψ = g_target·g₀** に一致 (§4.1) |
| **沈殿防止** | T_drift 周期で位置がループし、正味変位ゼロ。最大振幅 \|Δx\|_max がチャンバ半径より十分小なら蓄積なし (§5.3-4) |
| **遠心力混入** | サンプルが 2軸交点にあるため遠心力寄与は厳密にゼロ (§1) |
| **実装** | clinostat-calc の compute_tilt() は解析式と一致 (§6) |
| **適用範囲** | g_target ∈ [0, 1] · g₀, 細菌スケール (v_sed ≈ 1 μm/s) で T_drift 〜 20 min 設定下、チャンバ半径 5 mm 内で全範囲蓄積なし (§6) |

→ Tilt モードは **理論的に正しく目標重力を負荷し、かつ細菌培養スケールで沈殿問題を回避する**。

clinostat-calc が `python3 clinostat.py <g_target>` で返す (ψ, ω_out, ω_in) を実機に設定すれば、サンプル系で目標 g_target を知覚させつつ蓄積方向のない gravity field を実現できる。

---

## 付録 A. 解析の独立検算

特殊極限値でのチェック:

- **ψ = 0**: cos ψ = 1, 全成分 ⟨g_n⟩ → -g₀ × 1 = -g₀. 1g down ✓ (回転しなくても同じ)
- **ψ = π/2**: cos ψ = 0, ⟨g⃗⟩_T_bio = 0⃗. 微小重力 ✓
- **ψ = π/3 (0.5g)**: cos ψ = 0.5, |⟨g⟩| = 0.5 g₀ ✓
- **ω_in → 0 (内回転なし)**: β₀ 固定で ⟨g⃗⟩_T_bio = -g₀ cos ψ n̂_out (傾斜軸方向に静的重力) → 長時間平均も g₀ cos ψ (蓄積する) → 沈殿防止には ω_in > 0 必要 ✓
- **ω_in → ∞ (内も速い)**: β も時間平均されゼロ → ⟨g⃗⟩_T_bio = 0 → 微小重力になり目標を達成しない → 適切な速度分離が必要 ✓

これらは §4.1 の結果と整合し、解析の妥当性を独立に裏付ける。

## 付録 B. 参考文献的位置付け

- 1軸クリノローテーション: van Loon, J. J. (2007) "Some history and use of the random positioning machine, RPM, in gravity related research", *Advances in Space Research*
- 傾斜軸クリノスタット (slow tilted rotation): Hoson, T. et al. (1992) "Evaluation of the three-dimensional clinostat as a simulator of weightlessness", *Planta*
- 部分重力シミュレーション: Borst & van Loon (2009) "Technology and developments for the Random Positioning Machine"

本ツールの Tilt モードは「内軸を超低速ドリフトに使う傾斜クリノスタット」の変種で、上記文献群が扱う partial gravity simulator の理論的枠組みに収まる。
