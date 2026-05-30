# Claude プロンプト: 時系列重力シミュレーション実装依頼

> このファイルは「新規 Claude セッションに貼り付けて、最高品質のシミュレーションツール `simulate.py` を完成させるための仕様書兼依頼書」です。下の `---PROMPT START---` から `---PROMPT END---` までをそのまま貼ってください。

---PROMPT START---

# 依頼: clinostat-calc に最高品質の時系列重力シミュレーションを実装

あなたは 2軸クリノスタット研究用 CLI ツール `clinostat-calc` の開発を引き継ぎます。`simulate.py` を新規実装し、サンプル本体系で細胞が時刻 t に感じる重力ベクトル `g_sample(t)` を **物理的に厳密かつ数値的に高精度** にシミュレーションし、検証・可視化・データ出力までこなすツールを完成させてください。研究レポートや論文の図に直接使える成果物が目標です。

## まず読んでほしい資料 (順番に)

1. リポジトリ: https://github.com/GYRAVEC/clinostat-calc
2. `README.md` — プロジェクト全体像と4モードの位置づけ
3. `THEORY.md` — 物理導出 (§3〜§8 が tilt モード、§12 が centrifuge、§13 が hybrid、§14 が従来手法との比較)
4. `SIMULATION.md` — 既存の実装ガイド (大まかな設計案、これは下敷きにして良いが必須ではない)
5. `clinostat.py` — 各モードの計算式が定数と関数の形で実装済み。シミュレーションは **この関数を import して再利用** すること

これらを把握してから着手してください。理論的妥当性は `THEORY.md` の式に従い、シミュレーション結果が解析値と矛盾なく一致することが品質の前提条件です。

## 背景: 装置と4つの運用モード

2軸クリノスタットは:
- 内側フレーム半径: 8.5 cm
- 外側フレーム半径: 14.5 cm
- 用途: 細菌培養 (培養時間 24-48時間) で部分重力 (火星 0.378g など) を再現

4つのモード:

| モード | サンプル位置 | 重力源 | 沈殿 |
|--------|------------|--------|------|
| tilt | 2軸交点 (r≈0) | 地球重力射影 g₀·cos ψ | なし |
| centrifuge | 内側4角 (r=8.5cm) | 内側遠心 ω_in²·r | あり |
| hybrid | 任意 r | tilt と遠心の直交合成 | 主にtilt |
| switching | 2軸交点 | 1:1 同期 と 黄金比 の時間配分 | なし |

## 実装する成果物

### 必須ファイル

```
simulate.py                # メインのシミュレーション CLI
physics.py                 # 回転・力・座標変換の純物理関数
visualize.py               # 4枚プロット + 3Dアニメーション
tests/test_simulate.py     # pytest によるTC1〜TC8の自動テスト
requirements.txt           # 依存パッケージ
```

### 期待される CLI

```bash
# 基本: tilt モードで 30分シミュレーション、プロットと CSV を出力
python3 simulate.py --mode tilt --g-target 0.378 --duration 1800 --csv

# 細胞 200個の軌道と Brown 運動をシミュレーション、3Dアニメ出力
python3 simulate.py --mode centrifuge --g-target 0.378 \
    --duration 86400 --cells 200 --diffusion --animate

# 全モード検証テスト
python3 simulate.py --validate

# JSON で全パラメータ指定
python3 simulate.py --config experiment.json
```

### 全引数 (CLI)

| 引数 | 型 | 既定 | 説明 |
|------|----|----|-----|
| `--mode` | enum | (必須) | tilt/centrifuge/hybrid/switching/custom |
| `--g-target` | float | (必須) | 目標重力 (g 単位) |
| `--sample-r-cm` | float | 0.0 | サンプル位置の内軸からの半径 [cm] |
| `--duration` | float | 600 | シミュレーション時間 [秒] |
| `--dt` | float | (自動) | 時間ステップ [秒]、未指定なら ω 最大値の 1/20 |
| `--integrator` | enum | rk4 | 積分器 euler/rk4/symplectic |
| `--base-rpm` | float | 5.0 | 基準回転速度 |
| `--drift-min` | float | 20.0 | 位相ドリフト周期 [min] |
| `--cycle-min` | float | 60.0 | switching の切替周期 [min] |
| `--chamber-mm` | float | 5.0 | チャンバ半径 [mm] |
| `--cells` | int | 0 | 細胞数 (>0 で軌道シミュ) |
| `--diffusion` | flag | False | Brown 運動を含める |
| `--coriolis` | flag | False | Coriolis 力を含める |
| `--v-sed-um-s` | float | 1.0 | 細胞の自由沈降速度 [μm/s] |
| `--diffusion-D` | float | 1e-13 | 拡散係数 [m²/s] |
| `--animate` | flag | False | 3D アニメーション出力 |
| `--anim-fps` | int | 30 | アニメフレームレート |
| `--anim-format` | enum | mp4 | mp4/gif/webp |
| `--anim-speed` | float | 1.0 | アニメ再生速度倍率 |
| `--csv` | flag | False | CSV 出力 |
| `--no-plot` | flag | False | プロット抑制 |
| `--output-prefix` | str | "sim" | 出力ファイル接頭辞 |
| `--seed` | int | None | RNG シード (--cells や --diffusion 利用時) |
| `--validate` | flag | False | TC1〜TC8 全部を実行して結果を表で出力 |
| `--config` | str | None | JSON 設定ファイル |
| `--quiet` | flag | False | 進捗バー以外抑制 |
| `--verbose` | flag | False | 詳細ログ |

## 物理モデル仕様

### 座標系の定義

- 実験室系 (lab): 床固定。地球重力は `g_lab = (0, 0, -g₀)`、g₀ = 9.81 m/s²
- 外フレーム系 (outer): 外軸まわりに `α(t) = ω_out · t` で回転
- サンプル系 (body): 外フレーム上の内軸まわりに `β(t) = ω_in · t` で回転

### 回転の合成

`scipy.spatial.transform.Rotation` を使い、自前で行列を組まないこと:

```python
from scipy.spatial.transform import Rotation as R

n_out = np.array([np.sin(psi), 0.0, np.cos(psi)])      # 外軸 (lab 系)
n_in_outer = np.array([np.cos(psi), 0.0, -np.sin(psi)]) # 内軸 (外フレーム系で固定)

R_inner = R.from_rotvec(beta * n_in_outer)
R_outer = R.from_rotvec(alpha * n_out)
R_total = R_outer * R_inner                # body→lab の合成
```

### 力の各項

サンプル位置 `p_body` (body 系) における力は次の3項の合成:

**(1) 地球重力**
```python
g_lab = np.array([0, 0, -G_EARTH])
g_sample_from_earth = R_total.inv().apply(g_lab)
```

**(2) 内側遠心 (body 固定方向)**
```python
# サンプル位置が内軸から r 離れている場合の遠心力
r_inner_radial = np.linalg.norm(p_body[1:])  # 内軸 (x_body) に直交成分
# サンプル系で body 固定方向 (z_body) を仮定:
g_cent_inner_body = omega_in**2 * r_inner_radial * np.array([0, 0, 1])
```

**(3) 外側遠心 (時間変動)**
```python
p_lab = R_total.apply(p_body)
n_out_lab = np.array([np.sin(psi), 0, np.cos(psi)])
p_perp = p_lab - np.dot(p_lab, n_out_lab) * n_out_lab
R_perp_norm = np.linalg.norm(p_perp)
if R_perp_norm > 1e-12:
    F_cent_outer_lab = omega_out**2 * R_perp_norm * (p_perp / R_perp_norm)
else:
    F_cent_outer_lab = np.zeros(3)
g_cent_outer_sample = R_total.inv().apply(F_cent_outer_lab)
```

**合計**
```python
g_eff_sample = g_sample_from_earth + g_cent_inner_body + g_cent_outer_sample
```

### Coriolis 力 (`--coriolis` 指定時のみ)

細胞の自発速度 `v_cell_sample` (sample 系) に対して:

```
F_coriolis_sample = -2 m · ω_sample × v_cell_sample
```

`ω_sample` は両回転を sample 系で合成した角速度ベクトル。`R.from_rotvec` の `as_rotvec() / dt` で取れます。

### 細胞運動 (`--cells N` 指定時)

Stokes 抗力釣り合いで終端速度。質量小・粘性大の極限:

```
v_cell_sample = (v_sed / g₀) · g_eff_sample
```

ここに Brown 運動 (`--diffusion`) を加えるなら:

```python
v_cell_sample += np.random.normal(0, np.sqrt(2*D/dt), 3)
```

`D` は拡散係数 (`--diffusion-D` で指定、既定 1e-13 m²/s = 大腸菌)。Coriolis を加える場合は上式に `F_coriolis/μ` を追加。

### 数値積分

`--integrator` で選択:
- `euler`: 一次 (デフォルト dt 小さくする必要、検証用)
- `rk4`: 四次ルンゲクッタ (既定、汎用)
- `symplectic`: シンプレクティック (長時間の保存則を厳守したいとき)

時間ステップは `--dt` 未指定なら `1 / (20 × max(ω_in, ω_out)) [rad/s 単位]` で自動。

### 数値安定性

- 回転行列は scipy.spatial.transform 使用で自動で SO(3) に保たれる (浮動小数のドリフトなし)
- 細胞位置がチャンバ壁を越えたら反射 or 壁張り付き (`--reflect` で切替、既定: 張り付き)
- 加速度から速度に行く積分中は陽的に正規化テストを入れる

## 出力仕様

### CSV (`{prefix}_{mode}_{g}.csv`)

列 (`--csv` 指定時に出力):
```
time_s, alpha_rad, beta_rad,
g_sample_x, g_sample_y, g_sample_z, g_magnitude,
g_earth_x, g_earth_y, g_earth_z,
g_cent_inner_x, g_cent_inner_y, g_cent_inner_z,
g_cent_outer_x, g_cent_outer_y, g_cent_outer_z,
[--cells N 指定時: cell_i_x, cell_i_y, cell_i_z for i in 0..N-1]
```

`--duration` が長い場合は `--dt` を粗くするか、`--output-stride INT` で間引き保存できるようにしてください。

### プロット (PNG, 300 dpi)

`{prefix}_{mode}_{g}_overview.png`: 6パネル
- (a) g_sample 3成分の時系列 (時間軸は秒〜時間で自動スケール)
- (b) |g_sample| の時系列 + 理論値ライン (g_target × g₀)
- (c) 単位球上の g 方向トレース (3D scatter / line)
- (d) 短時間平均 (T_bio=600s 移動平均) の時系列
- (e) 周波数スペクトル (FFT) — 主要周波数の確認
- (f) g 成分の累積平均 (徐々に長時間平均値に収束する様子)

`{prefix}_{mode}_{g}_cells.png` (--cells > 0 のとき):
- (a) 細胞分布 2D ヒートマップ (xy 平面 + xz 平面)
- (b) 半径分布ヒストグラム (時刻 0, T/4, T/2, T で重ね描き)
- (c) 平均位置の時間発展
- (d) 個別細胞 5本の軌道線

すべてのプロット: タイトル/軸ラベル/凡例/グリッド/単位を完備。色は colorblind-safe (例: `tab10` パレット)。

### アニメーション (--animate 指定時)

`{prefix}_{mode}_{g}_animation.mp4`:
- 3D 立体に装置の外フレーム (グレーのワイヤー) + 内フレーム + サンプル位置 (青丸) + 細胞 (赤点) + 重力ベクトル (緑矢印、サンプル系)
- 左上に時刻表示、左下にパラメータ表示 (mode, g_target, ψ, ω 等)
- フレームレート `--anim-fps` (既定 30)、解像度 1920×1080
- 動画長さ = duration / anim-speed
- mp4 (ffmpeg 必須) または gif (Pillow) を `--anim-format` で選択
- 視点は左下から見上げる固定 (elev=20°, azim=-60°)

### サマリー JSON (`{prefix}_{mode}_{g}_summary.json`)

```json
{
  "mode": "tilt",
  "g_target": 0.378,
  "params": { /* 全 CLI 引数を再録 */ },
  "device_params": { /* psi, omega_in, omega_out, ... */ },
  "diagnostics": {
    "short_time_average_g_magnitude": 0.3778,
    "short_time_average_theoretical": 0.3780,
    "short_time_error_pct": 0.053,
    "long_time_average_g_vector": [1.2e-6, -3.4e-6, 8.9e-7],
    "long_time_average_magnitude": 4.0e-6,
    "drift_amplitude_simulated_mm": 0.1442,
    "drift_amplitude_theoretical_mm": 0.1444,
    "max_orthonormality_error": 1.2e-14,
    "max_energy_drift_rel": 3.4e-10
  },
  "tests_passed": ["TC1", "TC2", ...],
  "tests_failed": [],
  "runtime_s": 12.4,
  "n_steps": 60000,
  "integrator": "rk4"
}
```

## 必須検証テスト (TC1〜TC8)

`pytest tests/test_simulate.py` で全て pass すること。また `simulate.py --validate` で表形式で結果を出力:

| ID | ケース | 入力 | 期待出力 | 許容誤差 |
|----|--------|------|---------|---------|
| TC1 | 静置 | mode=tilt, psi=0, omega=0, t=10s | g_sample ≡ (0,0,-9.81) m/s² | < 1e-12 |
| TC2 | 水平 1軸回転 | mode=tilt, psi=π/2, ω_out=5 RPM, ω_in=0, t=60s | 周期 12s で振動、長期平均 \|⟨g⟩\| | < 0.01g |
| TC3 | 微小重力 RPM | tilt, ψ=π/2, ω_out=5RPM, ω_in=ω_out/φ, t=1h | 長時間平均 \|⟨g⟩\| → 0 | < 0.001g |
| TC4 | 火星 tilt | tilt, g_target=0.378, t=1h | 短期 \|⟨g⟩\|_T_bio = 0.378g、長期 \|⟨g⟩\| → 0 | 短期 1% / 長期 < 0.005g |
| TC5 | 火星 centrifuge | centrifuge, r=8.5cm, g=0.378 | サンプル系で 1成分が 3.708 m/s² 一定 | < 0.1% |
| TC6 | hybrid 2cm オフセット | hybrid, r=2cm, g=0.378 | tilt とほぼ同じ、誤差 | < 1e-5 g |
| TC7 | switching 月 | switching, g=0.166, t=2h | 60min サイクル内で 2モード切替、時間平均 0.166g | 1% |
| TC8 | 48時間長時間 | tilt, g=0.378, t=48h, integrator=rk4 | 直交性 (R^T R - I) ノルム < 1e-10 | 厳守 |

各テストの結果は標準出力に表で:

```
=== Validation Results ===
TC1: 静置                       [PASS]  error: 0.00e+00
TC2: 水平 1軸回転               [PASS]  error: 1.23e-04
TC3: 微小重力 RPM              [PASS]  error: 3.45e-04
...
=== Summary: 8/8 passed ===
```

## コード品質要件 (必須)

1. **型ヒント**: 全ての公開関数に Python 型注釈 (numpy.ndarray, Optional[X], TypedDict 等を活用)
2. **Docstring**: NumPy 形式で全公開関数に説明と引数・戻り値・例
3. **モジュール分割**: `physics.py` / `simulate.py` / `visualize.py` / `tests/` で責務分離
4. **エラー処理**: 不正入力には `ValueError` で明確なメッセージ (英語OK)
5. **再現性**: `--seed` 指定で乱数を含めて完全再現
6. **テスト**: `pytest -v tests/` で全 pass。カバレッジ 80% 以上目標
7. **依存**: `requirements.txt` に必要なパッケージとバージョン範囲を記載
   ```
   numpy>=1.24
   scipy>=1.10
   matplotlib>=3.6
   ffmpeg-python  # アニメーション用、オプション
   pytest>=7.0    # 開発用
   ```
8. **性能**: 24時間×dt=0.01 (8.6M ステップ) を 60秒以内で実行 (`numba.jit` 推奨)
9. **Python**: 3.9 以上 (`from __future__ import annotations` 推奨)
10. **コード規約**: black + ruff で整形・lint してから commit

## ドキュメント更新義務

1. `README.md`: `simulate.py` のセクション追加 (使い方の早見例)
2. `SIMULATION.md`: 「設計ガイド」から「実装後のリファレンス」に書き換え
3. `requirements.txt`: 新規作成
4. `THEORY.md`: 既存式とシミュレーション結果が一致する旨を §14.5 か新章で参照

## 推奨実装順 (3-5日)

```
Day 1: physics.py の中核 (回転合成、力の計算) + TC1, TC2 pass
Day 2: simulate.py 主ループ + CSV/JSON 出力 + TC3, TC4, TC5 pass
Day 3: visualize.py の 4-6 枚プロット + TC6, TC7, TC8 pass
Day 4: 細胞軌道シミュレーション + 3D アニメーション
Day 5: 性能最適化 (numba) + ドキュメント整備 + コードレビュー
```

各 Day 終了時に commit してください。Pull request にする必要はないので main に直接 push で OK (リポジトリの方針)。

## Git とコミット規約

- ブランチ: main に直接 commit (force push は禁止)
- author: `Rhizobium-gits <okay.bio.sato@gmail.com>` で固定
- Co-Authored-By: トレイラ **不要** (リポジトリ方針)
- コミットメッセージ: 1行サマリ + 詳細説明 (英語タイトル + 日本語詳細でも可、既存履歴に倣う)
- 1コミットで一つの論理単位、各コミットでテストが通ること

## 完了基準 (definition of done)

以下が全部 yes になったら完了:

- [ ] `python3 simulate.py --validate` で 8/8 pass
- [ ] `pytest -v tests/` で全 pass、カバレッジ 80% 以上
- [ ] `python3 simulate.py --mode tilt --g-target 0.378 --duration 1800 --animate --cells 100 --diffusion` がエラーなく完走し、CSV/PNG/MP4/JSON が全部生成される
- [ ] プロットが publication quality (タイトル/軸ラベル/凡例/単位/グリッド/colorblind-safe)
- [ ] シミュレーション結果と THEORY.md の解析値が 1% 以内一致 (短期平均) / 1e-5 g 以内 (長期平均)
- [ ] 48時間シミュレーションが 1分以内に完了 (numba 適用後)
- [ ] black + ruff の lint pass
- [ ] README/THEORY/SIMULATION 更新済み、新規ファイル `requirements.txt` 追加済み
- [ ] 各 commit でテストが通る粒度になっている

## 質問・前提確認

着手前に以下を確認してください:

1. リポジトリの既存ファイル (clinostat.py, THEORY.md, README.md, SIMULATION.md) を全部読んだか?
2. `compute_tilt()` などの既存関数の出力フォーマットを理解したか? (シミュレーションはこれを入力にする)
3. 物理モデルで不明な点があれば、推測せず質問してください (Coriolis の符号、座標系の向き、Stokes 抗力の係数 etc.)
4. もし不明確な仕様があったら自己解釈せずユーザに確認してください

## 補足: 最高品質のために

「動けば OK」ではなく以下を意識してください:

- **物理的厳密性**: 近似は明示。誤差オーダーを評価。
- **数値的正確性**: 浮動小数の誤差累積を診断 (orthonormality, energy)
- **検証可能性**: TC1〜TC8 だけでなく解析解との突合も適宜
- **再現性**: シード固定すれば完全再現、CI で連続実行可能
- **可視化品質**: 論文 figure として通用するレベル (Nature, PRL レベルの基準で)
- **拡張性**: 第5モード追加、別ハードウェア対応が容易
- **ドキュメント**: 6か月後に自分が読んでも理解できる

不明点や設計判断は遠慮なく質問してください。誤った前提で進むより、確認を優先してください。

実装、お願いします。

---PROMPT END---

## 使い方

1. 上の `---PROMPT START---` から `---PROMPT END---` までを、リポジトリのローカルクローンを開いた Claude Code セッションにそのまま貼り付ける
2. Claude が「読みました」と返したら「では実装してください」と短く促す
3. 完了基準のチェックリストを Claude に都度確認させる
4. PR にする必要なく main に直接 push (リポジトリ方針)

## カスタマイズ箇所

実情に応じて以下を変えてください:

- **培養生物**: 細菌 → 哺乳類細胞・酵母など (v_sed, D, T_bio, chamber size)
- **時間スケール**: 24h → 1週間など (--duration の最大想定値)
- **計算リソース**: numba 不可な環境なら性能要件を緩める
- **出力形式**: MP4 不可なら GIF だけに
- **ライセンス**: リポジトリに LICENSE がないので追加するか、現状維持か明示

## チェックリスト (依頼前に確認)

- [ ] リポジトリの最新の `main` を pull 済み
- [ ] `clinostat.py` と `THEORY.md` の最新版が手元にある
- [ ] `requirements.txt` がリポジトリにない場合、Claude が新規作成することを許可
- [ ] 開発 PC に Python 3.9+ と pip / venv 環境がある
- [ ] アニメーション機能を使うなら ffmpeg が入っている
