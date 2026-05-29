# clinostat-calc

2軸クリノスタット (内側 r = 8.5 cm / 外側 r = 14.5 cm) の回転速度を計算する CLI ツール。

## モデル

- **内側フレーム** … 遠心力で目標 g を生成
  `ω_in = √(g_target · g₀ / r_in)`
- **外側フレーム** … 黄金比で非共鳴クリノローテーション
  `ω_out = ω_in / φ`,  `φ = (1+√5)/2`

`g₀ = 9.81 m/s²`,  `r_in = 0.085 m`,  `r_out = 0.145 m`

## 使い方

```bash
# 単発計算 (数値 g 単位)
python3 clinostat.py 0.378

# プリセット名
python3 clinostat.py mars

# m/s² 単位入力
python3 clinostat.py "3.72 m/s2"

# プリセット一覧を一気に表示
python3 clinostat.py --all

# 対話モード (引数なし)
python3 clinostat.py
```

### 対話モードのコマンド

| 入力                 | 動作                                     |
|----------------------|------------------------------------------|
| `0.166` / `1.0`      | g 単位で計算                             |
| `mars` / `moon` …    | プリセット計算                           |
| `3.72 m/s2`          | m/s² 単位で計算                          |
| `list` / `presets`   | プリセット一覧表示                       |
| `h` / `help` / `?`   | プリセット名簡易表示                     |
| `q` / `quit` / `exit`| 終了                                     |

## プリセット

`earth, mars, moon, ceres, europa, titan, iss, ug`

## 注意

外側フレーム位置での遠心加速度が目標 g の 50% を超えると `[!] 外乱大` 警告を出します。
`g_target` が 1g に近いほど黄金比方式の外乱寄与は大きくなります。
