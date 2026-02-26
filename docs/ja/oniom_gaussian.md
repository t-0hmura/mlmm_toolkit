# Gaussian ONIOM モード（`oniom-export --mode g16`）

## 概要

> **要約:** Amber parm7 を用いて Gaussian ONIOM（`.com`/`.gjf`）入力を生成します。

### 概要
- **用途:** ML/MM で準備した系を Gaussian ONIOM に渡したい場合。
- **手法:** ParmEd で Amber parm7 トポロジを読み込み、モデル PDB で QM 層を定義、QM/MM 境界を検出しリンク原子を注釈。
- **出力:** メソッド、座標、レイヤーフラグ、結合情報を含む Gaussian `.com` または `.gjf`。
- **デフォルト:** `--method "B3LYP/6-31G(d,p)"`、`--nproc 8`、`--mem 16GB`、`--near 6.0`。
- **次のステップ:** Gaussian でエクスポートした入力を実行。

## 最小の例

```bash
mlmm oniom-export --mode g16 --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com -q 0 -m 1
```

## 出力チェックリスト

- `system.com` または `system.gjf`
- コンソールに QM 原子数や境界処理の要約

## よくある例

1. メソッドを明示して出力。

```bash
mlmm oniom-export --mode g16 --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com --method "B3LYP/6-31G(d,p)"
```

2. 元素順チェックを無効化。

```bash
mlmm oniom-export --mode g16 --parm7 real.parm7 -i pocket.xyz --model-pdb ml_region.pdb \
 -o system.gjf --no-element-check
```

3. 実行環境パラメータを調整。

```bash
mlmm oniom-export --mode g16 --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com --nproc 16 --mem 32GB --near 5.0
```

## ワークフロー

Gaussian モード（`mlmm oniom-export --mode g16`）は `--parm7` からトポロジ情報（ParmEd 経由）を読み、必要に応じて `-i/--input` の座標を使って Gaussian ONIOM 入力を書き出します。

1. parm7 から原子・結合・電荷情報を取得。
2. `-i/--input` がある場合は読み込み、`--element-check` で元素順を検証。
3. `--model-pdb` がある場合は QM 領域をマッピング。
4. QM/MM 境界を検出し、リンク原子情報を付与。
5. `%nproc`/`%mem`/method/座標/レイヤー/結合情報を含む入力を生成。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--parm7 PATH` | Amber parm7 トポロジ。 | 必須 |
| `-i, --input PATH` | 座標ファイル（PDB/XYZ）。原子順は parm7 と一致必須。 | _None_ |
| `--element-check / --no-element-check` | 入力と parm7 の元素順を検証。 | `True` |
| `--model-pdb PATH` | QM 領域原子を定義する PDB。 | _None_ |
| `-o, --output PATH` | 出力 Gaussian 入力（`.com`/`.gjf`）。 | 必須 |
| `--method TEXT` | QM メソッド/基底。 | `B3LYP/6-31G(d,p)` |
| `-q, --charge INT` | QM 領域電荷。 | 必須 |
| `-m, --multiplicity INT` | QM 領域多重度。 | `1` |
| `--near FLOAT` | 可動 MM 原子を決める距離カットオフ。 | `6.0` |
| `--nproc INT` | 使用コア数。 | `8` |
| `--mem TEXT` | メモリ指定。 | `16GB` |

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- `parmed` が必要です。
- 元素検証はベストエフォートです（原子順不変を仮定）。
- `--near` は可動 MM 原子の範囲指定に使われます。

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [oniom_orca](oniom_orca.md) -- ORCA モードガイド（`--mode orca`）
- [oniom_export](oniom_export.md) -- エクスポート全体ガイド
- [mm_parm](mm_parm.md) -- Amber トポロジ構築
- [define_layer](define_layer.md) -- レイヤー定義/確認
