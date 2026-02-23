# `oniom-gaussian`

## 概要

> **要約:** Amber parm7 を用いて Gaussian ONIOM（`.com`/`.gjf`）入力を生成します。

### クイックリファレンス
- **用途:** ML/MM で準備した系を Gaussian ONIOM に渡したい場合。
- **入力:** 必須 `--parm7`、任意 `-i/--input`（PDB/XYZ）、任意 `--model-pdb`。
- **出力:** レイヤー情報と結合情報を含む Gaussian 入力。
- **検証:** `--element-check`（デフォルト）で座標ファイルと parm7 の元素順を照合。

## 最小例

```bash
mlmm oniom-gaussian --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com -q 0 -m 1
```

## 出力の見方

- `system.com` または `system.gjf`
- コンソールに QM 原子数や境界処理の要約

## よくある例

1. メソッドを明示して出力。

```bash
mlmm oniom-gaussian --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com --method "B3LYP/6-31G(d,p)"
```

2. 元素順チェックを無効化。

```bash
mlmm oniom-gaussian --parm7 real.parm7 -i pocket.xyz --model-pdb ml_region.pdb \
 -o system.gjf --no-element-check
```

3. 実行環境パラメータを調整。

```bash
mlmm oniom-gaussian --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com --nproc 16 --mem 32GB --near 5.0
```

## 使用法

```bash
mlmm oniom-gaussian --parm7 real.parm7 [-i coords.pdb] [--model-pdb ml_region.pdb] \
 -o output.com [-q CHARGE] [-m MULT] [--method "B3LYP/6-31G(d,p)"] \
 [--near FLOAT] [--nproc INT] [--mem TEXT] \
 [--element-check|--no-element-check]
```

## 説明

`oniom-gaussian` は `--parm7` からトポロジ情報を読み、必要に応じて `-i/--input` の座標を使って Gaussian ONIOM 入力を書き出します。

主な処理:

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
| `-q, --charge INT` | QM 領域電荷。 | `0` |
| `-m, --mult INT` | QM 領域多重度。 | `1` |
| `--near FLOAT` | 可動 MM 原子を決める距離カットオフ。 | `6.0` |
| `--nproc INT` | 使用コア数。 | `8` |
| `--mem TEXT` | メモリ指定。 | `16GB` |

## 出力

```text
<output>.com / <output>.gjf
```

## 注意事項

- `parmed` が必要です。
- 元素検証はベストエフォートです（原子順不変を仮定）。
- `--near` は可動 MM 原子の範囲指定に使われます。

---

## 関連項目

- [oniom-orca](oniom_orca.md) -- ORCA QM/MM エクスポータ
- [oniom_export](oniom_export.md) -- エクスポート全体ガイド
- [mm_parm](mm_parm.md) -- Amber トポロジ構築
- [define_layer](define_layer.md) -- レイヤー定義/確認
