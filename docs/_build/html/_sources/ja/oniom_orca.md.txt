# `oniom-orca`

## 概要

> **要約:** Amber parm7 を用いて ORCA QM/MM（`.inp`）入力を生成します。

### クイックリファレンス
- **用途:** Amber 系のパラメータを保ったまま ORCA QM/MM を実行したい場合。
- **入力:** 必須 `--parm7`、任意 `-i/--input`（PDB/XYZ）、任意 `--model-pdb`。
- **出力:** QM 原子集合と ORCAFF 参照を含む ORCA 入力。
- **ORCAFF:** 未指定時は自動探索/生成。`--convert-orcaff`（デフォルト有効）。

## 最小例

```bash
mlmm oniom-orca --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp -q 0 -m 1
```

## 出力の見方

- `system.inp`
- `ORCAFF.prms`（既存利用または自動生成）

## よくある例

1. 基本出力。

```bash
mlmm oniom-orca --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp -q 0 -m 1
```

2. 全 QM+MM 系の総電荷/総多重度を明示。

```bash
mlmm oniom-orca --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp -q 0 -m 1 --total-charge -1 --total-mult 1
```

3. ORCAFF を明示指定し自動変換を無効化。

```bash
mlmm oniom-orca --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp --orcaff ./ORCAFF.prms --no-convert-orcaff
```

## 使用法

```bash
mlmm oniom-orca --parm7 real.parm7 [-i coords.pdb] [--model-pdb ml_region.pdb] \
 -o output.inp [-q CHARGE] [-m MULT] [--method "B3LYP D3BJ def2-SVP"] \
 [--total-charge INT] [--total-mult INT] [--nproc INT] [--near FLOAT] \
 [--orcaff PATH] [--convert-orcaff|--no-convert-orcaff] \
 [--element-check|--no-element-check]
```

## 説明

`oniom-orca` は `--parm7` からトポロジ情報を取得し、ORCA QM/MM 入力を書き出します。

主な処理:

1. parm7 から原子・結合・電荷情報を取得。
2. `-i/--input` がある場合は読み込み、`--element-check` で元素順を検証。
3. `--model-pdb` から QM 領域原子を同定。
4. `--orcaff` 指定または自動生成ルールで ORCAFF を解決。
5. method、並列設定、QM 選択、総電荷/総多重度を含む `.inp` を生成。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--parm7 PATH` | Amber parm7 トポロジ。 | 必須 |
| `-i, --input PATH` | 座標ファイル（PDB/XYZ）。原子順は parm7 と一致必須。 | _None_ |
| `--element-check / --no-element-check` | 入力と parm7 の元素順を検証。 | `True` |
| `--model-pdb PATH` | QM 領域原子を定義する PDB。 | _None_ |
| `-o, --output PATH` | 出力 ORCA 入力（`.inp`）。 | 必須 |
| `--method TEXT` | QM メソッド/基底。 | `B3LYP D3BJ def2-SVP` |
| `-q, --charge INT` | QM 領域電荷。 | 必須 |
| `-m, --multiplicity INT` | QM 領域多重度。 | `1` |
| `--total-charge INT` | 全 QM+MM 系の総電荷（`Charge_Total`）。 | トポロジ由来 |
| `--total-mult INT` | 全 QM+MM 系の総多重度（`Mult_Total`）。 | `--multiplicity` と同じ |
| `--nproc INT` | 使用コア数。 | `8` |
| `--near FLOAT` | ActiveAtoms 判定カットオフ（層タグなし時）。 | `6.0` |
| `--orcaff PATH` | ORCAFF.prms のパス。 | 自動 |
| `--convert-orcaff/--no-convert-orcaff` | ORCAFF 未検出時の `orca_mm -convff -AMBER` 実行可否。 | `True` |

## 出力

```text
<output>.inp
```

## 注意事項

- `parmed` が必要です。
- 元素検証はベストエフォートです（原子順不変を仮定）。
- ORCA のリンク原子は QM/MM 境界情報から自動処理されます。

---

## 関連項目

- [oniom-gaussian](oniom_gaussian.md) -- Gaussian ONIOM エクスポータ
- [oniom_export](oniom_export.md) -- エクスポート全体ガイド
- [mm_parm](mm_parm.md) -- Amber トポロジ構築
- [define_layer](define_layer.md) -- レイヤー定義/確認
