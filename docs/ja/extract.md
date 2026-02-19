# `extract`

## 概要

> **要約:** タンパク質-基質複合体から指定された基質残基周辺の活性部位ポケットを抽出します。生化学的に適切な切断（主鎖/側鎖キャッピング）を適用し、切断結合にリンク水素を付加できます。

`mlmm extract` は、指定された基質残基周辺の活性部位ポケット（結合ポケット）を自動抽出するツールです。単一構造およびアンサンブル（マルチモデルまたはファイルごと）に対応しています。QM/MM、ML/MM、およびクラスター QM モデルの構築に使用します。

### 概要
- **入力:** 1 つ以上の完全タンパク質-基質複合体 PDB ファイル。
- **出力:** 活性部位ポケット PDB。任意でリンク水素を含みます。
- **用途:** ML/MM 計算や QM/MM 計算の前処理として、活性部位周辺のモデル構造を準備。

### 残基包含規則
- 基質残基は常に含まれます。
- 標準カットオフ（`--radius`、デフォルト 2.6 A）:
  - `--exclude-backbone false` の場合: **任意の原子**がカットオフ内にある残基を含めます。
  - `--exclude-backbone true`（デフォルト）の場合: **アミノ酸残基**では、適格原子は**非主鎖**（N, H*, CA, HA*, C, O, OXT 以外）でなければなりません。非アミノ酸残基は任意の原子で適格となります。
- 独立したヘテロ-ヘテロ近接（`--radius-het2het`）: **基質ヘテロ原子（C/H 以外）** がカットオフ内の**タンパク質ヘテロ原子**に近接する残基を追加します。
- 水分子はデフォルトで含まれます（`--include-H2O true`; HOH/WAT/TIP3/SOL）。
- `--selected-resn` で残基を強制包含できます（鎖とインサーションコードに対応）。

## 使用法

```bash
mlmm extract -i INPUT.pdb [INPUT2.pdb ...] -c <substrate_spec> \
    [-o OUTPUT.pdb ...] [-r <A>] [--radius-het2het <A>] \
    [--include-H2O {true|false}] [--exclude-backbone {true|false}] \
    [--add-linkH {true|false}] [--selected-resn "CHAIN:RES" ...] \
    [--ligand-charge <number|"RES:Q,...">] [--verbose {true|false}]
```

### 例

```bash
# ID ベースの基質と明示的リガンド電荷によるミニマル実行
mlmm extract -i complex.pdb -c A:123 -o pocket.pdb --ligand-charge -3

# PDB で基質を指定; 残基名ごとの電荷マッピング（その他は 0）
mlmm extract -i complex.pdb -c substrate.pdb -o pocket.pdb \
    --ligand-charge "GPP:-3,MMT:-1"

# 名前ベースの基質選択（すべてのマッチを含む）
mlmm extract -i complex.pdb -c "GPP,MMT" -o pocket.pdb --ligand-charge -4

# マルチ構造から単一マルチモデル出力、ヘテロ-ヘテロ近接有効
mlmm extract -i complex1.pdb complex2.pdb -c A:123 \
    -o pocket_multi.pdb --radius-het2het 2.6 --ligand-charge -3 --verbose true
```

## ワークフロー

1. **入力構造の読み込み** -- 1 つ以上の完全タンパク質-基質複合体 PDB を読み込みます。
2. **基質の特定** -- `-c/--center` で指定された基質を PDB パス、残基 ID（例: `A:123` または `123,124`）、または残基名（例: `GPP,MMT`）で同定します。
3. **ポケット残基の選択** -- カットオフ距離に基づいて周囲の残基を選択します。`--exclude-backbone` が有効な場合、アミノ酸の主鎖原子は除外されます。
4. **リンク水素の付加** -- `--add-linkH true` の場合、切断された共有結合にリンク水素を付加します。
5. **電荷の導出** -- `--ligand-charge` が指定された場合、アミノ酸、イオン、リガンド電荷を合算してポケットの総電荷を算出します。
6. **出力の書き込み** -- ポケット PDB を指定されたパスに書き出します。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 入力 PDB ファイル（1 つ以上）。 | 必須 |
| `-c, --center TEXT` | 基質指定（PDB パス、残基 ID、または残基名）。 | 必須 |
| `-o, --output PATH...` | 出力ポケット PDB パス。 | _None_ |
| `-r, --radius FLOAT` | ポケット包含カットオフ (A)。 | `2.6` |
| `--radius-het2het FLOAT` | ヘテロ-ヘテロカットオフ (A)。 | `0.0` |
| `--include-H2O {true\|false}` | 水分子を含める（HOH/WAT/TIP3/SOL）。 | `true` |
| `--exclude-backbone {true\|false}` | 非基質アミノ酸の主鎖原子を除外。 | `true` |
| `--add-linkH {true\|false}` | 切断結合にリンク水素を付加。 | `true` |
| `--selected-resn TEXT` | 強制包含する残基。 | `""` |
| `--ligand-charge TEXT` | 総電荷または残基名ごとのマッピング（例: `"GPP:-3,MMT:-1"`）。 | _None_ |
| `--verbose {true\|false}` | 詳細ログを有効化。 | `true` |

## 出力

```
<output>.pdb                # 活性部位ポケット PDB（リンク水素を含む場合あり）
(stdout)                    # ポケット情報サマリー（残基数、原子数、電荷）
```

## 注意事項

- 基質は PDB パス、残基 ID（`A:123` や `123,124`）、または残基名（`GPP,MMT`）で指定できます。
- 名前ベースの指定では同名の残基がすべて選択されます。意図しない選択に注意してください。
- `--radius` または `--radius-het2het` に `0` を渡すと、内部で `0.001 A` にクランプされます。
- マルチ構造入力の場合、残基選択は全構造の和集合として統合されます。

---

## 関連項目

- [all](all.md) -- エンドツーエンドワークフロー（内部で `extract` を呼び出し）
- [mm_parm](mm_parm.md) -- 抽出後に AMBER トポロジーを構築
- [add-elem-info](add_elem_info.md) -- PDB 元素記号の修復（抽出前に必要な場合）
