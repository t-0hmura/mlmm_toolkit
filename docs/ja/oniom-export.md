# `oniom-export`

Amber トポロジーを持つ ML/MM 系を、外部 QM/MM 入力ファイル（Gaussian ONIOM = `--mode g16`、または ORCA QM/MM = `--mode orca`）へエクスポートします。`mlmm oniom-export` は Amber `parm7` トポロジーと座標ファイル、ML 領域（QM 領域）定義を読み込み、そのまま実行できる入力ファイルを 1 つ書き出します。QM 領域は `--model-pdb` で指定し、周囲の MM 環境は対象プログラムのネイティブ形式で、QM/MM 切断面にリンク原子注釈を付けて出力します。

## 実行例

```bash
# Gaussian ONIOM 入力
mlmm oniom-export --parm real.parm7 -i pocket.pdb --model-pdb ml.pdb \
 -o out.gjf --mode g16 -q 0 -m 1
```

```bash
# ORCA QM/MM 入力（.inp 拡張子からモード推定）
mlmm oniom-export --parm real.parm7 -i pocket.pdb --model-pdb ml.pdb \
 -o out.inp -q 0 -m 1
```

```bash
# メソッド/基底とリソースを指定した Gaussian 入力
mlmm oniom-export --parm real.parm7 -i pocket.pdb --model-pdb ml.pdb \
 -o out.gjf --mode g16 --method 'wb97xd/def2-svp' --nproc 16 --mem 32GB -q 0 -m 1
```

## 処理の流れ

1. **トポロジー + 座標** -- `parm7` と `-i` 座標ファイルを読み込みます（原子順序はトポロジーと一致が必須。`--element-check` が元素配列を検証）。
2. **QM 領域** -- `--model-pdb` で QM（ML 領域）原子を定義し、`--near` で可動/活性 MM のカットオフ（Å）を設定します。
3. **リンク原子** -- 切断された QM/MM 結合ごとに配置されます。`--link-atom-method scaled`（デフォルト）は Morokuma/Dapprich の g-factor（`MLMMCore` ランタイムと一致）、`fixed` は固定 1.09/1.01 Å を使用します。
4. **書き出し** -- `-o` に対象形式の入力ファイルを出力します。ORCA モードでは `ORCAFF.prms` のパスを特定します（`--convert-orcaff` が有効なら Amber から `orca_mm -convff -AMBER` で自動変換）。

## 出力

- `<output>.{gjf,com}`（g16）または `<output>.inp`（ORCA） -- QM/MM 入力ファイル
- ORCA モードでは出力ディレクトリの `<parm7_stem>.ORCAFF.prms`（力場パラメータ）も読み込み/生成します

## CLI オプション

全フラグの一覧は生成済みの[コマンドリファレンス](../reference/commands/index.md)にあります。以下の表は説明が必要なオプションを扱います。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--parm FILE` | Amber parm7 トポロジーファイル | 必須 |
| `-i, --input FILE` | 現在構造の座標ファイル（`.pdb` / `.xyz`） | _None_ |
| `--model-pdb FILE` | QM 領域原子を定義する PDB | _None_ |
| `-o, --output FILE` | 出力ファイルパス（g16 は `.gjf` / `.com`、ORCA は `.inp`） | 必須 |
| `--mode [g16\|orca]` | エクスポートモード。未指定時は `-o` 拡張子から推定 | _推定_ |
| `--method TEXT` | QM メソッドと基底関数 | モード依存 |
| `-q, --charge INT` | QM 領域の電荷 | 必須 |
| `-m, --multiplicity INT` | QM 領域の多重度 | `1` |
| `--near FLOAT` | 可動/活性 MM 原子の距離カットオフ（Å） | `6.0` |
| `--nproc INT` | プロセッサ数 | `8` |
| `--mem TEXT` | メモリ割り当て（g16 モード） | `16GB` |
| `--total-charge INT` / `--total-mult INT` | 全 QM+MM 系の総電荷/総多重度（ORCA `Charge_Total` / `Mult_Total`） | _None_ |
| `--orcaff PATH` | `ORCAFF.prms` のパス（ORCA モード）。未指定時は出力ディレクトリに生成 | _None_ |
| `--convert-orcaff / --no-convert-orcaff` | `ORCAFF.prms` 欠損時に `orca_mm -convff -AMBER` で自動変換（ORCA モード） | `True` |
| `--element-check / --no-element-check` | `--input` の元素配列を parm7 トポロジーと照合 | `True` |
| `--link-atom-method [scaled\|fixed]` | リンク H 配置: `scaled`（g-factor、ランタイム一致）または `fixed`（1.09/1.01 Å） | `scaled` |

`mlmm oniom-export --help` はコアオプション、`mlmm oniom-export --help-advanced` は全オプションを表示します。

## 注記

- モード選択: `--mode` が最優先です。`--mode` 未指定時は `-o` から推定します。
  - `.gjf` / `.com` -> `g16`
  - `.inp` -> `orca`
- `--mode` 未指定かつ `-o` が未知拡張子の場合はエラーになります。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [oniom_gaussian](oniom-gaussian.md) -- Gaussian モード詳細（`--mode g16`）
- [oniom_orca](oniom-orca.md) -- ORCA モード詳細（`--mode orca`）
- [oniom_import](oniom-import.md) -- ONIOM 入力から XYZ/層付き PDB を再構築
- [mm_parm](mm-parm.md) -- Amber トポロジー構築
- [define_layer](define-layer.md) -- レイヤー定義/確認
