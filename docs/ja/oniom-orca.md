# ORCA QM/MM モード（`oniom-export --mode orca`）

Amber parm7 を用いて ORCA QM/MM（`.inp`）入力を生成します。これは `oniom-export` の ORCA 詳細モード（`mlmm oniom-export --mode orca`）です。Amber parm7 でトポロジーを記述した ML/MM 系について、ORCA QM/MM（`.inp`）入力が必要なときに使います。parm7 ファイルからトポロジー情報を読み込み、model 領域の PDB を QM 原子へ対応付け、ORCAFF パラメータを解決し、ORCA 6.0 実行向けの単一 ORCA QM/MM 入力ファイルを書き出します。

## 実行例

基本的なエクスポート:

```bash
mlmm oniom-export --mode orca --parm real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp -q 0 -m 1
```

全 QM+MM 系の総電荷/総多重度を明示的に指定:

```bash
mlmm oniom-export --mode orca --parm real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp -q 0 -m 1 --total-charge -1 --total-mult 1
```

ORCAFF のパスを明示し自動変換を無効化:

```bash
mlmm oniom-export --mode orca --parm real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp --orcaff ./ORCAFF.prms --no-convert-orcaff
```

## 処理の流れ

ORCA モード（`mlmm oniom-export --mode orca`）は `--parm` からトポロジー情報を取得し、ORCA QM/MM 入力を書き出します。

1. parm7 から原子・結合・電荷情報を取得。
2. `-i/--input` がある場合は読み込み、`--element-check` で元素順を検証。
3. `--model-pdb` から QM 領域原子を同定。
4. `--orcaff` 指定または自動生成ルールで ORCAFF を解決。
5. method、並列設定、QM 選択、総電荷/総多重度を含む `.inp` を生成。

## 出力

- `system.inp`
- `ORCAFF.prms`（既存利用または自動生成）

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--parm PATH` | Amber parm7 トポロジー。 | 必須 |
| `-i, --input PATH` | 座標ファイル（PDB/XYZ）。原子順は parm7 と一致必須。 | _None_ |
| `--element-check / --no-element-check` | 入力と parm7 の元素順を検証。 | `True` |
| `--model-pdb PATH` | QM 領域原子を定義する PDB。 | _None_ |
| `-o, --output PATH` | 出力 ORCA 入力（`.inp`）。 | 必須 |
| `--method TEXT` | QM メソッド/基底。 | `B3LYP D3BJ def2-SVP` |
| `-q, --charge INT` | QM 領域電荷。 | 必須 |
| `-m, --multiplicity INT` | QM 領域多重度。 | `1` |
| `--total-charge INT` | 全 QM+MM 系の総電荷（`Charge_Total`）。 | トポロジー由来 |
| `--total-mult INT` | 全 QM+MM 系の総多重度（`Mult_Total`）。 | `--multiplicity` と同じ |
| `--nproc INT` | 使用コア数。 | `8` |
| `--near FLOAT` | ActiveAtoms 判定カットオフ (Å)（層タグなし時）。 | `6.0` |
| `--orcaff PATH` | ORCAFF.prms のパス。 | _None_（自動で `<parm7名>.ORCAFF.prms` を解決） |
| `--convert-orcaff/--no-convert-orcaff` | ORCAFF 未検出時の `orca_mm -convff -AMBER` 実行可否。 | `True` |

## 関連項目

- [oniom_gaussian](oniom-gaussian.md) -- Gaussian モードガイド（`--mode g16`）
- [oniom_export](oniom-export.md) -- エクスポート全体ガイド
- [mm_parm](mm-parm.md) -- Amber トポロジー構築
- [define_layer](define-layer.md) -- レイヤー定義/確認
- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- ORCA 6.0 マニュアル（QM/MM）: <https://www.faccts.de/docs/orca/6.0/manual/contents/typical/qmmm.html>
