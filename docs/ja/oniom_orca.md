# ORCA QM/MM モード（`oniom-export --mode orca`）

## 概要

> **要約:** Amber parm7 を用いて ORCA QM/MM（`.inp`）入力を生成します。

## 最小例

```bash
mlmm oniom-export --mode orca --parm real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp -q 0 -m 1
```

## 出力の見方

- `system.inp`
- `ORCAFF.prms`（既存利用または自動生成）

## よくある例

1. 基本出力。

```bash
mlmm oniom-export --mode orca --parm real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp -q 0 -m 1
```

2. 全 QM+MM 系の総電荷/総多重度を明示。

```bash
mlmm oniom-export --mode orca --parm real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp -q 0 -m 1 --total-charge -1 --total-mult 1
```

3. ORCAFF を明示指定し自動変換を無効化。

```bash
mlmm oniom-export --mode orca --parm real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp --orcaff ./ORCAFF.prms --no-convert-orcaff
```

## ワークフロー

ORCA モード（`mlmm oniom-export --mode orca`）は `--parm` からトポロジー情報を取得し、ORCA QM/MM 入力を書き出します。

1. parm7 から原子・結合・電荷情報を取得。
2. `-i/--input` がある場合は読み込み、`--element-check` で元素順を検証。
3. `--model-pdb` から QM 領域原子を同定。
4. `--orcaff` 指定または自動生成ルールで ORCAFF を解決。
5. method、並列設定、QM 選択、総電荷/総多重度を含む `.inp` を生成。

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
| `--near FLOAT` | ActiveAtoms 判定カットオフ（層タグなし時）。 | `6.0` |
| `--orcaff PATH` | ORCAFF.prms のパス。 | 自動 |
| `--convert-orcaff/--no-convert-orcaff` | ORCAFF 未検出時の `orca_mm -convff -AMBER` 実行可否。 | `True` |

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- `parmed` が必要です。
- 元素検証はベストエフォートです（原子順不変を仮定）。
- ORCA のリンク原子は QM/MM 境界情報から自動処理されます。

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [oniom_gaussian](oniom_gaussian.md) -- Gaussian モードガイド（`--mode g16`）
- [oniom_export](oniom_export.md) -- エクスポート全体ガイド
- [mm_parm](mm_parm.md) -- Amber トポロジー構築
- [define_layer](define_layer.md) -- レイヤー定義/確認
- ORCA 6.0 マニュアル（QM/MM）: <https://www.faccts.de/docs/orca/6.0/manual/contents/typical/qmmm.html>
