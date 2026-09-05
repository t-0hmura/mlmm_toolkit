# クイックスタート: `mlmm all`

## 目的

2 つの完全系 PDB から、一気通貫のワークフローを 1 回実行します。

以下の `R.pdb` / `P.pdb` は、自分で準備した水素付きの構造です。同じ原子を同じ順序で含め、リガンド電荷を水素数に合わせてください（[入力の準備](getting-started.md#概要)）。

## 最小コマンド

```bash
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --out-dir ./result_all
```

後処理（TS 最適化、熱化学、DFT）まで同時に実行する場合:

```bash
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft --out-dir ./result_all
```

## 出力の検証

- `result_all/summary.log`
- `result_all/summary.json` — エネルギーを解釈する前に[実行結果と理由](json-output.md#実行結果と科学的妥当性)を確認
- `result_all/mep.pdb`（bridge 入力では `mep.cif` もルートに移動）と生出力 `result_all/_work/path_opt/`（`--refine-path` 時は `_work/path_search/`）

## 補足

- `--dry-run` で引数と実行計画を確認してから実行できます。
- `mlmm all --help` は主要オプション、`mlmm all --help-advanced` は詳細オプションも含めた全オプションを表示します。
- 別の MLIP バックエンドを使用するには、`-b orb`（または `mace`、`aimnet2`）を追加します。デフォルトは `uma` です。

## 次のステップ

- 単一構造スキャン: [クイックスタート: `mlmm scan` + `-s`（YAML スペック）](quickstart-scan-spec.md)
- TS 検証: [クイックスタート: `mlmm tsopt` -> `mlmm freq`](quickstart-tsopt-freq.md)
- 全オプション: [all](all.md)
