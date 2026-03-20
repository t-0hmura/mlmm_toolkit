# クイックスタート: `mlmm all`

## 目的

2 つの完全系 PDB から、end-to-endのワークフローを 1 回実行します。

## 最小コマンド

```bash
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --out-dir ./result_all
```

後処理（TS 最適化、熱化学、DFT）まで同時に実行する場合:

```bash
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft --out-dir ./result_all
```

## まず確認する出力

- `result_all/summary.log`
- `result_all/summary.yaml`
- `result_all/path_search/mep.pdb`（または `result_all/path_search/seg_*/`）

## 補足

- `--dry-run`（`--help-advanced` に表示）で引数と実行計画を確認できます。
- `mlmm all --help` は主要オプション、`mlmm all --help-advanced` は `--dry-run` 等を含む全オプションを表示します。
- 別の MLIP バックエンドを使用するには、`-b orb`（または `mace`、`aimnet2`）を追加します。デフォルトは `uma` です。
- `--embedcharge` を追加すると、MM 環境から ML 領域への静電影響を考慮する xTB 点電荷埋め込み補正が有効になります。

## 次の導線

- 単一構造スキャン: [クイックスタート: `mlmm scan` + `-s`（YAML スペック）](quickstart_scan_spec.md)
- TS 検証: [クイックスタート: `mlmm tsopt` -> `mlmm freq`](quickstart_tsopt_freq.md)
- 全オプション: [all](all.md)
