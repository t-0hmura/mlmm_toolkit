# クイックスタート: `mlmm all`

## 目的

反応物・生成物のPDBから、ML領域の選択、MMトポロジー・層の準備、MEP探索を実行します。
MEPは単一パスの `path-opt` が既定で、`--refine-path` で再帰的な `path-search` を使います。
TS最適化・IRC・熱化学・DFTは必要に応じて追加できます。

## 事前に必要なもの

- 同じ原子を同じ順序で含む、水素付きの完全系PDBを2つ（R/P）。PyMOLで保存する場合は *Original atom order* を有効にします。
- `-l RES:CHARGE` は実際の水素数と整合させます。例えばSAMは水素23個で `SAM:1`、22個で `SAM:0`。不整合はantechamberの電子数エラーの原因になります。
- GPUを推奨します。MLIPは既定の `uma` のほか、`-b` で `orb` / `mace` / `aimnet2` を選べます。

## 最小コマンド

[公開サンプル](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples/toy_system)の
`r_complex.pdb` と `p_complex.pdb` を同じ作業ディレクトリに保存し、そこで実行します。

```bash
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 \
  --ligand-charge 'PRE:0' -q -1 -m 1 --out-dir ./result_all
```

TS最適化・IRC・熱化学・DFTまで追加する場合:

```bash
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 \
  --ligand-charge 'PRE:0' -q -1 -m 1 \
  --tsopt --thermo --dft --out-dir ./result_all
```

`-c` はML領域の中心残基、`-r` は抽出半径（Å）、`-q` / `-m` はML領域の電荷・多重度です。

## 出力の検証

エネルギーや結合変化を解釈する前に、`summary.json` の[実行結果と理由](json-output.md#実行と要求段階の完了状況)を確認します。
`summary.log` に結果の要約、出力ルートに `mep.pdb`（bridge入力では `mep.cif` も）と
`energy_diagram_MEP.png` を保存します。
生出力は `_work/path_opt/`（`--refine-path` 時は `_work/path_search/`）にあります。
全体の出力ツリーは [all](all.md)、ファイル名は [出力構造](output-layout.md) を参照してください。

## 補足

- `--dry-run` で引数と実行計画を確認できます。
- `mlmm all --help` は主要オプション、`mlmm all --help-advanced` は全オプションを表示します。

## 次のステップ

- 単一構造スキャン: [クイックスタート: scan](quickstart-scan-spec.md)
- TS検証: [クイックスタート: tsopt → freq](quickstart-tsopt-freq.md)
- エラー対処: [典型エラー別レシピ](recipes-common-errors.md) · [トラブルシューティング](troubleshooting.md)
