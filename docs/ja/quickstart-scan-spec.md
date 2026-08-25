# クイックスタート: scan

## 目的

単一構造を起点に、YAMLで定義した距離ターゲットへ拘束付きスキャンを行い、端点構造と軌跡を生成します。

## 事前に必要なファイル

- 入力構造: `pocket.pdb`
- MM トポロジー: `real.parm7`
- ML 領域定義: `ml_region.pdb`、明示的model index、または有効なB-factor layer

これらは通常、`mlmm all` / `mlmm extract` / `mlmm mm-parm` で生成します。

YAMLの1つの`stages`要素が1ステージです。同一要素内の複数距離tupleは
協奏的に駆動し、複数要素は多段階scanとして順次実行されます。2距離を独立なグリッド軸と
する場合は`scan2d`を使用します。

## 1. `scan.yaml` を作成

```yaml
one_based: true
stages:
 - [[12, 45, 2.20]]
 - [[10, 55, 1.35], [23, 34, 1.80]]
```

## 2. 実行

```bash
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan.yaml -o ./result_scan
```

## 出力の検証

- `result_scan/stage_01/result.pdb`
- `result_scan/stage_02/result.pdb`
- `result_scan/scan_trj.xyz`（常に出力）。`result_scan/scan.pdb` は変換と参照トポロジーを利用できる場合

## 補足

- 実行せずに（GPU 不要で）スペックを検証したい場合は `--print-parsed` を付けます。解析されたターゲットを表示して計算を実行する前に終了するため、上記のスキャン出力は生成されません。
- 詳細オプションは `mlmm scan --help-advanced` で確認できます。

## インラインリテラル入力（YAML ファイルなし）

YAML スペックファイルの代わりに、スキャンターゲットをコマンドラインで直接指定できます:

```bash
mlmm scan -i layered.pdb --parm system.parm7 -q 0 \
  --scan-lists '[(1,5,1.4)]' --no-preopt --no-endopt
```

PDB 原子セレクタも使用可能です:

```bash
mlmm scan -i layered.pdb --parm system.parm7 -q 0 \
  --scan-lists '[("TYR,285,CA","MMT,309,C10",2.20)]' --no-preopt --no-endopt
```

1-based の原子インデックスと PDB 原子名文字列の両方が使用できます。詳細は [scan.md](scan.md) を参照してください。

## 次の導線

- 経路精密化は [all](all.md) または [path-search](path-search.md) を参照してください。
