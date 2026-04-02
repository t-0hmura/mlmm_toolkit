# クイックスタート: scan

## 目的

単一構造に対して、YAML で定義した距離ターゲットへ段階的にスキャンします。

## 事前に必要なファイル

- 入力構造: `pocket.pdb`
- MM トポロジー: `real.parm7`
- ML 領域定義: `ml_region.pdb`

これらは通常、`mlmm all` / `mlmm extract` / `mlmm mm-parm` で生成します。

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
 -q 0 -s scan.yaml --print-parsed -o ./result_scan
```

## まず確認する出力

- `result_scan/stage_01/result.pdb`
- `result_scan/stage_02/result.pdb`
- `--dump` 指定時は `scan_trj.xyz` / `scan.pdb`

## 補足

- 詳細オプションは `mlmm scan --help-advanced` で確認できます。

## インラインリテラル入力（YAML ファイルなし）

YAML スペックファイルの代わりに、スキャンターゲットをコマンドラインで直接指定できます:

```bash
mlmm scan -i layered.pdb --parm system.parm7 -q 0 \
  --scan-lists "[(1,5,1.4)]" --no-preopt --no-endopt
```

PDB 原子セレクタも使用可能です:

```bash
mlmm scan -i layered.pdb --parm system.parm7 -q 0 \
  --scan-lists "[('TYR 285 CA','MMT 309 C10',2.20)]" --no-preopt --no-endopt
```

1-based の原子インデックスと PDB 原子名文字列の両方が使用できます。詳細は [scan.md](scan.md) を参照してください。

## 次の導線

- 経路精密化は [all](all.md) または [path-search](path-search.md) を参照してください。
