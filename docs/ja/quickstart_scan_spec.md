# クイックスタート: `mlmm scan` + `--spec`

## 目的

単一構造に対して、YAML で定義した距離ターゲットへ段階的にスキャンします。

## 事前に必要なファイル

- 入力構造: `pocket.pdb`
- MM トポロジ: `real.parm7`
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
mlmm scan -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --spec scan.yaml --print-parsed --out-dir ./result_scan
```

## まず確認する出力

- `result_scan/stage_01/result.pdb`
- `result_scan/stage_02/result.pdb`
- `--dump` 指定時は `scan_trj.xyz` / `scan.pdb`

## 補足

- 詳細オプションは `mlmm scan --help-advanced` で確認できます。

## 次の導線

- 経路精密化は [all](all.md) または [path-search](path_search.md) を参照してください。
