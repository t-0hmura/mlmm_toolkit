# クイックスタート: `mlmm tsopt` -> `mlmm freq`

## 目的

TS 候補を最適化し、振動解析で妥当性を確認します。

## 事前に必要なファイル

- TS 候補構造: `ts_guess.pdb`
- MM トポロジ: `real.parm7`
- ML 領域定義: `ml_region.pdb`

## 1. TS 最適化

```bash
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_tsopt
```

## 2. 最適化結果に対して振動解析

```bash
mlmm freq -i ./result_tsopt/final_geometry.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_freq
```

## まず確認する出力

- `result_tsopt/final_geometry.pdb`
- `result_freq/frequencies_cm-1.txt`
- `result_freq/mode_*_trj.xyz` と `result_freq/mode_*.pdb`

一次鞍点として妥当な TS では、虚振動（負の cm^-1）は 1 本になることが期待されます。

## 補足

- VRAM に余裕がある場合は `--hessian-calc-mode Analytical` を推奨します。
- 全オプションは `mlmm tsopt --help-advanced` と `mlmm freq --help-advanced` を参照してください。

## 次の導線

- 反応経路追跡は [irc](irc.md)、一括実行は [all](all.md) を参照してください。
