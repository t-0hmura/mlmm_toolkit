# クイックスタート: `mlmm tsopt`

## 目的

TS 候補を最適化し、一次鞍点（first-order saddle point）であることを確認します。

## 事前に必要なファイル

- TS 候補構造: `ts_guess.pdb`
- MM トポロジー: `real.parm7`
- ML 領域定義: `ml_region.pdb`

## 1. TS 最適化

```bash
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_tsopt
```

数値最適化が収束すると、`--skip-final-freq` 指定時を除き、`tsopt` が終端の Hessian と虚振動を確認します。コンソール出力で次の行を確認してください。

```
[Imaginary modes] n=1 ([-593.1])
```

## 出力の検証

- `result_tsopt/final_geometry.pdb` — 最適化済み TS 構造
- `result_tsopt/vib/` — 虚振動モード（変位ベクトル）の軌跡（`imag_*_trj.xyz`, `.pdb`）
- ターミナル出力: 一次鞍点の認定には虚振動が **n=1** であることが必要です。モード変位と IRC 接続性も確認してください。虚振動が複数残る場合は `--flatten` の適用を検討してください

## 2.（任意）個別の振動解析

全振動モードの一覧や熱化学補正（零点エネルギー (ZPE)、ギブズ自由エネルギーなど; `all` コマンドの `--thermo` に相当）が必要な場合は、別途 `freq` を実行してください。虚振動数の確認だけであれば、上記の `tsopt` の出力で十分です。

```bash
mlmm freq -i ./result_tsopt/final_geometry.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_freq
```

## 補足

- 対象系の予備計算で `Analytical` と `FiniteDifference` の実行時間とメモリ使用量を比較してください。
- 別の MLIP バックエンドを使用するには `-b orb`（または `mace`、`aimnet2`）を追加します。デフォルトは `uma` です。
- 全オプションは `mlmm tsopt --help-advanced` と `mlmm freq --help-advanced` を参照してください。

## 次の導線

- 反応経路追跡は [irc](irc.md)、一括実行は [all](all.md) を参照してください。
