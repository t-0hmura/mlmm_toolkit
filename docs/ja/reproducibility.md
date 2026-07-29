# 再現性と決定論性

GPU 上の MLIP 推論は、ハードウェアと software stack に依存する並列演算を含むため、デフォルトでビット単位の再現性は保証されません。対象 backend、model、hardware、software stack で数値感度を評価してください。

厳密な比較が必要な場合は、`--deterministic` で決定論的 algorithm を要求し、完全な対象 stack 上で生成物を比較してください。

## `--deterministic`

`--deterministic` はすべての計算系サブコマンド（`opt`、`tsopt`、`freq`、`irc`、`scan`、`scan2d`、`scan3d`、`path-opt`、`path-search`、`all`、`sp`）で受け付けられます。これは `torch.use_deterministic_algorithms` と `index_reduce_` shim を有効化し、mlmm-toolkit が制御する演算で決定論的 algorithm を要求します。

```bash
mlmm opt -i complex.pdb --parm enzyme.parm7 -q 0 --deterministic
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -q -1 --deterministic
```

- これは**プロセス全体に適用されます**。`all` に指定すると内部の全ステージへ伝播するため、ステージごとに指定する必要はありません。
- 使用する kernel が変わるため、performance が変化する場合があります。
- PyTorch が制御する演算で決定論的実装がない場合は例外になります。backend SDK や custom operation は別途検証が必要です。
- 環境変数 `MLMM_STRICT_DETERMINISTIC=1` は、CI や直接の Python API に対する同等のエントリポイントです。

### バックエンド対応

| ML バックエンド | `--deterministic` |
|---|---|
| `uma` | deterministic mode を受理。導入済み model/SDK と対象系で検証 |
| `orb` | deterministic mode を受理。導入済み model/SDK と対象系で検証 |
| `mace` | deterministic mode を受理。導入済み model/SDK と対象系で検証 |
| `aimnet2` | **未対応 — 拒否されます**（後述） |
| `custom`（`--calc-file`） | **未対応 — 拒否されます**。指定された計算器は mlmm-toolkit の制御外です |

MM low-level 層は CPU 上で実行されます。end-to-end の厳密な比較には、topology を含む入力、software version、hardware、backend 設定を固定してください。

## 精度と再現性

`--precision fp64` は数値精度を変更しますが、GPU 実行のビット単位同一性を保証しません。`--deterministic` は決定論的 algorithm を要求しますが、完全な対象 stack で厳密な再現性を確認してください。

`--precision fp64` と（内部的に常時有効な）fp64 Hessian（`H_double`）は独立した設定項目です。`--precision fp64` を渡すと、Hessian も追加で fp64 に強制され、オプティマイザの線形代数がモデルより低い精度で警告なく実行されることがないようにします。

精度の選択は対象 backend/model と系で検証してください。backend ごとの default は [デバイス設定 & HPC セットアップ](device-hpc.md) を参照してください。

## AIMNet2 の制限

AIMNet2 はこれらの機能には対応していません:

- **`--precision fp64`** — AIMNet2 のモデル入力は上流で float32 にキャストされるため、「fp64」実行は実際には fp64 になりません。
- **`--deterministic`** — AIMNet2 の custom CUDA force kernel は `torch.use_deterministic_algorithms` の制御外なので、この flag では厳密な力の再現性を強制できません。このため option は拒否されます。

UMA、Orb、MACE は flag を受理しますが、導入済み backend/model/SDK と対象 stack で厳密な再現性を確認してください。
