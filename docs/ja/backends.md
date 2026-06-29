# MLIP Backends

mlmm-toolkit は、あらゆる ML/MM ワークフローステージ（`opt`、`scan`、`tsopt`、`freq`、`irc`、`path-search`,...）を単一の `MLMMCore` ONIOM 結合オブジェクトを通じて実行します。`MLMMCore` は ML 領域を、private な `_create_ml_backend` ファクトリ経由でバックエンドごとのアダプタ（`_UMABackend` / `_OrbBackend` / `_MACEBackend` / `_AIMNet2Backend`）にディスパッチします。このページでは、バックエンドの選択方法、バックエンドごとの kwargs、新しいバックエンドの追加方法を説明します。

## Public surface

```python
from mlmm.backends.mlmm_calc import MLMMCore, MLMMASECalculator, mlmm

# Primary entry: ML/MM ONIOM core. Takes parm7 + layered PDB + model-PDB and
# returns energy / forces / (partial) Hessian on the full system.
core = MLMMCore(
    input_pdb="layered.pdb",
    real_parm7="real.parm7",
    model_pdb="model.pdb",
    backend="uma",          # one of: "uma", "orb", "mace", "aimnet2"
    model_charge=0, model_mult=1,
    uma_model="uma-s-1p1",
    uma_precision="fp32",   # or "fp64" (full-precision base inference)
)

# ASE adapter (DMF and other ASE-based stages)
ase_calc = MLMMASECalculator(core)

# pysisyphus Calculator adapter (opt / tsopt / freq / irc / path-search stages)
pysis_calc = mlmm(
    input_pdb="layered.pdb",
    real_parm7="real.parm7",
    model_pdb="model.pdb",
    backend="uma",
    model_charge=0, model_mult=1,
)
```

内部的には、`MLMMCore.__init__` が `_create_ml_backend(backend, ...)`（`mlmm/backends/mlmm_calc.py` 内の
private なファクトリ）を呼び出して適切なアダプタをインスタンス化します。このファクトリは未知のバックエンドに対して
`ValueError` を送出します。mlmm には `'auto'` フォールバックはありません。ワークフローコードが CLI で解決されたバックエンド名を渡します。

## File map

| file | role |
|------|------|
| `mlmm/backends/__init__.py` | `apply_precision_to_calc_cfg()` — 統一された `--precision fp32\|fp64` CLI フラグを各バックエンドのネイティブ kwarg（`uma_precision` / `orb_precision` / `mace_dtype`）にルーティングします |
| `mlmm/backends/mlmm_calc.py` | `MLMMCore`（ML/MM ONIOM 結合）+ `MLMMASECalculator`（ASE）+ `mlmm`（pysisyphus Calculator）+ バックエンドごとのアダプタ（`_UMABackend`、`_OrbBackend`、`_MACEBackend`、`_AIMNet2Backend`）+ private な `_create_ml_backend` ファクトリ + FD-Hessian の組み立て + 単位変換 |
| `mlmm/backends/xtb_embedcharge_correction.py` | MM→ML 環境効果のための xTB 点電荷埋め込み補正（`--embedcharge` フラグ） |

## Per-backend characteristics

| backend | install | model identifier | precision option |
|---------|---------|------------------|------------------|
| `uma` | `pip install fairchem-core` + HF auth | `uma-s-1p1` / `uma-m-1p1` | `uma_precision="fp32" \| "fp64"` |
| `orb` | `pip install orb-models` | `orb_v3_conservative_omol` | `orb_precision="float32-high" \| "float64"`（`"float32"` も別名として受理） |
| `mace` | 専用環境: `pip uninstall -y fairchem-core && pip install mace-torch`（`e3nn` の pin が UMA と競合） | `MACE-OMOL-0` | `mace_dtype="float32" \| "float64"` |
| `aimnet2` | `pip install aimnet` | `aimnet2` | n/a |

### UMA fp64

OMol で訓練された UMA をデフォルトの fp32 から fp64 に切り替えると、TSopt + Hessian に
無視できない影響を与える場合があります。次のように有効化します:

```bash
mlmm tsopt -i ts.pdb --parm real.parm7 -q 0 -m 1 --precision fp64...
mlmm freq -i opt.pdb --parm real.parm7 -q 0 -m 1 --precision fp64...
mlmm irc -i ts.pdb --parm real.parm7 -q 0 -m 1 --precision fp64...
```

統一された `--precision` フラグは、`mlmm/backends/__init__.py` の `apply_precision_to_calc_cfg`
によって各バックエンドのネイティブ kwarg（UMA は `uma_precision`、ORB は `orb_precision`、MACE は
`mace_dtype`）へルーティングされます。

統一された `--backend-model NAME` フラグも同様に、選択中の `--backend` のモデル変種を
上書きし、`apply_backend_model_to_calc_cfg` によってバックエンドのモデル kwarg
（`uma_model` / `orb_model` / `mace_model` / `aimnet2_model`）へルーティングされます。
未指定ならバックエンドデフォルトのモデルを使用します。

または YAML 設定経由（バックエンドごとの kwarg 名）:

```yaml
calc:
 uma_precision: fp64
```

`InferenceSettings` API のために `fairchem-core ≥ 2.0` が必要です。

## Add-a-backend recipe (5 steps)

`--backend xyz` として公開する新しいバックエンド `XYZModel` を追加するには:

1. **バックエンドアダプタを作成** — `mlmm/backends/mlmm_calc.py`（あるいは大きくなる場合は
 `mlmm/backends/xyz.py` のような新規ファイル）に、`_MLBackend`（ABC）を継承する
 `_XYZBackend(_MLBackend)` を実装します。ASE 経路が必要な場合は並行して `_XYZASEBackend` も実装します。
 どちらも共通 kwargs（`charge / spin / device /
 freeze_atoms / hessian_calc_mode / return_partial_hessian / hessian_double /
 print_timing / model`）と、バックエンド固有の kwargs（`precision`、
 `default_dtype` など）を受け取る必要があります。
2. **`_MLBackend` に準拠** — 抽象メソッド
 `eval(atoms, need_grad=True) -> (E_eV, F_eV, opaque)`（エネルギーは eV、力は
 eV/Å、加えてバックエンド固有の opaque オブジェクト）、`hessian_analytical(opaque, n_atoms,
 *, dtype) -> torch.Tensor`（Hessian を eV/Å² で返す）、および
 `supports_analytical_hessian` と `device` プロパティを実装します。`_MLBackend` を継承すると、
 汎用の有限差分 `hessian_fd(...)`（バックエンドが解析的 Hessian を持たない場合に使用）を
 そのまま利用できます。
3. **`_create_ml_backend` に登録** — `mlmm/backends/mlmm_calc.py` のファクトリを拡張して、
 `backend == "xyz"` を `_XYZBackend(...)` にディスパッチします。
 `MLMMCore` が転送できるよう、新しいバックエンドの kwargs を `_create_ml_backend(...)` の
 シグネチャに追加します。
4. **統一された `--precision` フラグを配線**（任意） — バックエンドが精度の設定項目を公開する場合は、
 `mlmm/backends/__init__.py` の `_PRECISION_DISPATCH` 内の `"fp32"`
 と `"fp64"` の両方に `"xyz": (kw_name, kw_value)` エントリを追加し、
 ユーザー向けの `--precision fp32|fp64` CLI フラグが正しくルーティングされるようにします。
5. **ドキュメント化 + smoke** — このページの file map / バックエンドごとのテーブルにエントリを追加し、
 model identifier + インストールコマンドを記載し、新しいバックエンドが end-to-end で
 動作確認されるよう `tests/smoke/run.sh` に `xyz` 行を追加します。

## VRAM invariant (ML/MM-specific)

ML/MM ONIOM ジョブでは、ML バックエンドが PySCF（DFT 補正）、
parmed（parm7）、MM 力場配列と同一デバイス上に共存します。`mlmm/backends/mlmm_calc.py` 内の
方向ごとの FD-Hessian ループは、この合計メモリ使用量に収まるよう設計されています。GPU smoke ゲート全体
（`tests/smoke/run.sh`）の再実行とピーク VRAM の監視を行わない限り、**方向ごとのループを
バッチ化テンソルにリファクタリングしないでください**。さもないと全タンパク質 ML/MM の `all` ジョブが OOM します。ステージランナーは
ステージ間で `del calc` を行い、`all` ワークフローはステージ境界で `gc.collect()` を実行します。
これは public contract の一部であり、ワークフローのリファクタリング時に
削除してはなりません。

## ONIOM coupling vs raw MLIP

`mlmm/backends/mlmm_calc.py` の MLIP アダプタは、**ML 領域のみ**を評価します。
減算的 ONIOM エネルギー式（`# CHEMISTRY-RULE:1`）、リンク原子 Hessian の
B 行列射影（`# CHEMISTRY-RULE:2`）、3 層 5 パスの partial Hessian
組み立て（`# CHEMISTRY-RULE:8`）は、同じファイル内に存在します。新しい MLIP を追加する
バックエンドの作成者は ONIOM 結合を知る必要はありません。ML 領域のエネルギー / 力 / Hessian を
正しい単位で返す Calculator を公開するだけで十分です。

## See Also

- [Python API](python-api.md) — `MLMMCore` / `MLMMASECalculator` / `mlmm`（pysisyphus Calculator）の public surface。
- [Architecture](architecture.md) — 6 層ディレクトリマップ + 依存方向。
- [CONTRIBUTING](https://github.com/t-0hmura/mlmm_toolkit/blob/main/CONTRIBUTING.md) — Recipe 3.2「Add an MLIP backend」（完全なゲートサイクル参照付き）。
