# トラブルシューティング

このページでは、`mlmm` でよく遭遇するエラーと対処法をまとめます。
エラーメッセージでページ内検索すると、該当するセクションを素早く見つけることができます。
症状から先に当たりを付けたい場合は、先に [典型エラー別レシピ](recipes-common-errors.md) を見てからこのページに戻ってください。

---

## 実行前チェックリスト

長時間の計算を実行する前に、以下を確認してください。

- `mlmm -h` でヘルプが表示される
- UMA のモデルがダウンロードできる（Hugging Face のログイン/トークンが利用可能）
- AmberTools（tleap, antechamber, parmchk2）が PATH 上にある（`mm-parm` を使う場合）
- hessian_ff の C++ 拡張がビルド済み
- 酵素系ワークフローでは、入力 PDB に **水素** と **元素記号（element column）** が入っている
- 複数の PDB を与える場合、**同じ原子が同じ順序** で並んでいる（座標だけが異なる）

---

## インストール / 環境の問題

### AmberTools が見つからない

症状:
- `mm-parm` 実行時に「AmberTools (tleap, antechamber, parmchk2) not found on PATH」と表示される。

対処:
- conda で AmberTools をインストールします:

  ```bash
  conda install -c conda-forge ambertools -y
  ```

- HPC クラスターでは環境モジュールでロードする場合もあります:

  ```bash
  module load ambertools
  ```

- AmberTools なしでも、`--real-parm7` を手動で用意すれば `opt`、`tsopt`、`path-search` 等は動作します。

---

### hessian_ff のビルドエラー

症状:
- `import hessian_ff` で ImportError が出る
- `make` が失敗する（コンパイラが見つからない、ヘッダが見つからない等）

対処:
- C++ コンパイラ（gcc 8+ または clang）がインストールされていることを確認:

  ```bash
  g++ --version
  ```

- PyTorch のヘッダが利用可能であることを確認:

  ```bash
  python -c "import torch; print(torch.utils.cmake_prefix_path)"
  ```

- ビルドを再実行:

  ```bash
  cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make clean && make
  ```

---

### UMA のダウンロード/認証エラー（Hugging Face）

症状:
- モデルをダウンロードできない、認証が必要、といったエラー。

対処:
- 環境/マシンごとに一度ログインします:

  ```bash
  huggingface-cli login
  ```

- HPC では、計算ノードから HF キャッシュ（ホームディレクトリ等）に書き込み可能か確認してください。

---

### CUDA / PyTorch の不整合

症状:
- GPU があるのに `torch.cuda.is_available()` が False
- import 時に CUDA runtime error が出る

対処:
- クラスタの CUDA バージョンに適合する PyTorch をインストールしてください:
- GPU が見えているか確認します:

  ```bash
  nvidia-smi
  python -c "import torch; print(torch.version.cuda, torch.cuda.is_available())"
  ```

---

### CUDA メモリ不足（OOM）

症状:
- `torch.cuda.OutOfMemoryError` が発生
- 「CUDA out of memory」メッセージ

対処:
- **ML ヘシアンモードを有限差分に切り替え**: `--hessian-calc-mode FiniteDifference`（解析ヘシアンは VRAM を多く消費）
- **ポケットサイズを縮小**: `--radius` を小さくしてクラスターモデルの原子数を減らす
- **3 層 + Hessian 対象制御を活用**: `hess_cutoff` と `movable_cutoff` を設定し、ヘシアン計算対象の原子数を制限:
  ```yaml
  calc:
    hess_cutoff: 3.6
    movable_cutoff: 8.0
  ```
- **`define-layer` で事前に層を定義**し、`use_bfactor_layers: true` で読み取る
- **GPU メモリが大きいカードに変更**（16 GB 以上推奨、500 原子以上の系）

---

### DMF モードが動かない（cyipopt がない）

DMF（`--mep-mode dmf`）を使うときに IPOPT/cyipopt の import エラーが出る場合:

対処:
- `mlmm_toolkit` をインストールする前に、conda-forge から `cyipopt` をインストールしてください:

  ```bash
  conda install -c conda-forge cyipopt
  ```

---

### 図のエクスポートが失敗する（Chrome がない）

Plotly/Chrome 系のエラーで静的画像が出ない場合:

対処:
- headless Chrome をインストールしてください:

  ```bash
  plotly_get_chrome -y
  ```

---

## 入力 / 抽出の問題

### 「Element symbols are missing ... add-elem-info を実行してください」

典型的なメッセージ:

```text
Element symbols are missing in '...'.
Please run `mlmm add-elem-info -i ...` to populate element columns before running extract.
```

対処:
- 次を実行して element 列（元素記号列）を補完します:

  ```bash
  mlmm add-elem-info -i input.pdb -o input_with_elem.pdb
  ```

- その後、`extract` / `all` を補完後の PDB で再実行します。

---

### 「Atom count mismatch」「Atom order mismatch」

典型的なメッセージ:

```text
[multi] Atom count mismatch between input #1 and input #2: ...
[multi] Atom order mismatch between input #1 and input #2.
```

対処:
- **すべて** の構造を同じ前処理ワークフロー（同じプロトン化ツール、同じ設定）で作り直します。
- MD 由来なら、同一のトポロジ/軌跡からフレーム抽出する方が安全です。

---

### 「ポケットが空っぽ / 必要な残基が落ちる」

症状:
- 抽出ポケットが想定より小さい
- 触媒残基が含まれない

対処の例:
- `--radius` を増やす（例: 2.6 → 3.5 Å）
- `--selected-resn` で残基を強制包含する（例: `--selected-resn 'A:123,B:456'`）
- 主鎖削除が強すぎる場合は `--no-exclude-backbone` を試す

---

## parm7 / MM の問題

### parm7 と入力 PDB の原子数不整合

症状:
- 「Coordinate shape mismatch for ... got (N, 3), expected (M, 3)」
- 「Atom count mismatch: parm7 has M atoms, but coordinate file has N atoms」

対処:
- parm7 を生成した PDB と、計算に使う PDB が同一の原子セットであることを確認
- `mm-parm` で parm7 を再生成する場合、出力される PDB（`<prefix>.pdb`）を計算入力として使用
- tleap が水素を追加/削除している可能性があるため、`mm-parm` 出力の PDB を入力に使う

---

### parm7 の元素順序が PDB と一致しない

症状:
- `oniom-gaussian` / `oniom-orca` で「Element sequence mismatch at atom index ...」

対処:
- `--no-element-check` で元素チェックを無効化（結果を手動で検証すること）
- 正しい対処は、parm7 生成時と同じ PDB を `-i` に指定すること

---

### antechamber が失敗する

症状:
- `mm-parm` 実行中に「antechamber failed (see log)」

対処:
- `--keep-temp` を付けて再実行し、`<resname>.antechamber.log` を確認:

  ```bash
  mlmm mm-parm -i input.pdb --ligand-charge 'LIG:-1' --keep-temp
  ```

- リガンドの電荷と多重度が正しいか確認（`--ligand-charge`, `--ligand-mult`）
- 水素が正しく付加されているか確認
- 入力 PDB の TER レコードが適切か確認

---

## 電荷 / スピンの問題

### 「電荷が必須」系のエラー

対処:
- 電荷を明示する:

  ```bash
  mlmm path-search -i R.pdb P.pdb --real-parm7 real.parm7 --model-pdb model.pdb -q 0
  ```

- あるいは（抽出ありの場合）残基名ごとの電荷マッピングを与える:

  ```bash
  mlmm -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
  ```

---

## 計算 / 収束の問題

### TS 最適化が収束しない

症状:
- TS 最適化が多くのサイクルを回しても収束しない
- 最適化後も複数の虚振動数が残る

対処の例:
- オプティマイザモードを切り替える: `--opt-mode light` (Dimer) または `--opt-mode heavy` (RS-I-RFO)
- 余分な虚モードのフラット化を有効にする: `--flatten-imag-mode`
- 最大サイクル数を増やす: `--tsopt-max-cycles 20000`
- `hess_cutoff` を調整して、ヘシアン計算に含む原子の範囲を広げる

---

### IRC が正常に終了しない

症状:
- IRC が明確な極小に到達する前に停止
- エネルギーが振動したり勾配が大きいまま

対処の例:
- ステップサイズを減らす: `--step-size 0.05`（デフォルトは 0.10）
- 最大サイクル数を増やす: `--max-cycles 200`
- IRC 実行前に TS 候補で虚数振動数が 1 本であることを確認

---

### MEP 探索（GSM/DMF）が失敗または予期しない結果

症状:
- 経路探索が有効な MEP なしで終了
- 結合変化が正しく検出されない

対処の例:
- `--max-nodes` を増やす（複雑な反応には 15 や 20 など）
- 端点の事前最適化を有効にする: `--preopt`
- 別の MEP 手法を試す: `--mep-mode dmf`（GSM が失敗した場合）またはその逆
- YAML で結合検出パラメータを調整（`bond.bond_factor`、`bond.delta_fraction`）

---

## パフォーマンス / 安定性のヒント

- **VRAM 不足**: `--radius` を減らす、`hess_cutoff`/`movable_cutoff` を設定する
- **解析ヘシアンが遅いまたは OOM**: `--hessian-calc-mode FiniteDifference` を使用。`Analytical` は十分な VRAM がある場合のみ推奨
- **MM ヘシアン計算が遅い**: `hess_cutoff` を設定して Hessian-MM 原子数を制限する
- **大規模系（1000 原子以上）**: `define-layer` で 3 層を適切に定義し、Frozen 層で計算コストを削減
- **ML と MM の並列実行**: デフォルトで ML（GPU）と MM（CPU）は並列実行されます。`mm_threads` で CPU スレッド数を調整可能

---

## 不具合報告のときに添えると助かる情報

- 実行したコマンド（コピペ可能な形）
- `summary.log`（またはコンソール出力）
- 再現する最小入力（可能なら）
- OS / Python / CUDA / PyTorch バージョン
- AmberTools / hessian_ff のバージョン
