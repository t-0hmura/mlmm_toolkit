# トラブルシューティング

このページでは、`mlmm` でよく遭遇するエラーと対処法をまとめます。
エラーメッセージでページ内検索すると、該当するセクションを素早く見つけることができます。
症状から先に当たりを付けたい場合は、先に [典型エラー別レシピ](recipes-common-errors.md) を見てからこのページに戻ってください。

---

## 実行前チェックリスト

長時間の計算を実行する前に、以下を確認してください。

- `mlmm -h` でヘルプが表示される
- MLIP モデルの重みがダウンロードできる（デフォルトの UMA バックエンドの場合、Hugging Face のログイン/トークンが必要。他のバックエンドは別のソースからダウンロードする場合がある）
- 酵素系ワークフローでは、入力 PDB に **水素** と **元素記号（element column）** が入っている
- 複数の PDB を与える場合、**同じ原子が同じ順序** で並んでいる（座標だけが異なる）
- **AmberTools** が conda チャンネル（またはソースビルド）で正しくインストールされ、`tleap` が利用可能（`mm-parm` を使う場合）
- hessian_ff の C++ ネイティブ拡張が正しくビルド済み（自動ビルドが失敗した場合は `cd hessian_ff/native && make` を実行）

---

## 入力 / 抽出の問題

### 「Element symbols are missing... add-elem-info を実行してください」

典型的なメッセージ:

```text
Element symbols are missing in '...'.
Please run `mlmm add-elem-info -i...` to populate element columns before running extract.
```

対処:
- 次を実行して element 列（元素記号列）を補完します:

 ```bash
 mlmm add-elem-info -i input.pdb -o input_with_elem.pdb
 ```

- その後、`extract` / `all` を補完後の PDB で再実行します。

なぜ発生するか:
- PDB によっては元素カラムが一貫して埋められていないことがあります。`extract` は正確な原子型判定のために元素記号を必要とします。

---

### 「Atom count mismatch」「Atom order mismatch」

典型的なメッセージ:

```text
[multi] Atom count mismatch between input #1 and input #2:...
[multi] Atom order mismatch between input #1 and input #2.
```

対処:
- **すべて** の構造を同じ前処理ワークフロー（同じプロトン化ツール、同じ設定）で作り直します。
- 水素付加を行う場合は、全フレームで一貫した原子順序になる方法を使用してください。

ヒント:
- MD 由来のアンサンブルでは、異なるツールで作った PDB を混在させるのではなく、**同一のトポロジー/軌跡** からフレーム抽出する方が安全です。

代替手段:
- 複数構造の入力を揃えることが困難な場合は、**単一構造スキャンワークフロー** を使用してください。1 つの PDB と `--scan-lists` を指定し、距離スキャンで端点を生成します。

---

### ML 領域に重要な残基が含まれない

症状:
- 抽出ポケットが想定より小さい
- 触媒残基が含まれない

対処の例:
- `--radius` を増やす（例: 2.6 → 3.5 Angstrom）
- `--selected-resn` で残基を強制包含する（例: `--selected-resn 'A:123,B:456'`）
- PyMOL 等の分子ビューアで活性部位の原子を選択・エクスポートし、ML 領域の PDB を手動で作成することもできます。`--model-pdb` でこの PDB を指定してください。

---

### エネルギー・障壁の計算値が不正確

症状:
- 計算されたエネルギーや反応障壁が不合理に見える
- モデルサイズを大きくすると結果が大幅に変わる

対処:
- 抽出されたポケットが小さすぎると、エネルギーや障壁の計算値が不正確になることがあります。抽出半径を大きくする（例: `-r 4.0` 以上）ことで、タンパク質環境をより多く含めて精度を改善できます:

 ```bash
 mlmm extract -i complex.pdb -c 'SUB' -o pocket.pdb -r 4.0
 ```

---

### 非標準残基が正しく切断されない

抽出されたポケットに非標準の3文字コードを持つ修飾アミノ酸残基（リン酸化セリン、メチル化リシンなど）が含まれている場合、デフォルトでは主鎖切断やリンク水素付加が適用されません。`--modified-residue` で登録してください:

```bash
mlmm extract -i complex.pdb -c PRE --modified-residue "SEP,TPO,MLY" -o pocket.pdb
```

同じフラグは `all` コマンドでも使用可能で、抽出ステージに転送されます。

`--modified-residue` で対応できない場合は、リンク水素を含むポケットモデルを手動で構築し、`--parm` と `--model-pdb` を使って下流コマンドに直接渡してください。

---

## 電荷 / スピンの問題

### 電荷解決の問題

計算系サブコマンドでは `-q/--charge` と `-m/--multiplicity` を明示する必要があります。
`all` では電荷は `-q/--charge` 上書き -> 抽出サマリー -> （抽出スキップ時）`--ligand-charge` フォールバック の順で解決されます。

対処:
- 電荷と多重度を明示的に指定します:

 ```bash
 mlmm path-search -i R.pdb P.pdb --parm real.parm7 --model-pdb model.pdb -q 0 -m 1
 ```

- あるいは、抽出経由で自動導出させるため `all` を使い、残基名ごとの電荷マッピングを与えます:

 ```bash
 mlmm -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
 ```

---

## AmberTools / mm-parm の問題

### tleap が見つからない

典型的なメッセージ:

```text
FileNotFoundError: tleap not found on PATH
```

または

```text
mm-parm requires AmberTools (tleap, antechamber, parmchk2).
```

対処:
- conda で AmberTools をインストールします:

 ```bash
 conda install -c conda-forge ambertools -y
 ```

- ソースからビルド（<https://ambermd.org/AmberTools.php>）するか、HPC クラスターでは環境モジュールでロードします:

 ```bash
 module load ambertools
 ```

- 利用可能か確認します:

 ```bash
 which tleap
 which antechamber
 which parmchk2
 ```

- AmberTools なしでも、`--parm` を手動で用意すれば `opt`、`tsopt`、`path-search` 等は動作します。

---

### antechamber が失敗する

症状:
- `mm-parm` 実行中にリガンドパラメータ化が失敗する
- 原子型割り当てや電荷計算に関するエラーが出る

対処の例:
- リガンドの PDB で元素記号と結合構造が正しいか確認する
- `--ligand-charge` が正しく指定されていることを確認: `-l 'GPP:-3,SAM:1'`
- `--keep-temp` を付けて再実行し、`<resname>.antechamber.log` を確認する:

 ```bash
 mlmm mm-parm -i input.pdb -l 'LIG:-1' --keep-temp
 ```

- 水素が正しく付加されているか、TER レコードが適切か確認する
- 非一重項リガンドには `--ligand-mult` を指定する（例: `--ligand-mult 'HEM:1,NO:2'`）。デフォルトのスピン多重度は 1（一重項）
- 抽出されたリガンド PDB に対して antechamber を手動実行して原因を切り分ける:

 ```bash
 antechamber -i ligand.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc -nc -3 -at gaff2
 ```

- より高精度の部分電荷が必要な場合は、HF/6-31G* 計算から RESP 電荷を算出し、カスタム `frcmod`/`lib` ファイルを用意することを検討してください（AM1-BCC に頼らない方法）。

---

### parm7/rst7 の不整合エラー

典型的なメッセージ:

```text
Atom count in parm7 (...) does not match input PDB (...)
```

または

```text
RuntimeError: parm7 topology does not match the input structure
```

または

```text
Coordinate shape mismatch for... got (N, 3), expected (M, 3)
```

対処:
- parm7 を生成した PDB と、計算に使う PDB が同一の原子セット（同じ順序）であることを確認します。
- `mm-parm` で parm7 を再生成してください。
- `mm-parm` 実行後に PDB の原子を編集・並べ替え **しない** でください。
- tleap が水素を追加/削除している可能性があるため、`mm-parm` 出力の PDB（`<prefix>.pdb`）を計算入力として使用してください。

---

### parm7 の元素順序が PDB と一致しない

症状:
- `oniom-export` で「Element sequence mismatch at atom index...」

対処:
- `--no-element-check` で元素チェックを無効化（結果を手動で検証すること）
- 正しい対処は、parm7 生成時と同じ PDB を `-i` に指定することです。

---

## hessian_ff ビルドの問題

### ビルドが失敗する（"make" エラー）

典型的な症状:
- `hessian_ff/native/` での `make` がコンパイルエラーを出す
- `ImportError: cannot import name 'ForceFieldTorch' from 'hessian_ff'`
- `RuntimeError: hessian_ff build attempts failed: ...`

対処の例:
- C++ コンパイラ（g++ >= 9）がインストールされていることを確認:

 ```bash
 g++ --version
 ```

- PyTorch のヘッダが利用可能であることを確認:

 ```bash
 python -c "import torch; print(torch.utils.cmake_prefix_path)"
 ```

- HPC ではコンパイラモジュールをロード:

 ```bash
 module load gcc/11
 ```

- クリーンしてリビルド:

 ```bash
 conda install -c conda-forge ninja -y
 cd hessian_ff/native && make clean && make
 ```

---

### hessian_ff の import エラー

典型的なメッセージ:

```text
ImportError: cannot import name 'ForceFieldTorch' from 'hessian_ff'
```

または:

```text
RuntimeError: hessian_ff build attempts failed: ...
To rebuild hessian_ff native extensions in this environment:
  conda install -c conda-forge ninja -y
  cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make clean && make
```

対処:
- C++ ネイティブ拡張を先にビルドする必要があります:

 ```bash
 cd hessian_ff/native && make
 ```

- `hessian_ff` パッケージが Python パス上にあることを確認してください（`pip install -e .` でインストールしていれば問題ありません）。

---

## B-factor レイヤー割り当ての問題

### レイヤー割り当てが想定と異なる

症状:
- 原子が想定外のレイヤーに割り当てられる
- ML 領域が小さすぎる、または大きすぎる

対処の例:
- B-factor エンコーディング: ML = 0.0、Movable-MM = 10.0、Frozen-MM = 20.0。
- レイヤーが割り当てられた PDB を分子ビューアで可視化（B-factor で色分け）する
- `--model-pdb` が正しく ML 領域の原子を定義しているか確認する
- `define-layer` の距離カットオフを調整する:
 - `--radius-freeze`（デフォルト 8.0 Angstrom）: Movable-MM/Frozen の境界を制御
- 必要に応じて、計算オプション（`hess_cutoff`, `hess_mm_atoms`）でヘシアン対象 MM を別途制御する
- YAML で `use_bfactor_layers: true` を使う場合、B-factor 値が期待されるエンコーディング（0.0, 10.0, 20.0; 許容差 1.0）と一致するか確認する

---

### B-factor 値が認識されない

典型的な症状:
- 計算機がすべての原子を Frozen または ML として扱う
- B-factor 値が {0.0, 10.0, 20.0} のいずれでもない

対処:
- `define-layer` を再実行して正しい B-factor エンコーディングを確保してください。
- 許容差 1.0 が適用されます: 0/10/20 に近い B-factor が ML/Movable/Frozen にマッピングされます。
- B-factor を任意の値に手動編集しないでください。

---

### `--detect-layer` が想定どおりに働かない

症状:
- B-factor からのレイヤー自動判定で、ML / Movable / Frozen の分割が想定と異なる
- `--detect-layer` 有効時に `--model-pdb` なしで実行すると失敗する

対処の例:
- 入力が PDB（または `--ref-pdb` 付き XYZ）であることを確認する
- `define-layer` で B-factor を明示的に再付与し、生成された PDB を使う
- 距離ベース制御を使う場合は `hess_cutoff` / `movable_cutoff` を指定し、必要なら `--no-detect-layer` に切り替える
- `--movable-cutoff` を与えると `--detect-layer` が無効化される点に注意する

---

## インストール / 環境の問題

### MLIP モデルのダウンロードエラー

症状:
- MLIP モデルの重みをダウンロードできない、認証が必要、といったエラー。デフォルトの UMA バックエンドでは、Hugging Face のログイン/トークンが不足していることが原因です。

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
- クラスタの CUDA バージョンに適合する PyTorch をインストールしてください。
- GPU が見えているか確認します:

 ```bash
 nvidia-smi
 python -c "import torch; print(torch.version.cuda, torch.cuda.is_available())"
 ```

---

### DMF モードが動かない（cyipopt がない）

DMF（`--mep-mode dmf`）を使うときに IPOPT/cyipopt の import エラーが出る場合:

対処:
- `mlmm` をインストールする前に、conda-forge から `cyipopt` をインストールしてください:

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

## 計算 / 収束の問題

### CUDA メモリ不足（OOM）

症状:
- `torch.cuda.OutOfMemoryError: CUDA out of memory`
- 「CUDA out of memory」メッセージ
- ヘシアン計算中にシステムがハングまたはクラッシュする

ML/MM 系は純粋な MLIP 計算よりも一般的に大きいため、VRAM の負荷が高くなります。

対処の例（優先度順）:
- **Frozen-MM 層を確認**: `define-layer` で Frozen-MM 原子（B=20.0）が正しく割り当てられているか確認する。Frozen-MM 領域が小さすぎると、Movable-MM 領域（ひいてはヘシアン）が不必要に大きくなる。`--radius-freeze` を小さくして Frozen 領域を拡大する。
- **ML 領域サイズを縮小**: `extract` の `--radius` を小さくするか、`--model-pdb` で手動定義した小さい ML 領域 PDB を指定する。
- **有限差分 ML ヘシアンを使用**: `--hessian-calc-mode FiniteDifference`（VRAM 消費が少ないが低速）。
- **`define-layer` で事前に層を定義** し、`use_bfactor_layers: true` で読み取る。
- **GPU メモリが大きいカードに変更**: 500 原子以上の ML 領域には 24 GB 以上推奨、1000 原子以上には 48 GB 以上推奨。

---

### TS 最適化が収束しない

症状:
- TS 最適化が多くのサイクルを回しても収束しない
- 最適化後も複数の虚振動数が残る

対処の例:
- オプティマイザモードを切り替える: `--opt-mode grad` (Dimer) または `--opt-mode hess` (RS-I-RFO)
- 余分な虚モードのフラット化を有効にする: `--flatten`
- 最大サイクル数を増やす: `--max-cycles 20000`
- より厳しい収束閾値を使う: `--thresh baker` または `--thresh gau_tight`
- `hess_cutoff` を調整して、ヘシアン計算に含む原子の範囲を広げる

---

### IRC が正常に終了しない

症状:
- IRC が明確な極小に到達する前に停止
- エネルギーが振動したり勾配が大きいまま

対処の例:
- ステップサイズを減らす: `--step-size 0.05`（デフォルトは 0.10）
- 最大サイクル数を増やす: `--max-cycles 200`
- IRC 実行前に TS 候補で虚振動数が 1 本であることを確認

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

- **VRAM 不足**: ML 領域サイズを縮小、ヘシアン対象 MM 領域を縮小、ノード数を削減（`--max-nodes`）、または軽量なオプティマイザ設定（`--opt-mode grad`）を使用。
- **解析ヘシアンが遅いまたは OOM**: `--hessian-calc-mode FiniteDifference` を使用。`Analytical` は十分な VRAM がある場合のみ推奨（ML 原子 300 以上には 24 GB 以上推奨）
- **MM ヘシアン**: `mm_fd: true`（デフォルト）は MM ヘシアンに有限差分を使用。解析 MM ヘシアン（`mm_fd: false`）は小規模系では高速だがメモリ消費が増える場合がある
- **MM ヘシアン計算が遅い**: `hess_cutoff` を設定して Hessian-MM 原子数を制限する
- **大規模系（2000 原子以上）**: `define-layer` で `--radius-freeze` を小さくして Frozen 層（B=20）を拡大し、可動自由度数を削減する
- **マルチ GPU**: ML を 1 つの GPU（`ml_cuda_idx: 0`）、MM を別の GPU（`mm_device: cuda`, `mm_cuda_idx: 1`）に配置可能
- **ML と MM の並列実行**: デフォルトで ML（GPU）と MM（CPU）は並列実行されます。`mm_threads` で CPU スレッド数を調整可能

---

## バックエンド固有の問題

### --backend orb/mace/aimnet2 で ImportError が出る

**症状:** `ImportError: orb-models is required for the ORB backend` のようなエラーが表示される

**対処:** 使用するバックエンドに対応するオプション依存パッケージをインストールします:
```bash
pip install "mlmm-toolkit[orb]"      # ORB バックエンド
pip install "mlmm-toolkit[aimnet]"  # AIMNet2 バックエンド
pip install --no-deps mace-torch      # MACE バックエンド
```

---

### 非 UMA バックエンドで CUDA メモリ不足になる

**症状:** ORB、MACE、AIMNet2 の使用時に `RuntimeError: CUDA out of memory` が発生する

**対処:** 非 UMA バックエンドは有限差分ヘシアンを使用するため、より多くの VRAM を消費します。以下の方法を試してください:
- `--hessian-calc-mode FiniteDifference` を明示的に指定し、`hess_cutoff` を小さめに設定する
- YAML で `ml_device: cpu` を指定する（遅くなるが VRAM 制限を回避できる）

---

### --embedcharge 使用時に xTB が見つからない

**症状:** `XTBEmbedError: xTB command not found` が表示される

**対処:** xTB をインストールし、`$PATH` 上にあることを確認します:
```bash
conda install -c conda-forge xtb
```

---

## 不具合報告のときに添えると助かる情報

- 実行したコマンド（コピペ可能な形）
- `summary.log`（またはコンソール出力）
- 再現する最小入力（可能なら）
- OS / Python / CUDA / PyTorch バージョン
- AmberTools / hessian_ff のバージョン
