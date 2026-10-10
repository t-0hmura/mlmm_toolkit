# トラブルシューティング

症状を早見表で探し、示された節で対処を読んでください。

(ja-troubleshooting-quick-table)=
## 早見表

| 症状 | 最初にやること | 詳細（節） |
| --- | --- | --- |
| **入力 / 抽出** | | |
| 元素欄が空で `extract` が止まる（`Element symbols are missing in '...'`）。`all` は空の元素欄を自分で埋め、割り当てられない原子が残ると止まる | 元の PDB に `add-elem-info` を適用してください | {ref}`入力 / 抽出の問題 <ja-input--extraction>` |
| `[multi] Atom count mismatch` / `Coordinate shape mismatch` / `Element sequence mismatch` | 同じ前処理ツール・同じ設定で全 PDB を作り直し、計算に使う構造からトポロジーを作り直してください。`mm-parm` の後に原子を並べ替えないでください | {ref}`入力 / 抽出の問題 <ja-input--extraction>`、{ref}`AmberTools / mm-parm の問題 <ja-ambertools--mm-parm>` |
| **電荷 / スピン** | | |
| `ML-region charge is unresolved` / `[all] ML-region charge could not be resolved` | `-q/--charge` または `-l/--ligand-charge` を明示してください | {ref}`電荷 / スピンの問題 <ja-charge--spin>` |
| 計算は通るが状態やエネルギーが不自然 | ML 領域の電荷と多重度を見直してください | {ref}`電荷 / スピンの問題 <ja-charge--spin>` |
| **計算 / 収束** | | |
| 実行時に CUDA のメモリ不足（`torch.cuda.OutOfMemoryError`） | 固定 MM の層を確かめる、ML 領域を小さくする（`--radius`）、Hessian の範囲を絞る（`--hessian-cutoff`）、`Analytical` を選んでいたらデフォルトの `FiniteDifference` に戻す、VRAM の大きい GPU に移る | {ref}`CUDA メモリ不足 <ja-cuda-oom>` |
| TS 最適化が収束しない（`TS optimization did not converge`）、または収束後の n_imag が 1 でない | まず TS 候補を確かめ、次にオプティマイザを切り替えてください（`tsopt --opt-mode` / `all --opt-mode-post`）。n_imag ≥ 2 なら `--flatten` を付けます | {ref}`TS 最適化 <ja-troubleshooting-ts>`、{ref}`TS が取れないとき <ja-ts-search-fails>` |
| IRC が正常に終了しない | まず最適化後の端点を確かめ、次にステップを小さくしてください：`irc --step-size` または `all --irc-step-size` | {ref}`IRC <ja-troubleshooting-irc>` |
| エネルギーが平坦なのに最適化が終わらない（MLIP のノイズフロアの可能性） | `--max-cycles` に任せるか、`--stop-plateau` で早期停止を有効にしてください。止まるのが早すぎる・遅すぎるときは `--stop-plateau-thresh` / `--stop-plateau-window` を調整します | {ref}`プラトーでの停止 <ja-optimizer-stalls-with-flat-energy--forces-just-above-threshold-mlip-force-noise-floor>` |
| **インストール / 環境** | | |
| UMA モデルで 401 / 403 / アクセス制限付きリポジトリのエラー（`huggingface_hub.errors.GatedRepoError`） | `hf auth login` でログインし、UMA モデルのライセンスに同意してください | {ref}`インストール / 環境の問題 <ja-installation-environment-problems>` |
| `orb-models is required for the ORB backend`（AIMNet2 / MACE も同様） | バックエンドの追加パッケージを入れてください：`pip install "mlmm-toolkit[orb]"` または `"mlmm-toolkit[aimnet]"`。MACE は別の環境に入れます | {ref}`バックエンド固有の問題 <ja-troubleshooting-backends>` |
| `mm-parm` が実行できない（`AmberTools preflight failed`。`tleap` / `antechamber` / `parmchk2` が無い） | 先に AmberTools を使えるようにしてください | {ref}`AmberTools / mm-parm の問題 <ja-ambertools--mm-parm>` |
| `hessian_ff` のビルドや import のエラー（`hessian_ff build attempts failed`） | C++20 のコンパイラを確かめ、ネイティブ拡張を作り直してください | {ref}`hessian_ff ビルドの問題 <ja-hessian_ff-build--import>` |
| DMF モードの import エラー（`DMF mode (--mep-mode dmf) requires ase, cyipopt, and pydmf>=1.2`） | `cyipopt`（conda-forge）と `pydmf[torch]>=1.2`（PyPI）を入れてください | {ref}`インストール / 環境の問題 <ja-installation-environment-problems>` |
| CUDA / GPU の実行時エラー | GPU、PyTorch のビルド、ドライバをまとめて確かめてください | {ref}`インストール / 環境の問題 <ja-installation-environment-problems>` |
| 図の出力に失敗する | `plotly_get_chrome -y` でヘッドレス Chrome を入れてください | {ref}`インストール / 環境の問題 <ja-installation-environment-problems>` |

## 実行前チェックリスト

長い計算を回す前に、次を確かめてください。

- `mlmm -h` でヘルプが表示される。
- デフォルトの UMA バックエンドのために、このマシンで Hugging Face にログインできている。
- 入力の PDB/mmCIF に **水素** と **元素記号** が入っている。
- 複数の PDB を与える場合、**同じ原子が同じ順序** で並んでいる。
- `tleap`、`antechamber`、`parmchk2` が `$PATH` にある。
- `hessian_ff` の C++ 拡張が初回の使用時にビルドできる。{ref}`hessian_ff ビルドの問題 <ja-hessian_ff-build--import>` を参照してください。

---

(ja-input--extraction)=
## 入力 / 抽出の問題

### `Element symbols are missing in '...'`

- **症状**：`extract` が ``Element symbols are missing in '...'. For PDB input, run `mlmm add-elem-info -i ... --overwrite`, or write a fixed PDB with `-o` and pass that file to extract; ...`` で止まる。`all` は抽出の前に空の元素欄を自分で埋め、割り当てられない原子が残ると同じメッセージで止まる。
- **原因**：PDB の元素欄（77–78 列）が空のことが多く、`extract` は原子の種類を決めるために元素記号を使います。mmCIF の入力では `_atom_site.type_symbol` が必要です。
- **対処**：`add-elem-info` で元素欄を埋め、新しいファイルで再実行してください。`[add-elem-info] WARNING: Could not confidently assign N atoms; left unchanged.` の後に並んだ原子は、元素記号を手で 77–78 列に右詰めで書いてください。

  ```bash
  mlmm add-elem-info -i input.pdb -o input_with_elem.pdb
  ```

### `[multi] Atom count mismatch` / `[multi] Atom order mismatch`

- **症状**：複数の入力を与えた実行が `[multi] Atom count mismatch between input #1 and input #2: ...` や `[multi] Atom order mismatch between input #1 and input #2.` で止まる。
- **原因**：構造ごとに別のツールや設定で前処理したか、プロトン化をやり直した後に原子の順序が変わりました。
- **対処**：**すべて** の構造を、同じプロトン化ツール・同じ設定で作り直してください。MD のスナップショットなら、同じトポロジーと軌跡からフレームを取り出します。複数の入力をそろえるのが難しい場合は、1 つの PDB から [`--scan-lists`](quickstart-scan.md) で経路を作れます。

### ML 領域が小さい・触媒残基が入らない

- **症状**：切り出した ML 領域が想定より小さい、または触媒残基が含まれない。
- **原因**：この部位には半径（`-r/--radius`、デフォルト 2.6 Å）が小さすぎます。
- **対処**：`--radius` を大きくするか（例：2.6 → 3.5 Å）、`--selected-resn 'A:TYR:44'` で残基を足してください。この残基からは距離の探索を始めません。`-c` に足すと、`-r` が 0 より大きければ残基が丸ごと残ります。詳しくは {ref}`モデルを広げる <ja-model-setup-larger>` を参照してください。指定できる形は {ref}`残基の指定 <ja-selected-resn-takes-ids>` にあります。chain の欄が空の PDB では、`'44'` のように名前か番号だけを使います。ML 領域の原子を自分で選び、その PDB を [`--model-pdb`](model-setup.md#自分で組んだモデルを使う) で渡すこともできます。

### エネルギーや障壁が ML 領域の大きさで変わる

[ML 領域が小さい・触媒残基が入らない](#ml-領域が小さい触媒残基が入らない) のとおりに ML 領域を広げ、結果が ML 領域の大きさと境界の位置でどう変わるかを確かめてください。

### 修飾残基が切断されない

- **症状**：`extract` が `[extract] WARNING: Residue ... may be an amino acid (has N, CA, C, O) but is not recognized as a standard residue name. Backbone truncation was not applied. ...` を出し、その残基の主鎖が切られない。
- **原因**：主鎖の切断とリンク水素の付加には、残基の表への登録が必要です。SEP、TPO、MLY などは登録済みです。
- **対処**：`--modified-residue "HD1:0"` のように、残基を整数の電荷とともに登録してください（`mlmm all` でも使えます）。名前だけで書くときの決まりは [extract](extract.md#5-非標準残基--modified-residue) にあります。主鎖のトポロジーが特殊な場合は、ML 領域を手で組み、`--parm7` と `--model-pdb` で下流のコマンドに直接渡してください。

---

(ja-charge--spin)=
## 電荷 / スピンの問題

`-q` は系全体ではなく ML 領域の電荷です。`-l/--ligand-charge` の各残基名が構造にあるかを確かめてください。決まりは {ref}`電荷の指定 <ja-charge-specification>` にあります。

### `ML-region charge is unresolved` / `ML-region charge could not be resolved`

- **症状**：個別のコマンドが `ML-region charge is unresolved. Provide -q/--charge or --ligand-charge.` で、`all` が `[all] ML-region charge could not be resolved. Provide -q/--charge, --ligand-charge, or calc.model_charge in YAML.` で止まる。
- **原因**：`-q/--charge` を省くと、電荷は ML 領域にある標準残基・イオン・`-l/--ligand-charge` の値の合計から、次に YAML の `calc.model_charge` から決まります。そのどれからも決まらなかったため止まりました。`--model-indices` では `-l` から電荷を出せません。
- **対処**：電荷と多重度を明示するか、抽出ありの場合は残基ごとの電荷を与えてください。導いた電荷は端末に `Total active site model charge` として出ます。

  ```bash
  mlmm path-search -i R.pdb P.pdb --parm7 real.parm7 --model-pdb model.pdb -q 0 -m 1
  mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
  ```

---

(ja-ambertools--mm-parm)=
## AmberTools / mm-parm の問題

### `AmberTools preflight failed`

- **症状**：`mm-parm` が `AmberTools preflight failed. Missing required command(s): ... Required: tleap, antechamber, parmchk2` で、`all` が `[preflight] Missing required command(s) for mm_parm (AmberTools): ...` で止まる。
- **対処**：`conda install -c conda-forge ambertools=24.8 "numpy>=2,<2.5" -y` で AmberTools を入れるか、HPC では `module load ambertools` で読み込むか、ソースからビルドしてください（<https://ambermd.org/AmberTools.php>）。`which tleap antechamber parmchk2` で確かめます。AmberTools が無くても、別に作ったトポロジーを `--parm7` で渡せば個別のコマンドは動きます。

### リガンドで `antechamber` が失敗する

- **症状**：`mm-parm` が `[<RES>] antechamber failed (see log).` で止まる、または antechamber の前に `[<RES>] electron-count check failed before antechamber: ...` で止まる。
- **原因**：電荷や多重度がリガンドの水素の数と合っていないか、元素記号・結合・TER レコードが正しくありません。
- **対処**：
  - リガンドの元素記号、水素、結合、TER レコードを確かめてください。
  - 形式電荷を `-l 'LIG:-1'` で、一重項でないリガンドの多重度を `--ligand-mult 'HEM:1,NO:2'`（`all` では `--auto-mm-ligand-mult`）で与えてください。
  - `--keep-temp`（`all` では `--auto-mm-keep-temp`）を付けて再実行すると作業ディレクトリ `parm7build_*` が残るので、その中の `<resname>.antechamber.log` を読んでください。
  - リガンドに antechamber を手で実行して切り分けてください：`antechamber -i ligand.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc -nc -3 -at gaff2`。
  - RESP 電荷やほかの独自パラメータを使うときは、tleap と自分の `frcmod` / `lib` ファイルで[トポロジーを作り](mm-parm.md#使用上の注意点)、`--parm7` で渡してください。

### `Coordinate shape mismatch for '...': got (N, 3), expected (M, 3)`

- **症状**：計算がこのメッセージで止まる。
- **原因**：構造の原子が `parm7` のトポロジーと一致していません。
- **対処**：計算に使う構造から `mm-parm` でトポロジーを作り直すか、`-o` か `--add-h` で `mm-parm` が書き出す PDB で計算してください。`mm-parm` の後に PDB の原子を並べ替えないでください。

### `oniom-export` の `Element sequence mismatch at atom index ...`

- **対処**：`parm7` を作ったときと同じ PDB を `-i` に渡してください。`--no-element-check` でこの確認を外せます（結果は手で確かめます）。[原子の数が違うとき](oniom-export.md#使用上の注意点)は、確認を外しても止まります。

---

(ja-hessian_ff-build--import)=
## `hessian_ff` ビルドの問題

- **症状**：`hessian_ff build attempts failed: ...` が出る、または計算が `native bonded extension is unavailable.`、`native nonbonded extension is unavailable; torch fallback is disabled.`、`analytical Hessian native extension is unavailable.` で止まる。
- **原因**：C++ 拡張は初回の使用時に `torch.utils.cpp_extension` でコンパイルされます。C++20 対応のコンパイラ（GCC 13.3 で検証済み）と `ninja` が必要です。`ninja` は `mlmm-toolkit` と一緒に入ります。
- **対処**：
  - `g++ -std=c++20 -x c++ -fsyntax-only /dev/null` で C++20 に対応しているか確かめてください。コンパイラは `conda install -c conda-forge cxx-compiler` で入れるか、HPC では計算機の C++20 対応のコンパイラのモジュールを読み込みます。
  - PyTorch のヘッダが見つかるかを `python -c "import torch; print(torch.utils.cmake_prefix_path)"` で確かめてください。
  - ビルドはデフォルトでローカルの一時ディレクトリを使います。ネットワーク上のディレクトリ（NFS/Lustre）では PyTorch のビルドロックで止まることがあるためです。別のローカルのパスを使うときは `TORCH_EXTENSIONS_DIR` を指定してください。
  - `mlmm` を実行する Python 環境から `hessian_ff` を import できるかを確かめてください。
  - 手動でクリーンビルドするときは次を使います。

  ```bash
  cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make clean && make
  ```

---

## B-factor による層の割り当ての問題

層は B-factor に入っています：ML = 0.0、可動 MM = 10.0、固定 MM = 20.0（許容差 ±1.0）。`--detect-layer`（デフォルトで有効）がこれを読みます。

### 層の割り当てが想定と異なる・ML 領域が小さすぎる・大きすぎる

- 層を付けた PDB を分子ビューアで開き、B-factor で色分けしてください。
- `--model-pdb` が狙った原子を選んでいるか確かめてください。
- 可動 MM と固定 MM の境界は `define-layer --movable-cutoff`（デフォルト 8.0 Å）で調整します。
- Hessian に入れる原子は別に、`--hessian-cutoff` か YAML の `calc.hess_cutoff` / `calc.hess_mm_atoms` で決めます。

### B-factor が層として読まれない

- **症状**：`all` が `[all] ... does not contain a valid ML/MM B-factor partition (both ML and MM atoms are required). ...` や `[all] Automatic layer detection requires a valid 0/10/20 B-factor partition with both ML and MM atoms when extraction is skipped and --model-pdb is absent.` で止まる。
- **原因**：B-factor が層として読まれるのは、ML（0）の原子と MM（10 か 20）の原子がそれぞれ 1 つ以上あり、原子の 80% 以上がこのどれかの値を持つときだけです。
- **対処**：`define-layer` をもう一度実行し、書き出された PDB を使ってください。B-factor を任意の値に手で書き換えないでください。

### `--detect-layer` が想定どおりに働かない

- **症状**：自動で読んだ層の分け方が想定と異なる、または `-c` を付けない `all` が `[all] Skipping extraction (no -c/--center) with B-factor layer detection disabled requires --model-pdb. ...` で止まる。
- **対処**：
  - 入力は PDB か、`--ref-pdb` 付きの XYZ にしてください。
  - `define-layer` で層を付け、書き出された PDB を使ってください。
  - 計算のコマンドに `--movable-cutoff` を与えると `--detect-layer` は無効になり、MM の層は B-factor ではなく距離で決まります。

---

(ja-installation-environment-problems)=
## インストール / 環境の問題

まず、使っている環境にオプションのパッケージが入っているか、PyTorch から GPU が見えるかを確かめてください。直した後は `mlmm --version` と `python -c "import torch; print(torch.cuda.is_available())"` で確かめ、本番の前に一度 `--dry-run` を付けて実行し、オプションと入力を確かめます。

| 症状 | 原因 | 対処 |
| --- | --- | --- |
| UMA のダウンロードに失敗する（`huggingface_hub.errors.GatedRepoError`、`401`、`403`） | Hugging Face にログインしていない、または UMA モデルのライセンスに同意していない | 環境・マシンごとに一度 `hf auth login` を実行し、Hugging Face の UMA モデルのページでライセンスに同意してください。HPC では、計算ノードから Hugging Face のキャッシュディレクトリに書き込めるかを確かめます |
| `torch.cuda.is_available()` が `False`、またはインポート時に CUDA の実行時エラー | PyTorch のビルドが計算ノードの GPU・ドライバに合っていない | `nvidia-smi`、`python -m torch.utils.collect_env`、`python -m pip check` で、割り当てられた GPU、入っている wheel、ドライバを確かめてください。`nvidia-smi` が示す `CUDA Version` はドライバが扱える最も新しい CUDA です。それ以下の CUDA の PyTorch wheel（`cu126`、`cu130`、`cu132`）を入れてください |
| `--mep-mode dmf` が `DMF mode (--mep-mode dmf) requires ase, cyipopt, and pydmf>=1.2` で止まる | `cyipopt` と `pydmf` は `mlmm-toolkit` と一緒には入らない（`ase` は入る） | {ref}`インストールの手順 3 <ja-step-by-step-installation>` のとおり、`conda install -c conda-forge cyipopt -y` と `pip install 'pydmf[torch]>=1.2'`（`--dmf-backend cpu` だけなら `pip install 'pydmf>=1.2'`）を実行してください |
| 図の出力に失敗する（Plotly / Chrome） | ヘッドレス Chrome が無い | `plotly_get_chrome -y` を一度実行してください。Chromium のバイナリをダウンロードするので、インターネット接続が必要です |

### DMF が IPOPT 内で極端に遅い

IPOPT/MUMPS が並列版 BLIS を使う環境では、入れ子の並列化で長い待ちが生じることがあります。ジョブスクリプトなどで、Python や CLI の起動前に `BLIS_NUM_THREADS=1` を設定してください。外側の OpenMP/MM のスレッド数は変更不要です。`BLIS_JC_NT`、`BLIS_PC_NT`、`BLIS_IC_NT`、`BLIS_JR_NT`、`BLIS_IR_NT` の手動設定はこの制限より優先されるため、そのジョブの設定から外してください。起動済みの Notebook は、設定変更後にカーネルを再起動します。詳しくは [BLIS のスレッド設定](https://github.com/flame/blis/blob/2.0/docs/Multithreading.md)を参照してください。

---

(ja-calculation--convergence)=
## 計算 / 収束の問題

まず TS 候補を確かめてください。TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます。ν < −5.00 cm⁻¹ のモードを虚振動として数えます（YAML の `freq.zero_cutoff_cm` で変更）。その変位と IRC の端点を確かめてください。

(ja-cuda-oom)=
### CUDA メモリ不足（OOM）

次を順に試してください。

1. **固定 MM を確かめる**：`define-layer` で遠い原子が B = 20.0 になっているか確かめてください。固定 MM が小さすぎると、可動 MM とその Hessian が大きくなります。[`--movable-cutoff`](model-setup.md#可動-mm-の殻を薄くする) を小さくすると固定 MM が広がります。
2. **ML 領域を小さくする**：`extract` の `--radius` を小さくするか、`--model-pdb` で[小さい ML 領域](model-setup.md#ml-領域を小さくする)を渡します。
3. **Hessian の範囲を絞る**：`opt`、`tsopt`、`freq`、`sp` の [`--hessian-cutoff`](model-setup.md#hessian-の範囲を絞る) で、Hessian に入る可動 MM の原子を減らします。
4. **Hessian の計算方式を比べる**：有限差分は ML の自動微分のメモリを抑えることが多いものの、どちらの方式も動く原子の密な Hessian を作ります。対象の系で実行時間とピークメモリを比べ、`Analytical` を選んでいたらデフォルトの `FiniteDifference` に戻してください。
5. **メモリの大きい GPU に移る**：本番の前に、同じモデル・Hessian の方式・動く範囲のまま、移る先の GPU で短く試してください。
6. **ML を CPU で動かす**：YAML で `ml_device: cpu` を指定すると、時間はかかりますが GPU のメモリの制限を避けられます。

(ja-optimizer-stalls-with-flat-energy--forces-just-above-threshold-mlip-force-noise-floor)=
### エネルギーが平坦なまま最適化が終わらない（MLIP の力のノイズフロア）

- **症状**：`opt` / `tsopt` が回り続けるが、直近のサイクルでエネルギーがほぼ平坦になり、最大・RMS の力が `gau` / `baker` の閾値のわずか上で下がらなくなる。
- **原因**：MLIP の力には数値精度によるノイズフロアがあります。大きな ML/MM 系では、これが標準の力の閾値を上回り、構造がほぼ止まっていても力が閾値を下回りません。
- **対処**：
  - 実行の上限は `--max-cycles`（デフォルト 100000）です。エネルギーが平坦になった時点で早く止めたいときは、`--stop-plateau`（`opt`、`tsopt`、`all`）を付けてください。`stalled`（未収束）として止まります。
  - 判定の調整は [YAML 設定の一覧](yaml-reference.md#opt) を参照してください。
  - この判定は GSM / DMF では行いません。`path-opt` / `path-search` の単一構造の事前最適化では、YAML の `opt.energy_plateau: true` で有効にします。

(ja-troubleshooting-ts)=
### TS 最適化が収束しない・虚振動が複数残る

- **症状**：TS 最適化が多くのサイクルを回しても収束しない（`summary.log` に `TS optimization did not converge. Review the TS trajectory.`）、または収束後に [n_imag](tsopt.md#ts-の判定) が 2 以上（`TS imaginary-mode validation found n_imag=N.`）や 0（`[tsopt] No imaginary mode detected. Try all --refine-path.`）になる。
- **最適化が収束しないときの対処**：止まった理由とモードの変位を確かめてから、次を順に試してください。
  1. オプティマイザを RS-P-RFO（デフォルト）と Dimer 法の間で切り替える：単独では `tsopt --opt-mode hess` / `dimer`、`all` では `--opt-mode-post hess` / `grad`（Dimer）。
  2. YAML でステップサイズを小さくする。[YAML 設定の一覧](yaml-reference.md#ts-最適化セクション) を参照してください。
  3. 経路のよりよい HEI（最高エネルギーのイメージ）など、別の候補から始める。{ref}`TS が取れないとき <ja-ts-search-fails>` を参照してください。
- **n_imag ≥ 2 が残るときの対処**：`--flatten` を付けて最適化し直すか、収束の基準をデフォルトの `baker` から `gau_tight` か `gau_vtight` に締めてください（単独では `tsopt --thresh`、`all` では `--thresh-post`）。Hessian の範囲を広げる、`--refine-path` などのほかの手は {ref}`TS が取れないとき <ja-ts-search-fails>` にまとめてあります。

(ja-troubleshooting-irc)=
### IRC が正常に終了しない

IRC が収束せずに止まっても、端点の最適化で狙った R と P に着けば使えます。まず[最適化後の端点](irc.md#irc-の成否の判定)を確かめてください。

- **症状**：IRC が明確な極小構造に着く前に止まる、またはエネルギーが振動し勾配が大きいままになる。
- **原因**：この曲面にはステップが大きすぎるか、開始構造に虚振動が 2 つ以上あります。
- **対処**：
  - 単独の `irc`：`--step-size 0.05`（デフォルト 0.10 bohr）。[irc の例 4](irc.md#4-小さいステップでの再試行) を参照してください。
  - `all`：`--irc-step-size 0.05`。
  - 開始構造が n_imag = 1 であることを確かめてください。
  - 物理的な停止条件を無視してサイクルの上限まで追うには、単独で `irc --never-stop`、`all` で `--irc-never-stop` を指定し、軌跡と端点を確かめてください。

(ja-troubleshooting-mep)=
### MEP 探索（GSM / DMF）が失敗する・結合変化を取りこぼす

- **症状**：最小エネルギー経路（MEP）の探索が使える経路を作らずに終わる（`MEP optimization did not converge. Review the MEP trajectory and convergence log.`）、または予想した結合変化が出ない。
- **対処**：
  - 複雑な反応では `--max-nodes`（デフォルト 20）を 30 や 40 に増やしてください。
  - 端点の事前最適化は有効のままにしてください（デフォルト）。`--no-preopt` を付けていたら外します。
  - 別の手法を試してください：`--mep-mode dmf` ↔ `gsm`。
  - YAML の `bond.bond_factor` と `bond.delta_fraction` で結合変化の検出を調整してください。

---

(ja-troubleshooting-performance)=
## パフォーマンス / 安定性のヒント

- **メモリ不足**：{ref}`CUDA メモリ不足 <ja-cuda-oom>` の手を試すか、`--max-nodes` を減らします。`opt` と `scan` では {ref}`--opt-mode grad <ja-opt-mode-semantics>`（L-BFGS、Hessian なし）のままにします。
- **解析的な ML Hessian**：評価の回数を減らせることがありますが、メモリはバックエンドと系で変わります。試しの計算で `FiniteDifference` と比べてください。
- **MM Hessian**：デフォルトの `mm_fd: true`（有限差分）は速さよりメモリを優先します。`mm_fd: false` は小さな系では速いものの、メモリを多く使います。
- **複数の GPU**：ML は 1 つのデバイス（`ml_cuda_idx: 0`）を使います。デフォルトの `hessian_ff` の MM バックエンドは CPU で動きます。MM を別の GPU に置くときは、`mm_backend: openmm`、`mm_device: cuda`、`mm_cuda_idx: 1` を指定します。
- **ML と MM の並列実行**：デフォルトで ML（GPU）と MM（CPU）は並列に動きます。CPU のスレッド数は `mm_threads` で指定します。

(ja-troubleshooting-backends)=
## バックエンド固有の問題

### バックエンドのパッケージが無い

- **症状**：`orb-models is required for the ORB backend. ...`、`aimnet is required for the AIMNet2 backend. ...`、`mace-torch is required for the MACE backend. ...` が出る。
- **対処**：
  - ORB：`pip install "mlmm-toolkit[orb]"`。AIMNet2：`pip install "mlmm-toolkit[aimnet]"`。
  - MACE：専用の環境で `pip uninstall -y fairchem-core && pip install mace-torch` を実行してください。`mace-torch` は `e3nn==0.4.4` に固定し、UMA（`fairchem-core`）は `e3nn>=0.5` を必要とします。
  - 追加パッケージを入れても ORB を import できないときは、`python -m pip check` を実行し、依存関係の解決か import のエラーが名指ししたパッケージを直してください。関係のない PyG のパッケージは入れないでください。

---

## 不具合報告のときに添えると助かる情報

実行したコマンド、`summary.log`（または端末の出力）、再現できる最小の入力、環境（OS / Python / CUDA / PyTorch）、AmberTools と `hessian_ff` が入って動くかどうかを添えてください。

## 関連ドキュメント

- [反応機構を調べるコツ](mechanism-tips.md) — TS が取れないときに試すこと
- [インストール](installation.md) — 環境の構築とオプションのバックエンド
- [MLIP バックエンド](backends.md) — バックエンドの選び方
- [ML 領域と層の組み方](model-setup.md) — ML 領域と層を確かめる・削る・広げる
