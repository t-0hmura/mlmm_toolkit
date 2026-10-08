# インストール

`mlmm-toolkit` は Linux 環境（ワークステーションや HPC）向けで、本番計算では通常 CUDA 対応 GPU を使用します。MM の計算には、トポロジーを作る **AmberTools**（`tleap`）と、`hessian_ff` のカーネルをビルドする **C++20 対応コンパイラー** が必要です。どちらも下の conda のコマンドで入ります。

## クイックスタート

`nvidia-smi` の右上に出る `CUDA Version` は、ドライバーが扱える最も新しい CUDA です。PyTorch の wheel は、それ以下の CUDA のもの（`cu126`、`cu130`、`cu132`）を選びます。以下のコマンドは推奨の `cu130` を使います。

### 必須

```bash
# 1) AmberTools・PDBFixer・C++ コンパイラーを入れた conda 環境を作成
# 2) CUDA 対応の PyTorch ビルドをインストール
# 3) mlmm-toolkit をインストール
# 4) Plotly 静的画像 (PNG) エクスポート用のヘッドレス Chrome をインストール
#    Chromium のバイナリをダウンロード（インターネット接続が必要）

conda create -n mlmm-toolkit python=3.12 -y
conda activate mlmm-toolkit
conda install -c conda-forge ambertools=24.8 "numpy>=2,<2.5" pdbfixer cxx-compiler -y
TORCH_INDEX=cu130  # 推奨。cu126 / cu132 も選べる
pip install 'torch==2.13.0' --index-url "https://download.pytorch.org/whl/${TORCH_INDEX}"
pip install mlmm-toolkit
plotly_get_chrome -y
```

最後に、UMA モデルをダウンロードできるように **Hugging Face Hub** にログインします。無料の HF アカウントと読み取り専用トークンが必要で、先に <https://huggingface.co/facebook/UMA> で FAIR Chemistry License v1 に同意します:

```bash
hf auth login
# またはスクリプト内でトークン指定する場合:
hf auth login --token '<YOUR_ACCESS_TOKEN>' --add-to-git-credential
```

これはマシン/環境ごとに 1 回だけ行う必要があります。最後に `mlmm --version` でインストールを確かめます。

### 任意

DMF を使う場合は、cyipopt と pydmf も導入してください（[下記の手順 3](#詳細なインストール手順)）。ORB・AIMNet2・MACE・DFT の導入は[手順 7](#詳細なインストール手順) を参照してください。

(ja-step-by-step-installation)=
## 詳細なインストール手順

環境を段階的に構築する場合:

1. **クラスターやビルドが必要とする場合のみ CUDA toolkit をロード**

    ビルド済みの PyTorch wheel に `nvcc` は不要です。依存パッケージをソースからビルドする場合は `module avail cuda` を確認し、クラスターが指定するコンパイラと CUDA toolkit の組み合わせをロードしてください:

    ```bash
    module load cuda/<your-version>   # 例: cuda/12.6 または cuda/12.9
    ```

2. **AmberTools を入れた conda 環境を作成**

    システムの AmberTools のモジュールがあるクラスターでは、conda の AmberTools と衝突しないように、先に `module unload amber` を実行してください。

    ```bash
    conda create -n <your-env> python=3.12 -y
    conda activate <your-env>
    conda install -c conda-forge ambertools=24.8 "numpy>=2,<2.5" pdbfixer cxx-compiler -y
    ```

3. **cyipopt と pydmf をインストール**
    MEP 探索で DMF 法（`--mep-mode dmf`）を使用する場合に必要です。どちらも `mlmm-toolkit` と一緒には入りません。GSM のみを使用する場合はスキップできます。それでも `--mep-mode dmf` がインポートのエラーで止まるときは、{ref}`インストール / 環境の問題 <ja-installation--environment>` を参照してください。

    ```bash
    conda install -c conda-forge cyipopt -y
    pip install 'pydmf[torch]>=1.2'   # --dmf-backend cpu だけなら pip install 'pydmf>=1.2'
    ```

4. **適切な CUDA ビルドの PyTorch をインストール**

    推奨の例（`cu130`）:

    ```bash
    pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130
    ```

    PyTorch 2.13.0 の公式 wheel には `cu126`、`cu132`、`cpu` もあります。上のクイックスタートにある `nvidia-smi` の見方で wheel を選び、手順 8 で GPU が使えるかを確かめてください。[PyTorch の版の対応表](https://pytorch.org/get-started/previous-versions/) を参照してください。

5. **`mlmm-toolkit` 本体と可視化用 Chrome をインストール**

    ```bash
    pip install mlmm-toolkit
    plotly_get_chrome -y
    ```

    `hessian_ff` のカーネルは初回の使用時に自動でビルドされます。ビルドが失敗したときの手動の再ビルドは、{ref}`hessian_ff ビルドの問題 <ja-hessian_ff-build--import>` を参照してください。

6. **Hugging Face Hub (UMA モデル) にログイン**

    ```bash
    hf auth login
    ```

    利用許諾と非対話型ログインは、上の「必須」を参照してください。

    詳細は上流プロジェクトを参照してください:

    - fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
    - Hugging Face トークンとセキュリティ: <https://huggingface.co/docs/hub/security-tokens>

7. **（任意）追加の MLIP バックエンドと追加パッケージをインストール**

    mlmm-toolkit はデフォルトで UMA を使用します。他のバックエンドは対応する追加パッケージを導入し、`-b/--backend`（例: `-b orb`）で選択します:

    **ORB**（Python 3.11／3.12 が必要、3.12 推奨）:

    ```bash
    pip install "mlmm-toolkit[orb]"
    ```

    **AIMNet2**:

    ```bash
    pip install "mlmm-toolkit[aimnet]"
    ```

    **MACE**: `mace-torch` が要求する `e3nn==0.4.4` が UMA の `fairchem-core` と衝突するので、手順 2〜6 で作った別の環境に入れます。

    ```bash
    pip uninstall -y fairchem-core
    pip install mace-torch
    ```

    **DFT**（`-b dft`、`--dft`、`mlmm dft`）: `[dft]` は Linux x86_64 で CUDA 13 用の GPU4PySCF を導入します（手順 4 の cu130 / cu132 の PyTorch 向け）。cu126 の wheel を使うときは、代わりに `[dft-cuda12]` を導入します。aarch64 では [GPU4PySCF](https://github.com/pyscf/gpu4pyscf) をソースからビルドしてください。

    ```bash
    pip install "mlmm-toolkit[dft]"
    ```

    DFT/MM を使う場面と、MLIP/MM の TS を DFT/MM で確かめる方法は [MLIP の TS を DFT で確かめる](dft-backend.md) を参照してください。ほかに、`[openmm]`（MM に OpenMM）、`[mcp]`（`mlmm-mcp` サーバー）、`[pdbfixer]`（pip で PDBFixer）の追加パッケージがあります。

8. **インストールの確認**

    ```bash
    mlmm --version
    mlmm -h
    hf auth whoami
    ```

    1 行目でインストールされたバージョンが、2 行目でサブコマンドの一覧が、3 行目で Hugging Face のユーザー名が表示されます。GPU アクセスの確認:

    ```bash
    python -c "import torch; print('CUDA:', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"
    ```

    `CUDA: False` の場合、バージョンを変える前に、インストールした wheel、スケジューラがジョブに GPU を見せているか、ドライバー、環境のライブラリを確認してください:

    ```bash
    python -m torch.utils.collect_env
    python -m pip check
    ```

## システム要件

**OS:** Linux を推奨します。ネイティブの Windows では AmberTools（`tleap`）が使えないため、対応していません。

**Python:** 3.12 を推奨します（最低 3.11）。ORB バックエンドには 3.11 か 3.12 が必要です。

**GPU / CUDA:** 選んだ wheel に対応するドライバーの NVIDIA GPU（クイックスタートを参照）。新しい GPU アーキテクチャでは新しい wheel が必要なことがあります。CPU のみでも実行できますが、通常は大幅に遅くなります。

**AmberTools とコンパイラー:** `mm-parm` と `all` のトポロジー作成に AmberTools（`tleap`）が必要です。前の計算の対応する `--parm7` を渡すと、この段を省けます。既定の MM バックエンド `hessian_ff` には C++20 対応コンパイラー（GCC 13.3 で検証）が必要です。PDBFixer が要るのは `mm-parm --add-h` のときだけです。

**VRAM・RAM・ディスク:** メモリはバックエンド、原子数、Hessian の計算方式とともに増え、ディスクには環境、モデルの重み、トポロジー、軌跡と Hessian が入ります。代表的な計算を 1 つ対象の計算ノードで流し、最大使用量を見てください。

## 次のステップ

- [はじめに](getting-started.md) — 最短の実行と、次に読むページ
- [クイックスタート: `mlmm all`](quickstart-all.md) — R と P から MEP を作る
- [クイックスタート: `mlmm all --scan-lists`](quickstart-scan.md) — 1 つの構造から経路を作る
- [クイックスタート: TS-only モード](quickstart-tsopt.md) — TS 候補を最適化して確かめる
- [MLIP の TS を DFT で確かめる](dft-backend.md) — TS を DFT/MM で詰めて確かめる
- [共通オプションと残基・原子の指定](cli-conventions.md) — 共通のオプションと、残基・原子の指定の書き方
- [デバイス設定 & HPC セットアップ](device-hpc.md) — クラスターでの GPU の設定とジョブスクリプト
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと試すこと
