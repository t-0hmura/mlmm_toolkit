# mlmm MCP サーバー

AI エージェントから MCP（Model Context Protocol）で mlmm-toolkit の 22 個のツールを呼ぶための、インストール、ツールの一覧、クライアントの設定をまとめたページです。サーバー `mlmm-mcp` は stdio 上の JSON-RPC でやりとりするので、[MCP](https://modelcontextprotocol.io/) に対応したどのクライアントからも使えます。Claude Desktop、Claude Code、Cursor、Codeium のほか、公式の Python や TypeScript の MCP SDK で作ったエージェントからも使えます。

## インストール

```bash
pip install "mlmm-toolkit[mcp]"
```

これにより `mcp[cli]` 依存関係が追加され、`mlmm-mcp` コンソールスクリプトが登録されます。

## ツール

22 個のツールがあり、それぞれが CLI サブコマンドに 1 対 1 で対応します。各ツールは次のフィールドを持つ構造化された dict を返します。

- `schema_version`: 結果の形式の版（`"2.0"`）。各レスポンスでこの値を読むと、どのフィールドがあるかが分かります。
- `execution_status`: `completed` | `failed`
- `scientific_status`: `success` | `partial` | `failed`
- `summary_status`: 中身のある `summary` が返るのは `ok` のときだけです
  - `ok`: この呼び出しの `summary` を読めた
  - `not_required`: `out_dir` を持たず要約を書かないツール
  - `summary_missing`: `out_dir` に `summary.json` が無い
  - `summary_parse_error`: `summary.json`、または `detect_bond_changes` で出力された JSON を JSON の object として読めない
  - `summary_run_mismatch`: ファイルが別の実行のもの、または隣の `result.json` と一致しない
- `exit_code`: CLI のプロセスの終了コード
- `out_dir`: ステージランナーとスキャン / 経路 / パイプラインのツールの出力ディレクトリ。ほかのツールでは null
- `summary`: 読み込んだ `summary.json`。`detect_bond_changes` では `mlmm bond-summary --json` が出力する JSON。`out_dir` を持たないほかのツールでは空の object
- `stderr_tail` / `stdout_tail`: プロセス出力の末尾約 60 行
- `hint`: CLI のエラーメッセージの末尾の `; recover: <hint>` にある対処のヒント（ある場合）
- `argv`: 実行したコマンドライン全体（再現性のため）
- `run_id`: この呼び出しの UUID

各ツールの必須の引数は下の表にあります。そのうち入力のパス（`input_pdb`・`reactant_pdb` など）はそのコマンドの `-i` に渡ります。`parm7` は系全体のトポロジー（`--parm7`）です。省略できる引数はそのコマンドの CLI オプションを指定します。例えば `charge` は `-q`、`ligand_charge` は `-l`、`max_cycles` は `--max-cycles` です。どのツールも、追加の CLI フラグを並べた文字列のリスト `extra_args` と `timeout_seconds` を受け付け、ステージランナーとスキャン / 経路 / パイプラインのツールは `out_dir` も受け付けます。すべての引数とその型は、クライアントがツールの一覧と一緒に受け取る入力スキーマにあります。

### 構造化されたエラーエンベロープ

ステージランナーとスキャン / 経路のツールが失敗すると、返された `summary` に次のエラーのフィールドが入ります。エージェントはテキストをパースせずに、エラーのクラスで場合分けできます。

- `error`: エラーメッセージ
- `error_type`: 例外クラス名
- `error_class_chain`: そのクラスと親クラスの名前を、具体的なものから順に並べたもの（例: `["OptimizationError", "RuntimeError", "Exception", "BaseException"]`）
- `error_module`: 例外クラスが定義されているモジュール
- `error_label`: 上位レベルの CLI ステージラベル

### トポロジー / 層の準備

| MCP ツール | 必須の引数 | CLI サブコマンド | 目的 |
|---|---|---|---|
| `prepare_amber_topology` | `input_pdb`, `output_prefix` | `mlmm mm-parm` | AmberTools で系全体の AMBER parm7/rst7 を作る |
| `define_layer` | `input_pdb`, `output_pdb`、および `model_pdb` か `model_indices` | `mlmm define-layer` | ML / Movable-MM / Frozen の層を B-factor に書く |
| `extract_pocket` | `complex_pdb`, `ligand_id`, `radius_angstrom`, `output_pdb` | `mlmm extract` | 活性部位モデル: `ligand_id`（`-c`）で指定した中心から `radius_angstrom` 以内の残基を切り出す |

### ステージランナー

| MCP ツール | 必須の引数 | CLI サブコマンド | 目的 |
|---|---|---|---|
| `optimize_geometry` | `input_pdb`, `parm7` | `mlmm opt` | ONIOM 構造最適化（デフォルトは L-BFGS。RFO ではマイクロイテレーションを使う） |
| `find_transition_state` | `input_pdb`, `parm7` | `mlmm tsopt` | ONIOM TS 探索（RS-P-RFO / Dimer / RS-I-RFO / TRIM） |
| `run_irc` | `input_pdb`, `parm7` | `mlmm irc` | TS 構造からの ONIOM IRC 積分 |
| `compute_frequencies` | `input_pdb`, `parm7` | `mlmm freq` | ONIOM 振動解析 + 熱化学 |
| `run_single_point_oniom` | `input_pdb`, `parm7` | `mlmm sp` | ONIOM 一点エネルギー + 原子間力（+ `do_hess` で Hessian） |

### スキャン / 経路 / パイプライン

| MCP ツール | 必須の引数 | CLI サブコマンド | 目的 |
|---|---|---|---|
| `scan_1d` / `scan_2d` / `scan_3d` | `input_pdb`, `parm7`, `scan_lists` | `mlmm scan` / `mlmm scan2d` / `mlmm scan3d` | 調和拘束による ONIOM スキャン |
| `optimize_path` | `reactant_pdb`, `product_pdb`, `parm7` | `mlmm path-opt` | 2 端点間の ONIOM MEP 最適化 |
| `search_paths` | `input_pdb`, `product_pdb`, `parm7` | `mlmm path-search` | 再帰的な ONIOM 反応経路探索 |
| `run_full_pipeline` | `reactant_complex_pdb` | `mlmm all` | エンドツーエンド: extract → MEP → TS → IRC → freq → DFT |
| `run_single_point_dft` | `input_pdb`, `parm7` | `mlmm dft` | ML 領域の一点 DFT を MM のエネルギーと組み合わせる（GPU4PySCF または PySCF） |

### ONIOM 入出力（Gaussian / ORCA）

| MCP ツール | 必須の引数 | CLI サブコマンド | 目的 |
|---|---|---|---|
| `export_oniom_input` | `input_layered_pdb`, `parm7`, `charge`, `multiplicity`, `output_file` | `mlmm oniom-export` | Gaussian g16 または ORCA の ONIOM 入力を書き出す（`format_engine`: デフォルトは `g16`、ほかに `orca`） |
| `import_oniom_input` | `input_file`, `output_prefix` | `mlmm oniom-import` | Gaussian / ORCA の ONIOM 入力を XYZ と層付き PDB に読み戻す |

### 構造 / I/O ヘルパー

| MCP ツール | 必須の引数 | CLI サブコマンド | 目的 |
|---|---|---|---|
| `add_element_info` | `input_pdb`, `output_pdb` | `mlmm add-elem-info` | PDB の元素列を修復 |
| `fix_altloc` | `input_pdb`, `output_pdb` | `mlmm fix-altloc` | PDB の代替位置（altloc）を解決 |
| `plot_trajectory` | `input_trj_xyz`, `output_png` | `mlmm trj2fig` | エネルギープロファイル図（PNG。JPEG/SVG/PDF/HTML/CSV も可） |
| `plot_energy_diagram` | `energies`, `output_png` | `mlmm energy-diagram` | 与えた値からの状態エネルギー図 |
| `detect_bond_changes` | `reactant_pdb`, `product_pdb` | `mlmm bond-summary` | 2 つの構造（XYZ / PDB / GJF）の間の結合変化 |

### 電荷と順序付き入力

`charge` は系全体ではなく ML 領域の電荷（`-q`）です。残基名から求めたいときは、`charge` を省き、`"SAM:1,GPP:-3"` のような残基名ごとの `ligand_charge` を渡します。和を取る範囲は ML 領域で、入力の B-factor の層か `--model-pdb` から決まります。`--model-pdb` は `run_single_point_oniom` では `model_pdb`、ほかのツールでは `extra_args` で渡します。詳しくは {ref}`ML/MM の共通オプション <ja-mlmm-options>` を参照してください。

`search_paths` では、反応物の `input_pdb` と `product_pdb` の間に入る中間体を、順番に並べて `intermediate_pdbs` に渡します。`scan_1d`・`scan_2d`・`scan_3d` は `scan_lists` を `--scan-lists` の 1 つの値として渡します。書き方は [`scan`](scan.md) のページにあります。`run_full_pipeline` に `reactant_complex_pdb` だけを渡すときは、1 つの入力の `mlmm all` と同じく、`do_tsopt=True` か、`extra_args` で渡すスキャン（`["--scan-lists", "…"]`）が必要です。

## IRC と TS 最適化の設定

IRC と TS の引数は同じ名前の CLI オプションです。意味とデフォルトは各コマンドのページにあります。

- `run_irc`（`step_size`・`irc_pos_def`）: [`irc`](irc.md)。`--irc-pos-def` は [自動生成のオプションの一覧（英語のみ）](../reference/commands/irc.md) にあります
- `find_transition_state`（`opt_mode`・`microiter`・`flatten`）: [`tsopt`](tsopt.md) の `--opt-mode`。デフォルトは `hess`（RS-P-RFO）です。`microiter=False` でマイクロイテレーションを切ります。{ref}`コマンドごとの --opt-mode <ja-opt-mode-semantics>` も参照してください
- `run_full_pipeline`（`refine_path`・`do_tsopt`・`do_thermo`・`do_dft`・`thresh_post`）: [`all`](all.md) の `--refine-path`・`--tsopt`・`--thermo`・`--dft`。`--thresh-post` は [自動生成のオプションの一覧（英語のみ）](../reference/commands/all.md) にあります

## クライアント設定

クライアントごとに設定スキーマは異なります。次のスニペットはトップレベルの `mcpServers` オブジェクトを受け付けるクライアント用です。設定ファイルとスキーマは各クライアントの MCP ドキュメントで確認してください。

- Claude Desktop — `~/Library/Application Support/Claude/claude_desktop_config.json`（macOS） / `%APPDATA%\Claude\claude_desktop_config.json`（Windows）
- Cursor — `~/.cursor/mcp.json`
- Claude Code — ファイルは編集せず、`claude mcp add mlmm -- mlmm-mcp` を実行します

クライアントがサーバーを起動すると、ツールの一覧に 22 個のツールが出ます。

```json
{
  "mcpServers": {
    "mlmm": {
      "command": "mlmm-mcp",
      "args": []
    }
  }
}
```

環境変数（PATH / AMBERHOME / CUDA_VISIBLE_DEVICES）を設定する完全な例は [`examples/mcp_client_config.json`](../../examples/mcp_client_config.json) を参照してください。

VS Code は `.vscode/mcp.json` で[トップレベルの `servers` オブジェクト](https://code.visualstudio.com/docs/agents/reference/mcp-configuration)を使います。

```json
{
  "servers": {
    "mlmm": {
      "command": "mlmm-mcp",
      "args": []
    }
  }
}
```

### カスタム Python MCP クライアント

```python
import asyncio

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

async def main():
    server_params = StdioServerParameters(command="mlmm-mcp")
    async with stdio_client(server_params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            result = await session.call_tool(
                "optimize_geometry",
                arguments={
                    "input_pdb": "r_complex_layered.pdb",
                    "parm7": "real.parm7",
                    "charge": 0,
                    "max_cycles": 50,
                },
            )
            print(result.content)

asyncio.run(main())
```

## サンドボックス / 安全性に関する注意

- 各ツールは、サーバーと同じ Python インタープリタで `mlmm`（`python -m mlmm`）をサブプロセスとして実行し、作業ディレクトリは呼び出した側のままです。このため入力の相対パスの意味は変わらず、`PATH` の前にある別の `mlmm` も使われません。
- サーバーは呼び出した側の PATH、conda 環境、CUDA の設定、AmberTools のパスを引き継ぎます。`prepare_amber_topology` には PATH 上の AmberTools（antechamber、parmchk2、tleap）が必要です。opt / tsopt / irc / scan のような時間のかかるツールは CLI をサブプロセスで実行するので、呼び出しごとに `timeout_seconds` を設定し、止まらない計算を打ち切ってください（既定は時間制限なし）。
- ステージランナーとスキャン / 経路 / パイプラインのツールの出力は `out_dir` の下に置かれます。指定しないときは呼び出しごとに別の一時ディレクトリ（例: `mlmm_mcp_opt_…`）を使うので、同時の呼び出しがぶつかりません。
- ほかのツールは `out_dir` を持たず、指定した出力パスに書きます。`extra_args` で CLI のフラグを追加できますが、型付きの出力パス、`--out-dir`、`--out-json/--no-out-json`、`detect_bond_changes` の `--json/--no-json` は、`--out-dir=…` のようにつなげた形も含めて上書きできず、そのような呼び出しは CLI を起動する前に拒否されます。コマンドに渡したパスは、返された `argv` ですべて確かめられます。
- サーバーはファイルシステムのサンドボックスではなく、外部プログラムやモデルのキャッシュは `out_dir` の外にも書き込むことがあります。入力の構造、parm7、MLIP の重みは前もってディスクに置き、ファイルへのアクセスを限りたいときは、クライアント側のパス制限か OS / コンテナによる隔離を使ってください。

## 使用上の注意点

- `run_full_pipeline`、`out_dir` を持たないツール、`summary.json` を書く前に止まった実行（`summary_missing`）では、エラーのフィールドの代わりに `stderr_tail` と `hint` を読んでください。
- `--model-indices` を使うとき、または `--model-pdb` 無しで `--no-detect-layer` を使うときは `ligand_charge` から電荷を求められないので、`charge` を渡してください。
- `charge` と `ligand_charge` の両方を渡すと、明示した `charge` が優先されます。

## 関連ドキュメント

* [JSON 出力の一覧](json-output.md) — ツールが返す状態の欄と `summary.json`
* [トラブルシューティング](troubleshooting.md) — 実行に失敗したときの対処
* [コマンドの一覧（英語のみ）](../reference/commands/index.md) — 各ツールの元の CLI オプション（`extra_args` 用）
