# mlmm MCP サーバー

`mlmm-mcp` は [MCP](https://modelcontextprotocol.io/) サーバーで、MCP に対応した任意のエージェントが、stdio 上の JSON-RPC を介してすべての `mlmm` CLI サブコマンドを実行できるようにします。対象には Claude Desktop / Claude Code / Cursor / Codeium のほか、公式の Python または TypeScript MCP SDK 上に構築された任意のカスタムエージェントが含まれます。

## インストール

```bash
pip install "mlmm[mcp]"
```

これにより `mcp[cli]` 依存関係が追加され、`mlmm-mcp` コンソールスクリプトが登録されます。

## ツール

22 個のツールがあり、CLI サブコマンドごとに 1 つ対応します。各ツールは、以下を含む構造化された dict（`mlmm.mcp._runner` の `SubcmdResultDict`）を返します:

- `schema_version`: エンベロープのバージョン。実際の値: `mlmm.mcp._runner.MCP_SUBCMD_RESULT_SCHEMA_VERSION`。値の上昇はフィールドセットや値の型の変更を示します。このドキュメント内のリテラルではなく定数に対してピン留めしてください。
- `status`: `ok` | `failed` | `summary_missing` | `summary_parse_error`
- `exit_code`: サブプロセスの終了コード
- `out_dir`: CLI が書き込んだ作業ディレクトリ
- `summary`: パースされた `summary.json`（CLI 出力スキーマ。ステージごとの形状は [JSON 出力リファレンス](json-output.md) を参照）
- `stderr_tail` / `stdout_tail`: プロセス出力の末尾 ~60 行
- `hint`: CLI エラーメッセージからパースされた `; recover: <hint>` サフィックス（存在する場合）
- `argv`: 実行された完全な argv（再現性のため）

型付き Python の利用者向けに、`mlmm.mcp._runner` は `SubcmdResultDict`（実行時ペイロードを反映する `TypedDict`）と `MCP_SUBCMD_RESULT_STATUSES`（許可される `status` 文字列の列挙タプル）も公開しています。

### 構造化エラーエンベロープ

サブコマンドが失敗すると、パースされた `summary`（または兄弟の `result.json`）に拡張エラーエンベロープが含まれ、エージェントがテキストをパースせずに例外クラス階層をパターンマッチできます:

- `error`: 元の例外の `str(exc)`
- `error_type`: 例外クラス名（例: `"OptimizationError"`）
- `error_class_chain`: MRO のクラス名（例: `["OptimizationError", "RuntimeError", "Exception", "BaseException"]`）
- `error_module`: 例外クラスを定義しているモジュール
- `error_label`: 高レベルの CLI ステージラベル（例: `"opt"`、`"tsopt-stage"`）

### トポロジー / 層の準備（mlmm 固有）

| MCP ツール | CLI サブコマンド | 目的 |
|---|---|---|
| `prepare_amber_topology` | `mlmm mm-parm` | AmberTools を介して AMBER parm7/rst7 を生成 |
| `define_layer` | `mlmm define-layer` | ML / MM-movable / MM-frozen の B-factor 層を割り当て |
| `extract_pocket` | `mlmm extract` | リガンド周辺の球を切り出して活性部位モデルを作成 |

### ステージランナー（ONIOM 対応）

| MCP ツール | CLI サブコマンド | 目的 |
|---|---|---|
| `optimize_geometry` | `mlmm opt` | ONIOM 構造最適化（microiter マクロ / マイクロ） |
| `find_transition_state` | `mlmm tsopt` | ONIOM TS 探索（RS-I-RFO / Dimer / TRIM / RS-P-RFO） |
| `run_irc` | `mlmm irc` | TS 構造からの ONIOM IRC 積分 |
| `compute_frequencies` | `mlmm freq` | ONIOM 振動解析 + 熱化学 |
| `run_single_point_oniom` | `mlmm sp` | ONIOM 一点エネルギー + 力（+ 任意の Hessian） |

### スキャン / 経路 / パイプライン

| MCP ツール | CLI サブコマンド | 目的 |
|---|---|---|
| `scan_1d` / `scan_2d` / `scan_3d` | `mlmm scan{,2d,3d}` | ONIOM 拘束スキャン |
| `optimize_path` | `mlmm path-opt` | 両端点 ONIOM MEP 最適化 |
| `search_paths` | `mlmm path-search` | 再帰的 ONIOM 経路探索 |
| `run_full_pipeline` | `mlmm all` | エンドツーエンド: extract → MEP → TS → IRC → freq → DFT |
| `run_single_point_dft` | `mlmm dft` | gpu4pyscf を介した ONIOM 埋め込み一点 DFT |

### ONIOM 入出力（Gaussian / ORCA）

| MCP ツール | CLI サブコマンド | 目的 |
|---|---|---|
| `export_oniom_input` | `mlmm oniom-export` | Gaussian g16 / ORCA の ONIOM 入力デックを書き出し |
| `import_oniom_input` | `mlmm oniom-import` | Gaussian / ORCA の ONIOM 入力を XYZ / 層付き PDB に読み戻し |

### 構造 / 入出力ヘルパー

| MCP ツール | CLI サブコマンド | 目的 |
|---|---|---|
| `add_element_info` | `mlmm add-elem-info` | PDB の元素カラムを修復 |
| `fix_altloc` | `mlmm fix-altloc` | PDB の代替位置標識を解決 |
| `plot_trajectory` | `mlmm trj2fig` | エネルギープロファイルの PNG / HTML / SVG / PDF |
| `plot_energy_diagram` | `mlmm energy-diagram` | カテゴリ別エネルギー図 |
| `detect_bond_changes` | `mlmm bond-summary` | 2 つの PDB 間の結合変化 diff |

## オプトインの IRC 収束ガード

`run_irc` は `irc_pos_def: bool` を受け付けます。これを指定すると、IRC の収束に加えて正定値の質量重み付き Hessian も必要となり、rms のみの基準が局所極小に到達する前に成功と判定してしまう IRC の「ショルダー」誤収束をブロックします。デフォルトは `None`（rms のみ、レガシー）です。

`find_transition_state` は代替の TS オプティマイザとして `opt_mode="trim"`（Helgaker 1991）/ `opt_mode="rsprfo"`（Banerjee 1985）を受け付けます。これらのモードは microiter に対応していないため、サーバーは自動的に `--no-microiter` を渡します。

## クライアント設定

すべての MCP クライアントは同じ `mcpServers` スキーマを取ります。以下のスニペットを、クライアントの MCP 設定ファイル（macOS では Claude Desktop の `~/Library/Application Support/Claude/claude_desktop_config.json`、Cursor の `~/.cursor/mcp.json` など）に貼り付けてください:

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

明示的な環境変数オーバーライド（PATH / AMBERHOME / CUDA_VISIBLE_DEVICES）を含む完全な例については、[`examples/mcp_client_config.json`](../../examples/mcp_client_config.json) を参照してください。

## サンドボックス / 安全性に関する注意

- MCP サーバーは、呼び出し元の環境の PATH、conda 環境、CUDA セットアップ、AmberTools のパスを継承します。各ツールは `mlmm` CLI をサブプロセスとして起動するため、長時間実行されるツール（opt / tsopt / irc / scan）はプロセス外で実行されます。各呼び出しで `timeout_seconds` を設定して上限を設けてください。
- 出力ファイルは `out_dir` kwarg の配下に配置されます（デフォルトは一意の `tempfile.mkdtemp("mlmm_mcp_<subcmd>_…")`）。
- サーバーは `~/.bashrc` / ログイン環境を変更したり、ソフトウェアをインストールしたり、`out_dir` の外に書き込んだりしません。必要な入力（PDB、parm7、ML 重み）はあらかじめディスク上に存在している必要があります。
