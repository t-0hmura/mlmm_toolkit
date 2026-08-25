# mlmm MCP サーバー

`mlmm-mcp` は [MCP](https://modelcontextprotocol.io/) サーバーで、MCP に対応した任意のエージェントが、stdio 上の JSON-RPC を介してすべての `mlmm` CLI サブコマンドを実行できるようにします。対象には Claude Desktop / Claude Code / Cursor / Codeium のほか、公式の Python または TypeScript MCP SDK 上に構築された任意のカスタムエージェントが含まれます。

## インストール

```bash
pip install "mlmm-toolkit[mcp]"
```

これにより `mcp[cli]` 依存関係が追加され、`mlmm-mcp` コンソールスクリプトが登録されます。

## ツール

22 個のツールがあり、CLI サブコマンドごとに 1 つ対応します。各ツールは、以下を含む構造化された dict（`mlmm.mcp._runner` の `SubcmdResultDict`）を返します:

- `schema_version`: エンベロープのバージョン。実際の値: `mlmm.mcp._runner.MCP_SUBCMD_RESULT_SCHEMA_VERSION`。値の上昇はフィールドセットや値の型の変更を示します。リテラル値をハードコードせず、定数を参照してください。
- `status`: `ok` | `failed` | `summary_missing` | `summary_parse_error` | `summary_run_mismatch`
- `exit_code`: サブプロセスの終了コード
- `out_dir`: CLI が書き込んだ作業ディレクトリ
- `summary`: パースされた `summary.json`（CLI 出力スキーマ。ステージごとの形状は [JSON 出力リファレンス](json-output.md) を参照）
- `stderr_tail` / `stdout_tail`: プロセス出力の末尾 ~60 行
- `hint`: CLI エラーメッセージからパースされた `; recover: <hint>` サフィックス（存在する場合）
- `argv`: 実行された完全な argv（再現性のため）
- `run_id`: このサブプロセス呼び出しに割り当てられた UUID

型付き Python のユーザー向けに、`mlmm.mcp._runner` は `SubcmdResultDict`（実行時ペイロードを反映する `TypedDict`）と `MCP_SUBCMD_RESULT_STATUSES`（許可される `status` 文字列の列挙タプル）も公開しています。

サーバーは `mlmm` を、MCP サーバーを実行中の Python interpreter と
import 済み module（`python -m mlmm`）へ固定します。その source root を
child `PYTHONPATH` の先頭に置き、working directory は変更しません。
各呼び出しは `MLMM_RUN_ID` で `run_id` を渡し、その ID と一致する
`summary.json` だけを返します。leaf command では同一世代・同一 byte の
`result.json` も必須です。aggregate の `all` と `path-search` は一つの
summary のみを発行します。

出力先は typed tool parameter が所有します。summary tool の `extra_args`
では `-o`, `--out-dir`, `--out-json`, `--no-out-json` を、utility では typed
output option を上書きできません。短縮 option の連結形と
`--option=value` 形も subprocess 起動前に拒否されます。

### 構造化エラーエンベロープ

サブコマンドが失敗すると、パースされた `summary`（または同じディレクトリの `result.json`）に拡張エラーエンベロープが含まれ、エージェントがテキストをパースせずに例外クラス階層をパターンマッチできます:

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
| `scan_1d` / `scan_2d` / `scan_3d` | `mlmm scan` / `mlmm scan2d` / `mlmm scan3d` | ONIOM 拘束スキャン |
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
| `add_element_info` | `mlmm add-elem-info` | PDB の元素列を修復 |
| `fix_altloc` | `mlmm fix-altloc` | PDB の代替位置標識を解決 |
| `plot_trajectory` | `mlmm trj2fig` | エネルギープロファイルの PNG / HTML / SVG / PDF |
| `plot_energy_diagram` | `mlmm energy-diagram` | カテゴリ別エネルギー図 |
| `detect_bond_changes` | `mlmm bond-summary` | 2 つの PDB 間の結合変化 diff |

## オプトインの IRC 収束ガード

`run_irc` は `irc_pos_def: bool` を受け付けます。これを指定すると、IRC の収束に加えて正定値の質量重み付き Hessian も必要となり、rms のみの基準が局所極小に到達する前に成功と判定してしまう IRC の「ショルダー」誤収束をブロックします。デフォルトは `None`（rms のみ、レガシー）です。

`find_transition_state` のデフォルト `opt_mode="hess"` は RS-P-RFO を選択します。Hessian TS オプティマイザは `rsprfo`、`rsirfo`、`trim` でも明示的に選択できます。3 種はいずれも microiteration に対応し、デフォルトで有効です。無効化するには `microiter=False` を渡します。

## クライアント設定

クライアントごとに設定スキーマは異なります。次のスニペットはトップレベルの
`mcpServers` オブジェクトを受け付けるクライアント用です。設定ファイルと
スキーマは各クライアントの MCP ドキュメントで確認してください。

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

VS Code は `.vscode/mcp.json` でトップレベルの `servers` オブジェクトを使います
（[VS Code MCP 設定リファレンス](https://code.visualstudio.com/docs/agents/reference/mcp-configuration)）。

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

## サンドボックス / 安全性に関する注意

- MCP サーバーは外部 tool 用の PATH、conda 環境、CUDA セットアップ、AmberTools のパスを継承します。各ツールは server interpreter から import 済み `mlmm` module をサブプロセスとして起動するため、長時間実行されるツール（opt / tsopt / irc / scan）はプロセス外で実行されます。各呼び出しで `timeout_seconds` を設定して上限を設けてください。
- 計算コマンドのステージ出力は `out_dir` 配下に配置され、デフォルトでは一意の一時ディレクトリを使います。準備・変換 tool は明示的な出力 path を受け取り、外部 program、model cache、一時ファイル library はステージ外へ書き込む場合があります。
- MCP server は filesystem sandbox ではありません。必要な入力（PDB、parm7、ML 重み）はあらかじめディスク上に用意し、filesystem の封じ込めが必要なら client 側の path 制限または OS/container isolation を使用してください。
