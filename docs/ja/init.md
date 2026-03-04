# `init`（削除済み）

> **注意:** `init` サブコマンドは削除されました。各サブコマンドで `--show-config` を使用して現在の設定を確認するか、[YAML リファレンス](yaml_reference.md)を参照して手動で YAML ファイルを作成してください。ML バックエンドの選択（`--backend`）や xTB 点電荷埋め込み補正（`--embedcharge`）などの新しいオプションについては、各サブコマンドのドキュメントを参照してください。

## 以前の動作

`init` コマンドはスターター YAML テンプレートを生成していました：

```text
mlmm init --out mlmm_all.config.yaml
mlmm all --config mlmm_all.config.yaml --dry-run
```

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-o, --out PATH` | 生成先 YAML パス。 | `mlmm_all.config.yaml` |
| `--force/--no-force` | 既存ファイルを上書き。 | `False` |

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- 生成される YAML は最小テンプレートです（完全スキーマではありません）。

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [YAML リファレンス](yaml_reference.md) -- 設定キーの全体像
