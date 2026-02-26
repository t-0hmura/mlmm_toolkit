# `init`

## 概要

> **要約:** `mlmm all` 用のスターター YAML 設定ファイルを生成します。

### 概要
- **用途:** `mlmm all` の再現性の高い設定ファイル中心のワークフローを始めたい場合。
- **手法:** すべての設定キーのデフォルト値を含む YAML テンプレートを生成。
- **出力:** `mlmm_all.config.yaml`（または `--out` で指定したパス）。
- **デフォルト:** 出力先 `mlmm_all.config.yaml`、`--force` なしでは上書きしない。
- **次のステップ:** 生成した YAML を編集し、`mlmm all --config mlmm_all.config.yaml` を実行。

## 最小例

```bash
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
