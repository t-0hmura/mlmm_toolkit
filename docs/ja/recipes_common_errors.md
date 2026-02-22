# 典型エラー別レシピ

症状は分かるが、どのコマンドページから見ればよいか迷うときの入口です。
詳細は [トラブルシューティング](troubleshooting.md) を並行して参照してください。

## 早見表

| 症状 | 最初にやること | 詳細 |
| --- | --- | --- |
| 元素カラム欠落で抽出が止まる | 元の PDB に `add-elem-info` を適用 | [トラブルシューティング](troubleshooting.md) |
| 「電荷が必須」系エラー | `-q/--charge` と `-m/--multiplicity` を明示指定 | [トラブルシューティング](troubleshooting.md) |
| 計算は通るが状態/エネルギーが不自然 | CLI 規約の電荷/多重度ルールを再確認 | [トラブルシューティング](troubleshooting.md) |
| `mm-parm` が動かない（`tleap`/`antechamber`/`parmchk2`） | AmberTools の導入状態を先に修正 | [トラブルシューティング](troubleshooting.md) |
| `hessian_ff` の import/build エラー | `hessian_ff/native` を再ビルド | [トラブルシューティング](troubleshooting.md) |
| DMF モードの import エラー（`cyipopt`） | 実行環境に `cyipopt` を導入 | [トラブルシューティング](troubleshooting.md) |
| TSOPT/IRC が収束しない | ステップを緩め、サイクル上限を増やす | [トラブルシューティング](troubleshooting.md) |
| CUDA/GPU 実行時エラー | `torch.cuda.is_available()` と CUDA 組み合わせを確認 | [トラブルシューティング](troubleshooting.md) |
| 図の出力失敗 | Plotly 用 Chrome ランタイムを導入 | [トラブルシューティング](troubleshooting.md) |

## レシピ 1: MEP 前に抽出で止まる

- 兆候:
  - 元素情報不足、原子数不一致、ポケット抽出結果が空に近い。
- 最初の確認:
  - 入力構造が同じ前処理フローで作られ、原子順が揃っているか。
  - `extract` / `all` 前に元素カラムが埋まっているか。
- 典型的な修正手順:
  - 元素修復 -> 抽出再実行 -> ポケットサイズ/残基選択を再確認。

## レシピ 2: 電荷/スピンの解決で止まる

- 兆候:
  - 特に非 `.gjf` 入力で総電荷未解決エラーが出る。
- 最初の確認:
  - 対象状態に対して総電荷・多重度が妥当か。
  - `--ligand-charge` の残基キーが入力構造と一致しているか。
  - 結果が物理的に不自然な場合は [CLI Conventions](cli_conventions.md) の解決ルールを再確認。
- 典型的な修正手順:
  - 重要な実行では `-q` / `-m` を明示し、scan/path/tsopt を再試行。

## レシピ 3: ビルド・環境依存で止まる

- 兆候:
  - `mm-parm` 必須ツール不在、`hessian_ff` import 失敗、CUDA 不整合。
- 最初の確認:
  - 実行環境に必要実行ファイル・拡張モジュールが存在するか。
  - GPU 可視性と PyTorch CUDA 互換性に問題がないか。
- 典型的な修正手順:
  - 先にツールチェイン/ビルドを修復し、`--dry-run` で再確認後に本実行。

## レシピ 4: 収束・後処理で止まる

- 兆候:
  - TSOPT が停滞、IRC が不安定、MEP 精密化が途中停止。
- 最初の確認:
  - TS 候補が支配的な虚数モード 1 本を持つか。
  - 最適化条件（閾値、ステップ、サイクル上限）が攻め過ぎていないか。
- 典型的な修正手順:
  - 小規模ケースで条件を詰め、安定化後に本番条件へ戻す。
