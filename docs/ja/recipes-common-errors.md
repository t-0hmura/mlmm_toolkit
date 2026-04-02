# 典型エラー別レシピ

症状は分かるが、どのコマンドページから見ればよいか迷うときの入口です。
詳細は [トラブルシューティング](troubleshooting.md) を並行して参照してください。

## クイックルーティング

| 症状 | まず確認 | 次に読む |
| --- | --- | --- |
| 元素カラム欠損 / 抽出が中断 | 元の PDB に対して `add-elem-info` を実行 | [入力 / 抽出](troubleshooting.md#入力--抽出の問題) |
| "Charge is required" エラー | `-q/--charge` と `-m/--multiplicity` を明示的に指定 | [電荷 / スピン](troubleshooting.md#電荷--スピンの問題) |
| 実行後にエネルギー/状態がおかしい | CLI 規約の電荷/多重度ポリシーを再確認 | [電荷 / スピン](troubleshooting.md#電荷--スピンの問題) |
| `mm-parm` が実行できない（`tleap`/`antechamber`/`parmchk2` が見つからない） | AmberTools の利用可能性を先に修正 | [AmberTools / mm-parm](troubleshooting.md#ambertools--mm-parm-の問題) |
| `hessian_ff` の import/ビルドエラー | ネイティブ拡張（`hessian_ff/native`）を再ビルド | [hessian_ff ビルド](troubleshooting.md#hessian_ff-ビルドの問題) |
| DMF モードの import エラー（`cyipopt`） | アクティブ環境に `cyipopt` をインストール | [DMF モード](troubleshooting.md#dmf-モードが動かないcyipopt-がない) |
| TSOPT/IRC が収束しない | ステップ長を縮小（RFO/RS-I-RFO では trust_radius、L-BFGS では max_step）、サイクル数を増やし、まず TS の品質を検証 | [計算 / 収束](troubleshooting.md#計算--収束の問題) |
| CUDA/GPU ランタイム不整合 | `torch.cuda.is_available()` と CUDA ビルドの組み合わせを確認 | [CUDA / PyTorch](troubleshooting.md#cuda--pytorch-の不整合) |
| プロット出力の失敗 | Plotly エクスポート用の Chrome ランタイムをインストール | [プロット出力](troubleshooting.md#図のエクスポートが失敗するchrome-がない) |

## レシピ 1: MEP 前に抽出で止まる

兆候:
 - 元素情報不足、原子数不一致、ポケット抽出結果が空に近い。
最初の確認:
 - 入力構造が同じ前処理フローで作られ、原子順が揃っているか。
 - `extract` / `all` 前に元素カラムが埋まっているか。
典型的な修正手順:
 - 元素修復 -> 抽出再実行 -> ポケットサイズ/残基選択を再確認。

## レシピ 2: 電荷/スピンの解決で止まる

兆候:
 - ログに不自然な電荷が表示される場合（例: タンパク質電荷が間違っている、総電荷が期待と異なる）、電荷解決ルールを確認する。
最初の確認:
 - 対象状態に対して総電荷・多重度が妥当か。
 - `--ligand-charge` の残基キーが入力構造と一致しているか。
 - 結果が物理的に不自然な場合は [CLI Conventions](cli-conventions.md) の解決ルールを再確認。
典型的な修正手順:
 - 重要な実行では `-q` / `-m` を明示し、scan/path/tsopt を再試行。

## レシピ 3: ビルド・環境依存で止まる

兆候:
 - `mm-parm` 必須ツール不在、`hessian_ff` import 失敗、CUDA 不整合。
最初の確認:
 - 実行環境に必要実行ファイル・拡張モジュールが存在するか。
 - GPU 可視性と PyTorch CUDA 互換性に問題がないか。
典型的な修正手順:
 - 先にツールチェイン/ビルドを修復し、`--help-advanced` で再確認後に本実行。

## レシピ 4: 収束・後処理で止まる

兆候:
 - TSOPT が停滞、IRC が不安定、MEP 精密化が途中停止。
最初の確認:
 - TS 候補が支配的な虚振動数モード 1 本を持つか。
 - ステップ長（trust_radius / max_step）を縮小し、サイクル上限を増やす。
典型的な修正手順:
 - 小規模ケースで条件を詰め、安定化後に本番条件へ戻す。
