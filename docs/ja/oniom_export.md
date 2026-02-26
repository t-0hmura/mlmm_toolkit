# ONIOM エクスポート（`oniom-gaussian` / `oniom-orca`）

## 概要

> **要約:** Amber トポロジを持つ ML/MM 系を、外部 QM/MM 入力形式（Gaussian ONIOM または ORCA QM/MM）へエクスポートします。

### 概要
- **用途:** mlmm_toolkit で準備した系に対して外部 QM/MM 計算を実行したい場合。
- **手法:** Amber parm7 トポロジを読み込み、モデル PDB から ML 領域をマッピングし、Gaussian `.com`/`.gjf` または ORCA `.inp` を出力。
- **出力:** 適切なレイヤー/結合情報を含む Gaussian ONIOM 入力または ORCA QM/MM 入力。
- **次のステップ:** Gaussian または ORCA でエクスポートした入力を実行。

## コマンド別ガイド

- [oniom-gaussian](oniom_gaussian.md) -- Gaussian 向けの詳細説明・オプション・例
- [oniom-orca](oniom_orca.md) -- ORCA 向けの詳細説明・オプション・例

## 使い分け

| 必要なもの | 推奨コマンド |
| --- | --- |
| リンク原子注釈付き Gaussian ONIOM 入力 | [`oniom-gaussian`](oniom_gaussian.md) |
| ORCAFF 連携込みの ORCA QM/MM 入力 | [`oniom-orca`](oniom_orca.md) |

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [oniom-gaussian](oniom_gaussian.md) -- Gaussian ONIOM エクスポータ
- [oniom-orca](oniom_orca.md) -- ORCA QM/MM エクスポータ
- [mm_parm](mm_parm.md) -- Amber トポロジ構築
- [define_layer](define_layer.md) -- レイヤー定義/確認
