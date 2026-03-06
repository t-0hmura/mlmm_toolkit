# ONIOM エクスポート（`oniom-export`）

## 概要

> **要約:** Amber トポロジーを持つ ML/MM 系を、外部 QM/MM 入力形式（Gaussian ONIOM または ORCA QM/MM）へエクスポートします。

### 早見表
- **用途:** mlmm_toolkit で準備した系に対して外部 QM/MM 計算を実行したい場合。
- **手法:** Amber parm7 トポロジーを読み込み、モデル PDB から ML 領域をマッピングし、`--mode`（未指定時は `-o` 拡張子推定）に応じて Gaussian `.com`/`.gjf` または ORCA `.inp` を出力。
- **出力:** 適切なレイヤー/結合情報を含む Gaussian ONIOM 入力または ORCA QM/MM 入力。
- **次のステップ:** Gaussian または ORCA でエクスポートした入力を実行。

## 統合コマンド

```bash
mlmm oniom-export --parm real.parm7 -i pocket.pdb --model-pdb ml.pdb \
 -o out.<gjf|com|inp> --mode <g16|orca> -q 0 -m 1
```

- `--mode` が最優先です。
- `--mode` 未指定時は `-o` から推定します。
  - `.gjf` / `.com` -> `g16`
  - `.inp` -> `orca`
- `--mode` 未指定かつ `-o` が未知拡張子の場合はエラーになります。

## 使い分け

| 必要なもの | 推奨コマンド |
| --- | --- |
| リンク原子注釈付き Gaussian ONIOM 入力 | `mlmm oniom-export --mode g16` |
| ORCAFF 連携込みの ORCA QM/MM 入力 | `mlmm oniom-export --mode orca` |

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [oniom_gaussian](oniom_gaussian.md) -- Gaussian モード詳細（`--mode g16`）
- [oniom_orca](oniom_orca.md) -- ORCA モード詳細（`--mode orca`）
- [oniom_import](oniom_import.md) -- ONIOM 入力から XYZ/層付き PDB を再構築
- [mm_parm](mm_parm.md) -- Amber トポロジー構築
- [define_layer](define_layer.md) -- レイヤー定義/確認
