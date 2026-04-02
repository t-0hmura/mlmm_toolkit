# `oniom-export`

## 概要

> **要約:** Amber トポロジーを持つ ML/MM 系を、外部 QM/MM 入力形式（Gaussian ONIOM または ORCA QM/MM）へエクスポートします。

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

---

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [oniom_gaussian](oniom-gaussian.md) -- Gaussian モード詳細（`--mode g16`）
- [oniom_orca](oniom-orca.md) -- ORCA モード詳細（`--mode orca`）
- [oniom_import](oniom-import.md) -- ONIOM 入力から XYZ/層付き PDB を再構築
- [mm_parm](mm-parm.md) -- Amber トポロジー構築
- [define_layer](define-layer.md) -- レイヤー定義/確認
