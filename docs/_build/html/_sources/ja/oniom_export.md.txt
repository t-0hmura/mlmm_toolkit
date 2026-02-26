# ONIOM エクスポート（`oniom-gaussian` / `oniom-orca`）

## 概要

Amber トポロジを持つ ML/MM 系を、外部 QM/MM 入力形式へ変換するためのサブコマンドです。

- `mlmm oniom-gaussian` -> Gaussian ONIOM（`.com` / `.gjf`）
- `mlmm oniom-orca` -> ORCA QM/MM（`.inp`）

## コマンド別ガイド

- [oniom-gaussian](oniom_gaussian.md) -- Gaussian 向けの詳細説明・オプション・例
- [oniom-orca](oniom_orca.md) -- ORCA 向けの詳細説明・オプション・例

## 使い分け

| 必要なもの | 推奨コマンド |
| --- | --- |
| リンク原子注釈付き Gaussian ONIOM 入力 | [`oniom-gaussian`](oniom_gaussian.md) |
| ORCAFF 連携込みの ORCA QM/MM 入力 | [`oniom-orca`](oniom_orca.md) |

---

## 関連項目

- [reference/commands/oniom_gaussian](../reference/commands/oniom_gaussian.md)
- [reference/commands/oniom_orca](../reference/commands/oniom_orca.md)
- [mm_parm](mm_parm.md)
- [define_layer](define_layer.md)
