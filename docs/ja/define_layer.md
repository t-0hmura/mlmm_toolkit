# `define-layer`

## 概要

> **要約:** ML 領域からの距離に基づいて 3 層 ML/MM システムを定義し、レイヤー割り当てを出力 PDB の B 因子としてエンコードします。

`mlmm define-layer` は、酵素系を ML 領域の周囲に 3 つのレイヤーに分割し、割り当てを PDB の B 因子として書き出します。ML 領域はモデル PDB、明示的な原子インデックス、またはその組み合わせで指定できます。

### 3 層システム

| レイヤー | 名前 | B 因子 | 説明 |
| --- | --- | --- | --- |
| 1 | ML | 0.0 | ML 領域の原子 |
| 2 | Movable-MM | 10.0 | `--radius-freeze` 以内の MM 原子/残基 |
| 3 | Frozen | 20.0 | `--radius-freeze` より遠い MM 原子/残基 |

### レイヤー割り当て戦略

- **ML 原子を含まない残基:** 残基全体が、任意の ML 原子から残基内の任意の原子までの最小距離に基づいて単一のレイヤーに割り当てられます。
- **ML 原子を含む残基:** 同じ残基内の非 ML 原子は距離に基づいて個別に分類されます。

## 使用法

```bash
mlmm define-layer -i INPUT.pdb --model-pdb MODEL.pdb [--model-indices "0,1,2,..."] \
    [--radius-partial-hessian FLOAT] [--radius-freeze FLOAT] \
    [-o OUTPUT.pdb] [--one-based|--zero-based]
```

### 例

```bash
# モデル PDB を使用して ML 領域を定義
mlmm define-layer -i system.pdb --model-pdb ml_region.pdb -o labeled.pdb

# 明示的な原子インデックスを使用（0 始まり）
mlmm define-layer -i system.pdb --model-indices "0,1,2,3,4" --zero-based -o labeled.pdb

# カスタム半径
mlmm define-layer -i system.pdb --model-pdb ml_region.pdb \
    --radius-freeze 10.0 -o labeled.pdb
```

## 説明

1. **ML 領域の同定** -- ML 領域は `--model-pdb`（入力 PDB との原子マッチング）または `--model-indices`（明示的な原子インデックス）で定義されます。`--model-indices` が指定された場合、`--model-pdb` より優先されます。
2. **距離計算** -- 各非 ML 原子（または残基）について、任意の ML 原子からの最小距離を計算します。
3. **レイヤー割り当て** -- 非 ML 原子/残基は `--radius-freeze` に基づいて Movable-MM または Frozen に割り当てられます。
4. **出力** -- 出力 PDB の B 因子がレイヤー値（0、10、20）に設定されます。レイヤー割り当てのサマリーがコンソールに出力されます。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 全系を含む入力 PDB ファイル。 | 必須 |
| `--model-pdb PATH` | ML 領域の原子を定義する PDB ファイル。 | _None_ |
| `--model-indices TEXT` | ML 領域のカンマ区切り原子インデックス（例: `"0,1,2,3"` または `"1-10,15,20-25"`）。`--model-pdb` より優先。 | _None_ |
| `--radius-partial-hessian FLOAT` | 3 層モードでは非推奨（後方互換のため残存）。 | `0.0` |
| `--radius-freeze FLOAT` | ML 領域からの Movable-MM の距離カットオフ (A)。これを超える原子は Frozen。 | `8.0` |
| `-o, --output PATH` | B 因子がレイヤー値に設定された出力 PDB ファイル。 | `<input>_layered.pdb` |
| `--one-based / --zero-based` | `--model-indices` を 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |

## 出力

```
<output>.pdb                # B 因子が 0 / 10 / 20 に設定された PDB
(stdout)                    # レイヤー割り当てと原子数のサマリーテーブル
```

## 注意事項

- 距離は任意の ML 原子から残基内の任意の原子までの最小距離として計算されます。
- `--radius-partial-hessian` は後方互換のため受け付けますが、3 層モードでは無視されます。
- 出力 PDB は元のすべての原子レコードを保持します。変更されるのは B 因子カラムのみです。

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [mm_parm](mm_parm.md) -- レイヤー定義前に AMBER トポロジー（parm7/rst7）を構築
- [opt](opt.md) -- レイヤー化されたシステムを使用した単一構造最適化
- [all](all.md) -- 自動レイヤー定義を含むエンドツーエンドワークフロー
