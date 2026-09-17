# Gaussian ONIOM モード（`oniom-export --mode g16`）

Amber parm7 トポロジーを用いて、ML/MM システムを Gaussian ONIOM（`.com`/`.gjf`）入力へエクスポートします。`--parm` からトポロジーを、`-i/--input` の MLMM layered PDB から可動/固定原子を読み取ります。

入力 `parm7` は CMAP を含まない必要があります。Gaussian ONIOM は
CMAP 項を忠実に表現できないため、CMAP が存在する場合は出力前に停止します。

## 実行例

電荷・多重度を明示した最小構成のエクスポート:

```bash
mlmm oniom-export --mode g16 --parm real.parm7 -i pocket_layered.pdb --model-pdb ml_region.pdb \
 -o system.com -q 0 -m 1
```

```bash
# メソッドを明示して出力。
mlmm oniom-export --mode g16 --parm real.parm7 -i pocket_layered.pdb --model-pdb ml_region.pdb \
 -o system.com -q 0 -m 1 --method "wB97XD/def2-TZVPD"
```

```bash
# 元素順チェックを無効化。
mlmm oniom-export --mode g16 --parm real.parm7 -i pocket_layered.pdb --model-pdb ml_region.pdb \
 -o system.gjf -q 0 -m 1 --no-element-check
```

```bash
# 実行環境パラメータを調整。
mlmm oniom-export --mode g16 --parm real.parm7 -i pocket_layered.pdb --model-pdb ml_region.pdb \
 -o system.com -q 0 -m 1 --nproc 16 --mem 32GB
```

## 処理の流れ

1. parm7 から原子・結合・電荷情報を取得。
2. `-i/--input` の B-factor から可動/固定 layer を読み、`--element-check` で元素順を検証。
3. `--model-pdb` がある場合は QM 領域をマッピング。
4. QM/MM 境界を検出し、リンク原子情報を付与。
5. `%nprocshared`/`%mem`/method/座標/レイヤー/結合情報を含む入力を生成。

## 出力

- `system.com` または `system.gjf`
- コンソールに QM 原子数や境界処理の要約

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--parm PATH` | Amber parm7 トポロジー。 | 必須 |
| `-i, --input PATH` | MLMM layered PDB。原子順は parm7 と一致必須。 | 必須 |
| `--element-check / --no-element-check` | 入力と parm7 の元素順を検証。 | `True` |
| `--model-pdb PATH` | QM 領域原子を定義する PDB。 | _None_ |
| `-o, --output PATH` | 出力 Gaussian 入力（`.com`/`.gjf`）。 | 必須 |
| `--method TEXT` | QM メソッド/基底。 | `wB97XD/def2-TZVPD` |
| `-q, --charge INT` | QM 領域電荷。 | 必須 |
| `-m, --multiplicity INT` | QM 領域多重度。 | `1` |
| `--nproc INT` | 使用コア数。 | `8` |
| `--mem TEXT` | メモリ指定。 | `16GB` |

## 関連項目

- [oniom-orca](oniom-orca.md) -- ORCA モードガイド（`--mode orca`）
- [oniom-export](oniom-export.md) -- エクスポート全体ガイド
- [mm-parm](mm-parm.md) -- Amber トポロジー構築
- [define-layer](define-layer.md) -- レイヤー定義/確認
- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
