# `oniom-export`

Amber トポロジーを持つ ML/MM 系を、外部 QM/MM 入力ファイル（Gaussian ONIOM = `--mode g16`、または ORCA QM/MM = `--mode orca`）へエクスポートします。`parm7` と MLMM の layer 情報を持つ PDB を読み込み、可動/固定原子は PDB の B-factor から決定します。`--model-pdb` を指定した場合は、その原子を QM 領域として使用します。

両モードとも CMAP を含まない `parm7` が必要です。Gaussian ONIOM は
CMAP を忠実に表現できず、ORCA の MM エンジンも CMAP を適用しないため、
CMAP を含むトポロジーでは出力前に停止します。これはエクスポート形式の
制約であり、通常の mlmm 計算では両 MM 層の CMAP を有効にできます。

## 実行例

```bash
# Gaussian ONIOM 入力
mlmm oniom-export --parm real.parm7 -i pocket_layered.pdb --model-pdb ml.pdb \
 -o out.gjf --mode g16 -q 0 -m 1
```

```bash
# ORCA QM/MM 入力（.inp 拡張子からモード推定）
mlmm oniom-export --parm real.parm7 -i pocket_layered.pdb --model-pdb ml.pdb \
 -o out.inp -q 0 -m 1
```

```bash
# メソッド/基底とリソースを指定した Gaussian 入力
mlmm oniom-export --parm real.parm7 -i pocket_layered.pdb --model-pdb ml.pdb \
 -o out.gjf --mode g16 --method 'wb97xd/def2-svp' --nproc 16 --mem 32GB -q 0 -m 1
```

## 処理の流れ

1. **トポロジー + layer** -- `parm7` と `-i` の layered PDB を読み込みます（原子順序はトポロジーと一致が必須）。B-factor 0/10/20 が ML、可動 MM、固定 MM を表します。
2. **QM 領域** -- `--model-pdb` を指定しなければ B-factor の ML layer を使用します。可動/固定原子は常に layered PDB から取得します。
3. **QM/MM 境界** -- Gaussian は `--link-atom-method scaled`（デフォルトの Morokuma/Dapprich g-factor）または `fixed`（1.09/1.01 Å）で link H を配置します。ORCA は `QMAtoms`/`ORCAFF` から cap を生成し、出力中の link 座標は診断コメントです。
4. **書き出し** -- `-o` に対象形式の入力ファイルを出力します。ORCA モードでは `ORCAFF.prms` のパスも特定します。`--convert-orcaff` が有効なら `orca_mm -convff -AMBER` による変換を試みます。変換が無効または利用不能でも `.inp` は書き出され、ORCA 実行前に用意すべきパラメータパスを報告します。

## 出力

- `<output>.{gjf,com}`（g16）または `<output>.inp`（ORCA） -- QM/MM 入力ファイル
- ORCA モードでは `<parm7_stem>.ORCAFF.prms` を参照します。既存ファイルを再利用し、自動変換が有効かつ利用可能な場合だけ生成します。生成された `.inp` を実行する前に、この参照先が実在することを確認してください

## CLI オプション

全フラグの一覧は生成済みの[コマンドリファレンス](../reference/commands/index.md)にあります。以下の表は説明が必要なオプションを扱います。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--parm FILE` | Amber parm7 トポロジーファイル | 必須 |
| `-i, --input FILE` | MLMM layered PDB。原子順序は parm7 と一致し、B-factor が可動/固定原子を定義 | 必須 |
| `--model-pdb FILE` | QM 領域原子を定義する PDB | _None_ |
| `-o, --output FILE` | 出力ファイルパス（g16 は `.gjf` / `.com`、ORCA は `.inp`） | 必須 |
| `--mode [g16\|orca]` | エクスポートモード。未指定時は `-o` 拡張子から推定 | _推定_ |
| `--method TEXT` | QM メソッドと基底関数 | モード依存 |
| `-q, --charge INT` | QM 領域の電荷 | 必須 |
| `-m, --multiplicity INT` | QM 領域の多重度 | `1` |
| `--nproc INT` | プロセッサ数 | `8` |
| `--mem TEXT` | メモリ割り当て（g16 モード） | `16GB` |
| `--total-charge INT` / `--total-mult INT` | 全 QM+MM 系の総電荷/総多重度（ORCA `Charge_Total` / `Mult_Total`） | トポロジー由来 / `--multiplicity` と同じ |
| `--orcaff PATH` | `ORCAFF.prms` のパス（ORCA モード）。未指定時は派生パスを参照し、条件を満たす場合に自動生成を試行 | _None_ |
| `--convert-orcaff / --no-convert-orcaff` | `ORCAFF.prms` 欠損時に `orca_mm -convff -AMBER` で自動変換（ORCA モード） | `True` |
| `--element-check / --no-element-check` | `--input` の元素配列を parm7 トポロジーと照合 | `True` |
| `--link-atom-method [scaled\|fixed]` | Gaussian の link H 配置。ORCA では対応座標を診断コメントとして記録し、cap は `QMAtoms`/`ORCAFF` から生成 | `scaled` |

`mlmm oniom-export --help` はコアオプション、`mlmm oniom-export --help-advanced` は全オプションを表示します。

## 注記

- モード選択: `--mode` が最優先です。`--mode` 未指定時は `-o` から推定します。
  - `.gjf` / `.com` -> `g16`
  - `.inp` -> `orca`
- `--mode` 未指定かつ `-o` が未知拡張子の場合はエラーになります。
- PDB/ENT 入力では `MLMM_REF_PDB_ORDER_V1_SHA256=<digest>` を埋め込みます。
  座標・occupancy・B-factor は除外し、固定フィールドの原子名・残基・chain・
  insertion code・元素 identity を対象にします。`oniom-import --ref-pdb` は
  positional metadata を復元する前にこの marker を検証します。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [oniom-gaussian](oniom-gaussian.md) -- Gaussian モード詳細（`--mode g16`）
- [oniom-orca](oniom-orca.md) -- ORCA モード詳細（`--mode orca`）
- [oniom-import](oniom-import.md) -- ONIOM 入力から XYZ/層付き PDB を再構築
- [mm-parm](mm-parm.md) -- Amber トポロジー構築
- [define-layer](define-layer.md) -- レイヤー定義/確認
