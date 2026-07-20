# mmCIF と大規模構造

`mlmm-toolkit` は `.cif` / `.mmcif` と、PDB 固定幅の原子通し番号・残基番号の
上限を超える PDB を受け付けます。数値計算は内部 PDB で行います。入力
bridge が安全に再採番した一時 PDB を作成し、元の atom-site metadata を保持し、
出力時に元の chain・残基 ID を `.cif` companion へ復元します。

## 保持される情報

- 原子順、元素、原子名・残基名、座標、occupancy、B-factor、formal charge
- 複数文字の chain ID、10,000 以上の残基番号、insertion code
- 最初の座標 model（複数 model 入力では最初だけを使用）
- 残基ごとに occupancy で選んだ一貫した altloc conformer

内部 PDB の ID は実装上の一時値です。chain・残基 ID を利用する場合は、
bridge 入力に対して出力される CIF companion を使用してください。

## 対応コマンド

bridge は `all`, `extract`, `define-layer`, `sp`, `opt`, `tsopt`, `freq`,
`irc`, `dft`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search` と、
各コマンドの `--ref-pdb` で使われます。standalone `mm-parm` の入力は PDB
です。`all` は自動 parameterization の前に CIF を内部 PDB へ変換します。

```bash
mlmm all -i reactant.cif product.cif \
    -c 'enzyme_A:SAM:10001,enzyme_A:GPP:10002' \
    -l 'SAM:1,GPP:-3' --tsopt --thermo -o result

mlmm tsopt -i hei.xyz --ref-pdb full_system.mmcif \
    --parm full_system.parm7 -q -2 -o result_tsopt
```

## 一意な selector

残基名や残基番号が繰り返される構造では chain を含めて指定します。

| 用途 | 形式 | 例 |
|---|---|---|
| `extract` / `all -c` の ID 指定 | `CHAIN:RESSEQ[ICODE]` | `enzyme_A:10001B` |
| 残基名と ID の指定 | `CHAIN:RESNAME:RESSEQ[ICODE]` | `enzyme_A:SAM:10001B` |
| scan 原子指定 | `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` | `enzyme_A:SAM:10001B:CS1` |

`CHAIN:RESNAME` は同じ chain 内の該当残基をすべて選びます。production run
では一意な selector を使い、CLI が表示する選択結果を確認してください。

## Amber topology の契約

Amber `parm7` は位置対応です。全系入力と原子数・原子順が一致しなければ
なりません。model PDB は ML 領域を選ぶための、全系の原子順を保った部分集合
です。`mlmm-toolkit` は座標を割り当てる前に原子数と既知元素の並びを検証し、
不一致を検出すると停止します。同元素どうしの入れ替えはこの検証だけでは判別
できないため、R/IM/P で原子名・残基 ID・chain・原子順を保持してください。

## 上限とエラー

- 内部 bridge は最大 619,938 残基（62 internal chain × 9,999 残基）に対応し、
  超過時は ID を切り詰めず error にします。
- 座標値は内部 PDB の固定幅範囲に入る必要があります。必要なら構造全体を
  原点付近へ平行移動してください。
- decimal overflow / hybrid-36 の serial・残基番号を持つ PDB は自動正規化します。
- mmCIF の `_atom_site.type_symbol` 欠落、非有限座標、反応構造間の原子数不一致は
  hard error です。

ML 領域の境界は [概念](concepts.md) を、電荷と
scan selector は [CLI 規約](cli-conventions.md) を参照してください。
