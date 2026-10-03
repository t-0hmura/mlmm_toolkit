"""Regression tests: every bond cut by ``extract`` ends in a kept atom with a link H."""

from __future__ import annotations

import sys

import numpy as np
import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm requires Python >= 3.11",
)

# NMET-ALA-ASP-PRO-LYS-GLY-PRO-SER-CYS-GLY-HYP-GLU-ALA-CLYS (ff14SB names), built
# with tleap and minimized in vacuum so that distance-inferred bonds match the topology.
# Columns: resseq, resname, atom, x, y, z.
_PEPTIDE = """\
1 MET N 4.092 4.217 -0.308
1 MET H1 3.649 3.308 -0.248
1 MET H2 4.372 4.382 -1.262
1 MET H3 3.451 4.937 -0.005
1 MET CA 5.278 4.214 0.552
1 MET HA 5.791 5.155 0.430
1 MET CB 4.917 4.084 2.034
1 MET HB2 3.967 3.563 2.155
1 MET HB3 5.688 3.506 2.542
1 MET CG 4.859 5.473 2.690
1 MET HG2 5.436 5.425 3.613
1 MET HG3 5.356 6.197 2.050
1 MET SD 3.228 6.142 3.086
1 MET CE 3.805 7.670 3.892
1 MET HE1 4.303 8.308 3.159
1 MET HE2 2.959 8.206 4.325
1 MET HE3 4.516 7.425 4.685
1 MET C 6.234 3.133 0.091
1 MET O 5.882 1.960 0.150
2 ALA N 7.369 3.554 -0.465
2 ALA H 7.507 4.567 -0.472
2 ALA CA 8.559 2.773 -0.821
2 ALA HA 8.794 2.083 -0.009
2 ALA CB 8.282 1.971 -2.101
2 ALA HB1 8.011 2.650 -2.911
2 ALA HB2 9.180 1.427 -2.388
2 ALA HB3 7.478 1.257 -1.935
2 ALA C 9.732 3.769 -0.973
2 ALA O 10.259 3.981 -2.061
3 ASP N 10.002 4.518 0.100
3 ASP H 9.556 4.250 0.968
3 ASP CA 10.391 5.936 0.021
3 ASP HA 10.148 6.285 -0.983
3 ASP CB 9.495 6.724 0.999
3 ASP HB2 9.895 6.617 2.008
3 ASP HB3 9.478 7.780 0.741
3 ASP CG 8.064 6.198 0.978
3 ASP OD1 7.350 6.420 -0.029
3 ASP OD2 7.720 5.382 1.858
3 ASP C 11.900 6.203 0.253
3 ASP O 12.626 5.304 0.686
4 PRO N 12.399 7.430 -0.025
4 PRO CD 11.670 8.566 -0.573
4 PRO HD2 11.182 9.104 0.242
4 PRO HD3 10.938 8.264 -1.323
4 PRO CG 12.722 9.467 -1.212
4 PRO HG2 12.403 10.509 -1.230
4 PRO HG3 12.950 9.110 -2.217
4 PRO CB 13.928 9.255 -0.303
4 PRO HB2 13.821 9.873 0.591
4 PRO HB3 14.863 9.482 -0.815
4 PRO CA 13.822 7.767 0.058
4 PRO HA 14.342 7.191 -0.708
4 PRO C 14.458 7.472 1.419
4 PRO O 13.962 7.883 2.470
5 LYS N 15.613 6.795 1.393
5 LYS H 15.967 6.505 0.492
5 LYS CA 16.357 6.351 2.575
5 LYS HA 16.187 7.065 3.383
5 LYS CB 15.803 4.978 3.001
5 LYS HB2 14.750 5.096 3.268
5 LYS HB3 15.864 4.288 2.156
5 LYS CG 16.556 4.373 4.193
5 LYS HG2 17.583 4.164 3.892
5 LYS HG3 16.562 5.087 5.018
5 LYS CD 15.928 3.059 4.660
5 LYS HD2 14.914 3.250 5.017
5 LYS HD3 15.893 2.362 3.820
5 LYS CE 16.794 2.483 5.787
5 LYS HE2 17.819 2.378 5.414
5 LYS HE3 16.817 3.205 6.610
5 LYS NZ 16.281 1.174 6.260
5 LYS HZ1 16.273 0.504 5.501
5 LYS HZ2 16.870 0.813 7.000
5 LYS HZ3 15.338 1.271 6.615
5 LYS C 17.860 6.305 2.284
5 LYS O 18.266 5.955 1.181
6 GLY N 18.669 6.609 3.297
6 GLY H 18.246 6.859 4.175
6 GLY CA 20.123 6.423 3.289
6 GLY HA2 20.441 5.945 2.361
6 GLY HA3 20.599 7.402 3.346
6 GLY C 20.618 5.565 4.465
6 GLY O 19.812 4.899 5.127
7 PRO N 21.936 5.574 4.740
7 PRO CD 22.987 6.113 3.889
7 PRO HD2 23.099 7.181 4.081
7 PRO HD3 22.787 5.936 2.832
7 PRO CG 24.253 5.373 4.311
7 PRO HG2 25.152 5.955 4.102
7 PRO HG3 24.296 4.402 3.816
7 PRO CB 24.036 5.183 5.811
7 PRO HB2 24.348 6.091 6.332
7 PRO HB3 24.587 4.322 6.193
7 PRO CA 22.519 4.992 5.950
7 PRO HA 22.289 3.926 5.971
7 PRO C 22.004 5.650 7.243
7 PRO O 21.470 6.756 7.223
8 SER N 22.210 4.980 8.377
8 SER H 22.696 4.094 8.339
8 SER CA 21.978 5.513 9.726
8 SER HA 22.298 6.556 9.755
8 SER CB 20.485 5.469 10.070
8 SER HB2 19.921 6.003 9.302
8 SER HB3 20.145 4.432 10.108
8 SER OG 20.260 6.090 11.320
8 SER HG 19.316 6.203 11.463
8 SER C 22.819 4.731 10.746
8 SER O 23.324 3.656 10.422
9 CYS N 23.002 5.270 11.952
9 CYS H 22.463 6.095 12.190
9 CYS CA 23.817 4.694 13.025
9 CYS HA 23.718 3.607 13.008
9 CYS CB 25.286 5.061 12.754
9 CYS HB2 25.536 4.807 11.722
9 CYS HB3 25.423 6.134 12.900
9 CYS SG 26.396 4.144 13.859
9 CYS HG 27.547 4.691 13.453
9 CYS C 23.326 5.201 14.398
9 CYS O 22.642 6.223 14.469
10 GLY N 23.674 4.509 15.482
10 GLY H 24.294 3.714 15.383
10 GLY CA 23.340 4.900 16.853
10 GLY HA2 23.587 5.952 16.983
10 GLY HA3 22.272 4.767 17.022
10 GLY C 24.112 4.075 17.895
10 GLY O 24.821 3.144 17.512
11 HYP N 24.008 4.409 19.193
11 HYP CD 23.314 5.562 19.750
11 HYP HD22 22.295 5.277 20.021
11 HYP HD23 23.299 6.410 19.064
11 HYP CG 24.098 5.930 21.009
11 HYP HG 23.459 6.410 21.751
11 HYP OD1 25.185 6.776 20.682
11 HYP HD1 24.923 7.696 20.784
11 HYP CB 24.597 4.574 21.512
11 HYP HB2 23.843 4.147 22.179
11 HYP HB3 25.549 4.664 22.037
11 HYP CA 24.736 3.707 20.252
11 HYP HA 25.789 3.649 19.973
11 HYP C 24.205 2.288 20.508
11 HYP O 23.047 1.990 20.228
12 GLU N 25.020 1.454 21.153
12 GLU H 25.935 1.784 21.421
12 GLU CA 24.621 0.171 21.750
12 GLU HA 23.530 0.107 21.766
12 GLU CB 25.158 -0.998 20.892
12 GLU HB2 24.865 -0.836 19.853
12 GLU HB3 26.247 -1.017 20.950
12 GLU CG 24.593 -2.354 21.341
12 GLU HG2 24.897 -2.550 22.367
12 GLU HG3 23.503 -2.302 21.324
12 GLU CD 25.065 -3.546 20.505
12 GLU OE1 24.227 -4.420 20.177
12 GLU OE2 26.290 -3.760 20.348
12 GLU C 25.114 0.123 23.211
12 GLU O 26.056 0.831 23.567
13 ALA N 24.491 -0.701 24.056
13 ALA H 23.761 -1.301 23.700
13 ALA CA 24.946 -1.003 25.413
13 ALA HA 25.987 -0.690 25.519
13 ALA CB 24.104 -0.213 26.420
13 ALA HB1 23.058 -0.509 26.334
13 ALA HB2 24.456 -0.422 27.430
13 ALA HB3 24.203 0.853 26.221
13 ALA C 24.894 -2.520 25.662
13 ALA O 24.054 -3.210 25.079
14 LYS N 25.809 -3.017 26.500
14 LYS H 26.403 -2.376 27.016
14 LYS CA 26.100 -4.427 26.796
14 LYS HA 25.180 -4.986 26.945
14 LYS CB 26.933 -5.033 25.647
14 LYS HB2 27.454 -4.235 25.114
14 LYS HB3 27.699 -5.693 26.061
14 LYS CG 26.076 -5.863 24.682
14 LYS HG2 26.146 -6.912 24.973
14 LYS HG3 25.030 -5.565 24.750
14 LYS CD 26.531 -5.692 23.230
14 LYS HD2 26.263 -4.684 22.905
14 LYS HD3 27.613 -5.819 23.154
14 LYS CE 25.830 -6.730 22.352
14 LYS HE2 26.295 -7.702 22.530
14 LYS HE3 24.780 -6.790 22.654
14 LYS NZ 25.906 -6.365 20.921
14 LYS HZ1 26.570 -5.605 20.757
14 LYS HZ2 26.047 -7.135 20.293
14 LYS HZ3 25.059 -5.860 20.641
14 LYS C 26.859 -4.525 28.120
14 LYS O 27.430 -3.484 28.517
14 LYS OXT 26.840 -5.636 28.687
"""

_RADII = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05}
_CAP_BOND = 1.09


def _atom_line(record, serial, name, resname, resseq, xyz, element):
    x, y, z = xyz
    return (
        f"{record:<6}{serial:>5} {name:<4} {resname:>3} A{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{0.00:>6.2f}"
        f"          {element:>2}\n"
    )


def _write_peptide(path):
    lines = []
    for serial, row in enumerate(_PEPTIDE.splitlines(), start=1):
        resseq, resname, name, x, y, z = row.split()
        xyz = (float(x), float(y), float(z))
        lines.append(_atom_line("ATOM", serial, name, resname, int(resseq), xyz, name[0]))
    lines.append("TER\n")
    lines.append(_atom_line("HETATM", len(lines), "C1", "LIG", 200, (40.0, 40.0, 40.0), "C"))
    path.write_text("".join(lines) + "END\n", encoding="utf-8")


def _read_model(path):
    """Return ({(resseq, atom): (element, xyz)} for the peptide, [link-H xyz])."""
    atoms, caps = {}, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        resname = line[17:20].strip()
        if resname == "LKH":
            caps.append(xyz)
        elif resname != "LIG":
            atoms[(int(line[22:26]), line[12:16].strip())] = (line[76:78].strip(), xyz)
    return atoms, caps


def _boundary_defects(source, output, add_linkh=True):
    """List cut bonds that are not closed by a link H on a kept carbon."""
    ref, _ = _read_model(source)
    kept, caps = _read_model(output)
    defects = [f"{key} moved" for key, (_, xyz) in kept.items() if np.linalg.norm(xyz - ref[key][1]) > 1e-3]
    unused = list(range(len(caps)))
    keys = list(ref)
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            (elem_a, xyz_a), (elem_b, xyz_b) = ref[a], ref[b]
            if np.linalg.norm(xyz_a - xyz_b) > 1.15 * (_RADII[elem_a] + _RADII[elem_b]):
                continue
            if (a in kept) == (b in kept):
                continue
            inner, outer = (a, b) if a in kept else (b, a)
            if ref[inner][0] != "C" or ref[outer][0] == "H":
                defects.append(f"{inner}-{outer} cut without a link H")
                continue
            if not add_linkh:
                continue
            vec = ref[outer][1] - ref[inner][1]
            target = ref[inner][1] + _CAP_BOND * vec / np.linalg.norm(vec)
            match = [k for k in unused if np.linalg.norm(caps[k] - target) < 0.01]
            if match:
                unused.remove(match[0])
            else:
                defects.append(f"{inner}-{outer} has no link H")
    defects += [f"link H #{k + 1} closes no cut bond" for k in unused]
    return defects


CASES = {
    "default_contacts": dict(center="A:LYS:5,A:CYS:9", radius=2.6),
    "default_isolated_center": dict(center="A:LYS:5", radius=0.0),
    "exclude_isolated_center": dict(center="A:LYS:5", radius=2.6, exclude_backbone=True),
    "exclude_adjacent_centers": dict(center="A:SER:8,A:CYS:9", radius=2.6, exclude_backbone=True),
    "exclude_pro_neighbors": dict(
        center="A:LIG", selected_resn="A:PRO:4,A:PRO:7", radius=0.0, exclude_backbone=True
    ),
    "default_hyp_neighbor": dict(center="A:LIG", selected_resn="A:HYP:11", radius=0.0),
    "exclude_hyp_neighbor": dict(
        center="A:LIG", selected_resn="A:HYP:11", radius=0.0, exclude_backbone=True
    ),
    "default_terminal_centers": dict(center="A:MET:1,A:LYS:14", radius=0.0),
    "exclude_terminal_centers": dict(center="A:MET:1,A:LYS:14", radius=0.0, exclude_backbone=True),
}


def _extract(tmp_path, n_inputs=1, **options):
    from mlmm.workflows.extract import extract_api

    options.setdefault("add_linkh", True)

    sources = [tmp_path / f"state_{index}.pdb" for index in range(n_inputs)]
    outputs = [tmp_path / f"model_{index}.pdb" for index in range(n_inputs)]
    for source in sources:
        _write_peptide(source)
    result = extract_api(
        [str(path) for path in sources], output=[str(path) for path in outputs], **options
    )
    return sources[0], outputs, result


@pytest.mark.parametrize("n_inputs", [1, 2])
@pytest.mark.parametrize("case", sorted(CASES))
def test_model_boundary_is_kept_or_capped(tmp_path, case, n_inputs):
    source, outputs, _ = _extract(tmp_path, n_inputs, **CASES[case])
    for output in outputs:
        assert _boundary_defects(source, output) == []


@pytest.mark.parametrize("case", sorted(CASES))
def test_model_reports_no_uncapped_boundary(tmp_path, capsys, case):
    _, _, result = _extract(tmp_path, **CASES[case])
    assert result["uncapped_boundaries"] == []
    assert "outside the carbon cut positions" not in capsys.readouterr().err


def test_pro_n_side_neighbors_keep_alpha_hydrogens(tmp_path):
    _, outputs, _ = _extract(tmp_path, **CASES["exclude_pro_neighbors"])
    atoms, _ = _read_model(outputs[0])
    names = {resseq: {name for (seq, name) in atoms if seq == resseq} for resseq in (3, 6)}
    assert {"CA", "HA", "C", "O"} <= names[3] and "N" not in names[3]
    assert {"CA", "HA2", "HA3", "C", "O"} <= names[6] and "N" not in names[6]


@pytest.mark.parametrize("case", ["default_isolated_center", "exclude_isolated_center"])
def test_isolated_center_keeps_side_chain_only(tmp_path, capsys, case):
    _, outputs, _ = _extract(tmp_path, **CASES[case])
    atoms, _ = _read_model(outputs[0])
    lys = {name for (seq, name) in atoms if seq == 5}
    assert {"CB", "NZ", "HZ1"} <= lys
    assert not lys & {"N", "H", "CA", "HA", "C", "O"}
    assert "A:LYS:5" in capsys.readouterr().out


def test_exclude_backbone_keeps_main_chain_between_adjacent_centers(tmp_path):
    _, outputs, _ = _extract(tmp_path, **CASES["exclude_adjacent_centers"])
    atoms, _ = _read_model(outputs[0])
    ser = {name for (seq, name) in atoms if seq == 8}
    cys = {name for (seq, name) in atoms if seq == 9}
    assert {"CA", "HA", "C", "O", "OG"} <= ser and "N" not in ser
    assert {"N", "H", "CA", "HA", "SG"} <= cys and "C" not in cys


@pytest.mark.parametrize("exclude_backbone", [False, True])
@pytest.mark.parametrize(
    ("center", "resseq", "terminal_atoms", "protein_charge"),
    [
        # NH3+ of Met1 counts +1.
        ("A:MET:1", 1, ("N", "H1", "H2", "H3", "CA"), 1),
        # Lys14 has its side chain (+1) and COO- (-1).
        ("A:LYS:14", 14, ("C", "O", "OXT", "CA"), 0),
    ],
)
def test_terminal_centers_keep_terminal_groups_and_charge(
    tmp_path, center, resseq, terminal_atoms, protein_charge, exclude_backbone
):
    _, outputs, result = _extract(
        tmp_path, center=center, radius=0.0, exclude_backbone=exclude_backbone
    )
    atoms, _ = _read_model(outputs[0])
    assert {(resseq, name) for name in terminal_atoms} <= set(atoms)
    assert result["charge_summary"]["protein_charge"] == protein_charge


def test_glycine_center_without_neighbors_is_reported(tmp_path, capsys):
    _, outputs, _ = _extract(tmp_path, center="A:GLY:6", radius=0.0, exclude_backbone=True)
    atoms, _ = _read_model(outputs[0])
    assert not [key for key in atoms if key[0] == 6]
    assert "A:GLY:6" in capsys.readouterr().err


def test_no_add_linkh_leaves_only_carbon_cuts(tmp_path, capsys):
    options = dict(CASES["exclude_adjacent_centers"], add_linkh=False)
    source, outputs, result = _extract(tmp_path, **options)
    assert _boundary_defects(source, outputs[0], add_linkh=False) == []
    assert _read_model(outputs[0])[1] == []
    assert result["uncapped_boundaries"] == []
    assert "outside the carbon cut positions" not in capsys.readouterr().err
