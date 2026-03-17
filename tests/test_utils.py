"""Unit tests for mlmm.utils module."""

import sys
from pathlib import Path
import tempfile
import pytest
import numpy as np

# Skip on Python < 3.11
pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm requires Python >= 3.11",
)


def test_parse_indices_string():
    """Test parse_indices_string utility."""
    from mlmm.utils import parse_indices_string

    # Simple comma-separated list (1-based → 0-based by default)
    assert parse_indices_string("1,2,3") == [0, 1, 2]

    # Range notation
    assert parse_indices_string("1-5") == [0, 1, 2, 3, 4]

    # Mixed
    assert parse_indices_string("1,3-5,7") == [0, 2, 3, 4, 6]

    # Empty or None
    assert parse_indices_string("") == []
    assert parse_indices_string(None) == []

    # Whitespace handling
    assert parse_indices_string(" 1 , 2 , 3 ") == [0, 1, 2]

    # 0-based mode
    assert parse_indices_string("0,1,2", one_based=False) == [0, 1, 2]
    assert parse_indices_string("1,2,3", one_based=False) == [1, 2, 3]


def test_distance_A_from_coords():
    """Test distance calculation utility (input in Bohr, output in Å)."""
    from mlmm.utils import distance_A_from_coords
    from pysisyphus.constants import ANG2BOHR

    # Simple 3-atom system in Bohr units
    coords_bohr = np.array([
        [0.0, 0.0, 0.0],
        [1.0 * ANG2BOHR, 0.0, 0.0],  # 1 Å in Bohr
        [0.0, 1.0 * ANG2BOHR, 0.0],  # 1 Å in Bohr
    ])

    # Distance between atoms 0 and 1 (should be 1 Å)
    d01 = distance_A_from_coords(coords_bohr, 0, 1)
    assert abs(d01 - 1.0) < 1e-10

    # Distance between atoms 0 and 2 (should be 1 Å)
    d02 = distance_A_from_coords(coords_bohr, 0, 2)
    assert abs(d02 - 1.0) < 1e-10

    # Distance between atoms 1 and 2 (should be sqrt(2) Å)
    d12 = distance_A_from_coords(coords_bohr, 1, 2)
    assert abs(d12 - np.sqrt(2.0)) < 1e-10


def test_distance_tag():
    """Test distance_tag formatting utility (integer tag with default digits=2)."""
    from mlmm.utils import distance_tag

    # Default: digits=2, pad=3
    # 1.234 Å → 123 (1.234 × 100 = 123.4 → 123)
    assert distance_tag(1.234) == "123"

    # 2.567 Å → 257 (2.567 × 100 = 256.7 → 257)
    assert distance_tag(2.567) == "257"

    # Zero
    assert distance_tag(0.0) == "000"

    # Large value: 12.345 Å → 1234 (12.345 × 100 = 1234.5 → 1234, Python rounds to even)
    assert distance_tag(12.345) == "1234"

    # Custom digits
    assert distance_tag(1.234, digits=3, pad=4) == "1234"  # 1.234 × 1000


def test_values_from_bounds():
    """Test values_from_bounds grid generation."""
    from mlmm.utils import values_from_bounds

    # Simple case
    values = values_from_bounds(0.0, 1.0, 0.25)
    expected = [0.0, 0.25, 0.5, 0.75, 1.0]
    np.testing.assert_allclose(values, expected, atol=1e-10)

    # Single point
    values = values_from_bounds(1.5, 1.5, 0.1)
    assert len(values) == 1
    assert abs(values[0] - 1.5) < 1e-10

    # Reversed bounds (should work)
    values = values_from_bounds(2.0, 1.0, 0.5)
    assert len(values) == 3


def test_load_yaml_dict():
    """Test YAML loading utility."""
    from mlmm.utils import load_yaml_dict

    # None input
    assert load_yaml_dict(None) == {}

    # Create a temporary YAML file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write("geom:\n  coord_type: cart\ncalc:\n  ml_device: cpu\n")
        yaml_path = f.name

    try:
        config = load_yaml_dict(yaml_path)
        assert "geom" in config
        assert "calc" in config
        assert config["geom"]["coord_type"] == "cart"
        assert config["calc"]["ml_device"] == "cpu"
    finally:
        Path(yaml_path).unlink()


def test_apply_yaml_overrides():
    """Test YAML override application (first matching path is used)."""
    from mlmm.utils import apply_yaml_overrides

    # Base config
    base_cfg = {"a": 1, "b": 2, "c": 3}

    # YAML overrides
    yaml_cfg = {
        "section1": {"a": 10, "d": 4},
    }

    # Apply overrides (only section1 is used)
    apply_yaml_overrides(
        yaml_cfg,
        [(base_cfg, (("section1",),))]
    )

    # Check results
    assert base_cfg["a"] == 10  # Overridden from section1
    assert base_cfg["b"] == 2   # Unchanged
    assert base_cfg["c"] == 3   # Unchanged
    assert base_cfg["d"] == 4   # Added from section1

    # Test multiple paths (first matching is used)
    base_cfg2 = {"x": 1}
    yaml_cfg2 = {"fallback": {"x": 99}, "primary": {"x": 42}}

    apply_yaml_overrides(
        yaml_cfg2,
        [(base_cfg2, (("primary",), ("fallback",)))]  # primary first
    )
    assert base_cfg2["x"] == 42  # primary is used


def test_ensure_dir():
    """Test directory creation utility."""
    from mlmm.utils import ensure_dir

    with tempfile.TemporaryDirectory() as tmpdir:
        test_dir = Path(tmpdir) / "subdir" / "nested"

        # Directory should not exist initially
        assert not test_dir.exists()

        # Create it
        ensure_dir(test_dir)

        # Should exist now
        assert test_dir.exists()
        assert test_dir.is_dir()

        # Should be idempotent
        ensure_dir(test_dir)
        assert test_dir.exists()


def test_pretty_block():
    """Test pretty_block formatting utility."""
    from mlmm.utils import pretty_block

    data = {"key1": "value1", "key2": 123, "key3": True}
    output = pretty_block("test_section", data)

    # Should contain header
    assert "test_section" in output

    # Should contain all keys
    assert "key1" in output
    assert "key2" in output
    assert "key3" in output


def test_pretty_block_with_numpy_scalars():
    """pretty_block should accept NumPy scalar/list inputs without YAML errors."""
    from mlmm.utils import pretty_block

    data = {"freeze_atoms": [np.int64(0), np.int64(2)], "ratio": np.float64(1.5)}
    output = pretty_block("freeze_atoms (effective)", data)

    assert "freeze_atoms" in output
    assert "- 0" in output
    assert "- 2" in output
    assert "ratio: 1.5" in output


def test_resolve_charge_spin_or_raise_requires_charge():
    """Charge remains unresolved by default and should raise."""
    from mlmm.utils import PreparedInputStructure, resolve_charge_spin_or_raise

    prepared = PreparedInputStructure(
        source_path=Path("dummy.pdb"),
        geom_path=Path("dummy.pdb"),
    )
    with pytest.raises(Exception, match="Total charge is unresolved"):
        resolve_charge_spin_or_raise(prepared, charge=None, spin=None)


def test_resolve_charge_spin_or_raise_accepts_explicit_charge():
    """Explicit charge with omitted spin should resolve using spin default."""
    from mlmm.utils import PreparedInputStructure, resolve_charge_spin_or_raise

    prepared = PreparedInputStructure(
        source_path=Path("dummy.pdb"),
        geom_path=Path("dummy.pdb"),
    )
    charge, spin = resolve_charge_spin_or_raise(prepared, charge=-1, spin=None)
    assert charge == -1
    assert spin == 1


def test_merge_freeze_atom_indices():
    """Test freeze atom index merging (requires geom_cfg dict)."""
    from mlmm.utils import merge_freeze_atom_indices

    # Empty geom_cfg
    geom_cfg = {}
    result = merge_freeze_atom_indices(geom_cfg)
    assert result == []
    assert geom_cfg["freeze_atoms"] == []

    # geom_cfg with existing freeze_atoms
    geom_cfg = {"freeze_atoms": [1, 2, 3]}
    result = merge_freeze_atom_indices(geom_cfg, [4, 5, 6])
    assert result == [1, 2, 3, 4, 5, 6]
    assert geom_cfg["freeze_atoms"] == [1, 2, 3, 4, 5, 6]

    # Merge with duplicates
    geom_cfg = {"freeze_atoms": [1, 2, 3]}
    result = merge_freeze_atom_indices(geom_cfg, [2, 3, 4, 5])
    assert result == [1, 2, 3, 4, 5]

    # Merge and sort (unordered input)
    geom_cfg = {"freeze_atoms": [5, 3, 1]}
    result = merge_freeze_atom_indices(geom_cfg, [4, 2])
    assert result == [1, 2, 3, 4, 5]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
