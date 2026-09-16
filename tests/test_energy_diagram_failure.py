"""Optional rendering must preserve the computed energy payload."""
import pytest
from mlmm.workflows import all as workflow


@pytest.mark.parametrize("failure", [None, "build", "layout", "write"])
def test_rendering_failure_retains_numeric_results(tmp_path, monkeypatch, failure):
    def check(stage):
        if failure == stage:
            raise RuntimeError(f"{stage} unavailable")

    class Figure:
        def update_layout(self, **kwargs):
            check("layout")

    def build(**kwargs):
        check("build")
        return Figure()

    def write(fig, path, **kwargs):
        check("write")
        path.write_bytes(b"image fixture")

    monkeypatch.setattr(workflow, "build_energy_diagram", build)
    monkeypatch.setattr(workflow, "write_plotly_image", write)
    energies = [-1.0, -0.9, -1.1]
    (tmp_path / "energy.png").write_bytes(b"stale prior image")
    payload = workflow._write_segment_energy_diagram(
        tmp_path / "energy", ["R", "TS", "P"], energies, "test",
    )
    assert payload["energies_au"] == energies
    assert payload["energies_kcal"] == pytest.approx(
        [0, .1 * workflow.AU2KCALPERMOL, -.1 * workflow.AU2KCALPERMOL]
    )
    assert payload["image_written"] is (failure is None)
    if failure:
        assert payload["image"] is None
        assert payload["image_error"] == f"{failure} unavailable"
    else:
        assert payload["image"] == str(tmp_path / "energy.png")
        assert "image_error" not in payload
