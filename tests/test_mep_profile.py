"""Fast local-peak and preceding-basin contracts, without model imports."""
import pytest
from pysisyphus.constants import AU2KCALPERMOL


@pytest.mark.parametrize("energies,peak,basin,barrier", [
    ([4.0, 0.0, 2.0, 1.0, 5.0], 2, 1, 2.0),
    ([0.0, 2.0, 1.0, 4.0, 3.0], 3, 2, 3.0),
    ([0.0, 2.0, 2.0, 1.0, 3.0], 1, 0, 2.0),
    ([0.0, 1.0, 1.0, 2.0], 3, 0, None),
    ([0.0, 2.0, 2.0], 2, 0, None),
    ([2.0, 1.0, 0.0], 0, 0, None),
])
def test_mep_barrier_uses_the_basin_before_the_selected_peak(energies, peak, basin, barrier):
    from mlmm.domain.mep_profile import mep_barrier_kcal, select_mep_peak

    assert select_mep_peak(energies) == (peak, basin)
    if barrier is None:
        assert mep_barrier_kcal(energies) is None
    else:
        assert mep_barrier_kcal(energies) == pytest.approx(barrier * AU2KCALPERMOL)
