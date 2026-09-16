"""Read actual solver status through the production adapter and result builder."""
import pytest

from mlmm.workflows.path_opt import _dmf_solver_outcome, _build_dmf_result_data, DMFMepResult


@pytest.mark.parametrize("status,converged", [(0,True),(1,True),(-1,False),(None,False),("bad",False)])
def test_solver_status_drives_public_dmf_result(status, converged):
    actual, parsed, reason = _dmf_solver_outcome(([], {"status":status}))
    assert actual is converged
    result = DMFMepResult(images=(), energies=(-1.,-.5,-1.2), hei_idx=1,
                          converged=actual, ipopt_status=parsed, reason=reason)
    payload = _build_dmf_result_data(result, {"model_charge":0,"model_mult":1})
    assert payload["converged"] is converged
    assert payload["status"] == ("converged" if converged else "not_converged")
