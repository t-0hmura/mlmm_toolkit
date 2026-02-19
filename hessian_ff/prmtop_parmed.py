from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Union


def read_prmtop_with_parmed(prmtop_path: Union[str, Path]) -> Dict[str, List[Any]]:
    """Read a prmtop using ParmEd.

    Returns a dict mapping prmtop %FLAG names to the underlying raw arrays.
    This package requires ParmEd and always uses this route.
    """
    try:
        from parmed.amber import AmberParm  # type: ignore
    except Exception as e:  # pragma: no cover
        raise ImportError(
            "ParmEd is required but not installed. Install it (pip install parmed)."
        ) from e

    prmtop_path = Path(prmtop_path)
    parm = AmberParm(str(prmtop_path))
    # parm.parm_data is a dict: flag -> numpy array/list
    return {k: list(v) for k, v in parm.parm_data.items()}
