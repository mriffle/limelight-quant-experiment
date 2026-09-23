"""pytest bootstrap for scratch tests: put ``scripts/scratch`` (the ``analysis``
package and runner scripts) and ``scripts/promoted`` (``loaders``, ``common``) on
``sys.path``.
"""

from __future__ import annotations

import sys
from pathlib import Path

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"

for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
