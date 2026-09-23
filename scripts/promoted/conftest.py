"""pytest bootstrap: put both ``scripts/scratch`` (the ``loaders`` package) and
``scripts/promoted`` (the existing ``common`` package: ``common.hashing``,
``common.design``) on ``sys.path`` so tests can import either without a name
collision (see ``loaders/__init__.py`` for the scheme).
"""

from __future__ import annotations

import sys
from pathlib import Path

_SCRATCH = Path(__file__).resolve().parent
_PROMOTED = _SCRATCH.parent / "promoted"

for _p in (str(_SCRATCH), str(_PROMOTED)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
