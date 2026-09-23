"""Project loader/QC package for the raloxifene/HLM LFQ experiment (Stage 3).

Lives at ``scripts/promoted/loaders/`` (promotion: ``git mv`` to
``scripts/promoted/loaders/``, unchanged import path) rather than as new modules
inside ``scripts/promoted/common/`` because pre-review code may not be added to an
already-promoted package. Entry scripts put **both** ``scripts/scratch`` (for this
``loaders`` package) and ``scripts/promoted`` (for the existing ``common`` package:
``common.hashing``, ``common.design``) on ``sys.path``, so the two packages coexist
under distinct top-level names with no import collision.
"""

from __future__ import annotations
