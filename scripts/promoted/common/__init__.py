"""Project-shared modules used by both ``metadata_characterize.py`` and
``metadata_figures.py`` (and any later Stage-1+ script): hashing/data-version
helpers (``common.hashing``), the shared Stage-1 design rules both scripts must
agree on (``common.design``), and the figure machinery (``common.figures``).

This package is a plain sibling of the scripts that import it: each script does
``sys.path.insert(0, str(Path(__file__).resolve().parent))`` then
``import common.xxx``, so promoting a script is a pure
``git mv scripts/scratch/X scripts/promoted/X`` (script + this ``common/``
directory move together; no import needs to change).
"""
