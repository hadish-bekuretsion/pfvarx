"""PfVarX — Plasmodium falciparum SNP/INDEL analysis platform.

RESEARCH USE ONLY. See ``pfvarx.release`` for the machine-readable release
status of every subsystem. No component of this package is certified for
clinical decision-making.
"""

from pfvarx.release import (
    RELEASE_STATUS,
    get_status,
    overall_status,
    status_banner,
)

__version__ = "0.9.0"
__all__ = [
    "__version__",
    "RELEASE_STATUS",
    "get_status",
    "overall_status",
    "status_banner",
]
