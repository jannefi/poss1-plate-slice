from __future__ import annotations
from typing import Dict


def init_buckets() -> Dict[str, int]:
    return {
        'spikes_rejected': 0,
        'morphology_rejected': 0,
        'hpm_objects': 0,
    }


def finalize(b: Dict[str, int]) -> Dict[str, int]:
    # ensure counts are present; could add sanity checks here
    return b
