"""Row assignment for overlapping features."""
from __future__ import annotations

import numpy as np


def pack_rows(starts, ends, gap: float = 0.0) -> np.ndarray:
    """Assign intervals to rows so that none overlap within a row.

    Greedy first fit over intervals sorted by start: each goes to the top-most
    row it fits in.

    Parameters
    ----------
    starts, ends : array-like
    gap : float
        Minimum distance between neighbours in one row.

    Returns
    -------
    numpy.ndarray of int
        Row index per interval, in input order.
    """
    starts = np.asarray(starts, dtype=float)
    ends = np.asarray(ends, dtype=float)
    rows = np.zeros(len(starts), dtype=int)
    row_ends: list = []
    for i in np.argsort(starts, kind="stable"):
        for r, row_end in enumerate(row_ends):
            if starts[i] >= row_end + gap:
                rows[i] = r
                row_ends[r] = ends[i]
                break
        else:
            rows[i] = len(row_ends)
            row_ends.append(ends[i])
    return rows
