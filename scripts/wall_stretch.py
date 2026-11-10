#!/usr/bin/env python3
"""Wall-normal element boundaries for case/channel.box.

The mesh runs from y_lo + z1 to y_hi - z1 (z1 = distance from each virtual
wall to its mesh boundary). In the lower half, element heights grow
geometrically by `ratio` from the boundary to the midplane (the first height
is chosen so the sum closes exactly on the midplane); the upper half is the
mirror image. Prints the boundaries (12 significant digits, wrapped for
genbox's 132-character line limit) and the GLL node heights in wall units.

    python3 scripts/wall_stretch.py --z1 0.0192307692307692 --ne-half 10 --ratio 1.15
"""
from __future__ import annotations

import argparse

import numpy as np


def half_boundaries(z1, y_lo, y_hi, ne_half, ratio):
    mid = 0.5 * (y_lo + y_hi)
    length = mid - (y_lo + z1)
    if ratio == 1.0:
        steps = np.full(ne_half, length / ne_half)
    else:
        h1 = length * (ratio - 1.0) / (ratio ** ne_half - 1.0)
        steps = h1 * ratio ** np.arange(ne_half)
    lo = (y_lo + z1) + np.concatenate([[0.0], np.cumsum(steps)])
    lo[-1] = mid  # close exactly on the midplane
    hi = (y_lo + y_hi) - lo[::-1][1:]
    return np.concatenate([lo, hi])


def gll_nodes(order):
    from numpy.polynomial import legendre as L
    c = np.zeros(order + 1)
    c[-1] = 1.0
    interior = np.real(L.legroots(L.legder(c)))
    return np.concatenate([[-1.0], np.sort(interior), [1.0]])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--z1", type=float, required=True)
    ap.add_argument("--ne-half", type=int, required=True)
    ap.add_argument("--ratio", type=float, required=True)
    ap.add_argument("--y-lo", type=float, default=0.0)
    ap.add_argument("--y-hi", type=float, default=2.0)
    ap.add_argument("--order", type=int, default=7)
    ap.add_argument("--re-tau", type=float, default=5200.0)
    ap.add_argument("--nodes", type=int, default=2, help="lower-half elements to tabulate")
    a = ap.parse_args()
    b = half_boundaries(a.z1, a.y_lo, a.y_hi, a.ne_half, a.ratio)
    vals = [f"{v:.12g}" for v in b]
    print(f"# .box: nely = {len(b) - 1} (positive: explicit boundaries follow)")
    for i in range(0, len(vals), 6):
        print("  ".join(vals[i:i + 6]))
    g = gll_nodes(a.order)
    h = 0.5 * (a.y_hi - a.y_lo)
    print("\n# lower-half GLL nodes: element, node, y, distance from virtual wall, y+ (u_tau = 1)")
    for e in range(min(a.nodes, a.ne_half)):
        for k, xi in enumerate(g):
            y = b[e] + 0.5 * (xi + 1.0) * (b[e + 1] - b[e])
            d = y - a.y_lo
            print(f"# e={e} k={k} y={y:.9f} dist={d:.9f} y+={d * a.re_tau / h:.3f}")


if __name__ == "__main__":
    main()
