#!/usr/bin/env python3
"""Turn the wall-model diagnostics in a nekRS log into a CSV time series.

Reads the `[logbnd]` and `[logbnd-wall]` lines that case/channel.udf prints
every diag_every steps and writes one row per print:

  step, t, Ub          bulk velocity over the simulated region (mesh only)
  maxU                 max |u| over all nodes
  utau1, utau2         unweighted mean of the per-node u_tau, top / bottom wall
  utauF1, utauF2       sqrt(-Fx/A): u_tau from the area-integrated delivered traction
  budget_res           H dUb/dt - (G H + (Fx1 + Fx2)/A), with H = mesh height,
                       G = dpdx; dUb/dt by central differences (NaN at the ends)

    python3 scripts/extract_log.py run.log [--out timeseries.csv] [--dpdx 1]
"""
from __future__ import annotations

import argparse
import csv
import math
import re

KV = re.compile(r"(\w+)=(\S+)")


def kv(line):
    out = {}
    for k, v in KV.findall(line):
        try:
            out[k] = float(v)
        except ValueError:
            out[k] = v
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log")
    ap.add_argument("--out", default="timeseries.csv")
    ap.add_argument("--dpdx", type=float, default=1.0)
    a = ap.parse_args()

    steps, walls, setup = {}, {}, {}
    with open(a.log, errors="replace") as f:
        for ln in f:
            if ln.startswith("logbnd-setup:"):
                setup = kv(ln)
            elif ln.startswith("[logbnd-wall]"):
                d = kv(ln)
                walls[(int(d["step"]), int(d["wall"]))] = d
            elif ln.startswith("[logbnd]"):
                d = kv(ln)
                steps[int(d["step"])] = d
    if not steps:
        raise SystemExit("no [logbnd] lines found")
    H = setup.get("meshYHi", float("nan")) - setup.get("meshYLo", float("nan"))

    rows = []
    for s in sorted(steps):
        w1, w2 = walls.get((s, 1)), walls.get((s, 2))
        if w1 is None or w2 is None:
            continue
        rows.append({
            "step": s, "t": steps[s]["t"], "Ub": steps[s]["Ubulk_res"], "maxU": steps[s]["maxU"],
            "utau1": w1["uTau_mean"], "utau2": w2["uTau_mean"],
            "utauF1": math.sqrt(max(-w1["Fx"] / w1["A"], 0.0)),
            "utauF2": math.sqrt(max(-w2["Fx"] / w2["A"], 0.0)),
            "Fx1": w1["Fx"], "Fx2": w2["Fx"], "A": w1["A"],
        })
    for i, r in enumerate(rows):
        if 0 < i < len(rows) - 1:
            dudt = (rows[i + 1]["Ub"] - rows[i - 1]["Ub"]) / (rows[i + 1]["t"] - rows[i - 1]["t"])
            r["budget_res"] = H * dudt - (a.dpdx * H + (r["Fx1"] + r["Fx2"]) / r["A"])
        else:
            r["budget_res"] = float("nan")
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    last = rows[-1]
    print(f"{len(rows)} rows -> {a.out}; last t={last['t']:.4f} Ub={last['Ub']:.4f} "
          f"utau(mean)={0.5 * (last['utau1'] + last['utau2']):.4f} "
          f"utau(traction)={0.5 * (last['utauF1'] + last['utauF2']):.4f}")


if __name__ == "__main__":
    main()
