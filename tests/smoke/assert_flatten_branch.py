#!/usr/bin/env python3
"""Require execution of the higher-order-saddle flatten branch."""

from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path


log_path = Path(sys.argv[1])
result_path = Path(sys.argv[2])
text = log_path.read_text(encoding="utf-8")

initial_match = re.search(r"^\[Imaginary modes\] n=(\d+)", text, flags=re.MULTILINE)
if initial_match is None:
    raise SystemExit("flatten smoke omitted the initial imaginary-mode count")
initial = int(initial_match.group(1))
if initial <= 1:
    raise SystemExit(f"flatten fixture is not a higher-order candidate: n_imag={initial}")
if "[flatten] Extra imaginary modes detected; starting RS-I-RFO flatten loop." not in text:
    raise SystemExit("flatten branch did not start")
if "[flatten] RS-I-RFO iteration 1/1" not in text:
    raise SystemExit("flatten branch did not execute its requested iteration")
if "skipping flatten loop" in text.lower() or "No eligible modes to flatten" in text:
    raise SystemExit("flatten branch was skipped after entry")

post_counts = [
    int(value)
    for value in re.findall(
        r"^\[Imaginary modes:(?:primary|alternate)\] n=(\d+)",
        text,
        flags=re.MULTILINE,
    )
]
if not post_counts:
    raise SystemExit("flatten branch omitted post-iteration imaginary-mode counts")
if min(post_counts) > initial:
    raise SystemExit(
        f"flatten branch increased saddle order: initial={initial}, post={post_counts}"
    )

payload = json.loads(result_path.read_text(encoding="utf-8"))
final_count = payload.get("n_imaginary_modes")
if final_count is None or not math.isfinite(float(final_count)):
    raise SystemExit(f"flatten result omitted final saddle order: {final_count!r}")
if int(final_count) > initial:
    raise SystemExit(
        f"flatten result increased saddle order: initial={initial}, final={final_count}"
    )
