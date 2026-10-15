#!/usr/bin/env python3
"""Run fixed-model native likelihood parity against EvoOnline's Julia exports.

First run EvoOnline/models/diffubar/scripts/generate-julia-reference.jl to
produce PREFIX.model.tsv, PREFIX.categories.tsv and PREFIX.log-likelihoods.tsv.
Then invoke this script with the SAME ORIGINAL, unscaled tagged Newick tree.
NumPy is required; Julia is only needed to produce the reference exports.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True, help="hyphy3-diffubar executable")
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--tree", type=Path, required=True, help="original, unscaled tagged Newick")
    parser.add_argument("--reference-prefix", type=Path, required=True)
    parser.add_argument("--foreground-grid", type=int, default=6)
    parser.add_argument("--background-grid", type=int, default=4)
    parser.add_argument("--tolerance", type=float, default=1e-8)
    parser.add_argument("--groups", help="two comma-separated tags, in reference order")
    args = parser.parse_args()
    if not np.isfinite(args.tolerance) or args.tolerance <= 0:
        parser.error("--tolerance must be finite and positive")
    prefix = str(args.reference_prefix)
    values = np.loadtxt(prefix + ".model.tsv").reshape(-1)
    # The supplied Julia generator uses universal_code and exports 61 equilibrium entries.
    if len(values) != 3 + 6 + 12 + 61 or not np.isfinite(values).all():
        raise ValueError("Expected a universal-code Julia model row: alpha,beta,LL,6 rates,12 F3x4,61 equilibrium")
    model = {
        "schema": "hyphy3-diffubar-fitted-model-v1", "tree_branch_units": "input-tree",
        "exchangeability_order": ["AC", "AG", "AT", "CG", "CT", "GT"],
        "nucleotide_order": ["A", "C", "G", "T"], "genetic_code": "Universal",
        "alpha": float(values[0]), "beta": float(values[1]), "codon_log_likelihood": float(values[2]),
        "exchangeabilities": values[3:9].tolist(), "f3x4": values[9:21].reshape(3, 4).tolist(),
    }
    reference_categories = np.loadtxt(prefix + ".categories.tsv", ndmin=2)
    reference_logs = np.loadtxt(prefix + ".log-likelihoods.tsv", ndmin=2)
    if reference_categories.shape[0] != reference_logs.shape[0]:
        raise ValueError("Reference category and likelihood dimensions disagree")
    with tempfile.TemporaryDirectory(prefix="diffubar-parity-") as temporary:
        temp = Path(temporary)
        (temp / "model.json").write_text(json.dumps(model, allow_nan=False))
        command = [str(args.binary.resolve()), "--alignment", str(args.alignment.resolve()),
                   "--tree", str(args.tree.resolve()), "--fitted-model", str(temp / "model.json"),
                   "--foreground-grid", str(args.foreground_grid), "--background-grid", str(args.background_grid),
                   "--iterations", "2", "--burn-in", "1", "--seed", "1", "--threads", "1",
                   "--output", str(temp / "result.json"), "--grid-likelihoods", str(temp / "native.tsv"), "--quiet"]
        if args.groups:
            command += ["--groups", args.groups]
        subprocess.run(command, check=True)
        native = np.loadtxt(temp / "native.tsv", ndmin=2)
        result = json.loads((temp / "result.json").read_text())
        categories = np.asarray(result["grid"]["categories"], dtype=float)
    if categories.shape != reference_categories.shape or native.shape != reference_logs.shape:
        raise ValueError("Dimension mismatch; check grid densities, group ordering, and root/background conventions")
    # Julia's Cartesian-product order can differ; explicitly match category tuples.
    order = []
    for category in reference_categories:
        matches = np.flatnonzero(np.all(np.isclose(categories, category, rtol=1e-12, atol=1e-12), axis=1))
        if len(matches) != 1:
            raise ValueError(f"No unique native match for Julia category {category.tolist()}")
        order.append(int(matches[0]))
    reordered = native[order, :]
    if not np.array_equal(np.isneginf(reordered), np.isneginf(reference_logs)):
        raise AssertionError("Impossible-pattern/category masks differ")
    valid = np.isfinite(reference_logs) & np.isfinite(reordered)
    if np.isnan(reference_logs).any() or np.isnan(reordered).any() or not valid.any():
        raise AssertionError("Invalid or entirely non-finite likelihood grid")
    error = np.abs(reordered[valid] - reference_logs[valid])
    worst = float(error.max())
    print(f"Fixed-model Julia/native grid comparison: max |delta lnL| = {worst:.12g}; "
          f"{reordered.shape[0]} categories x {reordered.shape[1]} sites")
    if worst > args.tolerance:
        raise AssertionError(f"Likelihood mismatch exceeds tolerance {args.tolerance:g}")
    print("PASS. This comparison does not test the native global optimizer or cross-language MCMC streams.")

if __name__ == "__main__":
    main()
