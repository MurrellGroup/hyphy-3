#!/usr/bin/env python3
"""Regenerate independent fixed-model likelihood fixture (requires NumPy/SciPy).

Uses EvoOnline's six-taxon demo, a declared non-fitted parameter fixture, an
independent MG94 generator, dense scipy.linalg.expm, and direct pruning. No
hyphy-3 code, eigenvectors, probabilities, or fitted output are used.
"""
from pathlib import Path
from itertools import product
import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests/data/diffubar"
RATES = np.array([[0, 1.2, .9, .7], [1.2, 0, 1.1, 2.1], [.9, 1.1, 0, .8], [.7, 2.1, .8, 0]])
F3X4 = np.array([[.3, .2, .2, .3], [.2, .3, .3, .2], [.25, .25, .25, .25]])
GLOBAL_ALPHA = .8
BASES = "ACGT"
# Universal translation in the lexically ordered 64-codon state space.
AA64 = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"
CODONS = ["".join(c) for c in product(BASES, repeat=3)]
TRANSLATION = dict(zip(CODONS, AA64, strict=True))
CODONS = [c for c in CODONS if TRANSLATION[c] != "*"]
INDEX = {c: i for i, c in enumerate(CODONS)}
BASE_INDEX = {c: i for i, c in enumerate(BASES)}
PI = np.array([np.prod([F3X4[p, BASE_INDEX[c[p]]] for p in range(3)]) for c in CODONS])
PI /= PI.sum()


def axis(lower: float, upper: float, below: int) -> np.ndarray:
    start = np.log10(lower + .05)
    step = (np.log10(1.05) - start) / below
    n = int(np.floor((np.log10(upper + .05) - start) / step + 1e-12)) + 1
    result = 10 ** (start + step * np.arange(n)) - .05
    result[below] = 1
    return result


def generator(alpha: float, omega: float) -> np.ndarray:
    q = np.zeros((len(CODONS), len(CODONS)))
    for i, source in enumerate(CODONS):
        for j, target in enumerate(CODONS):
            positions = [k for k in range(3) if source[k] != target[k]]
            if len(positions) != 1:
                continue
            k = positions[0]
            a, b = BASE_INDEX[source[k]], BASE_INDEX[target[k]]
            q[i, j] = alpha * RATES[a, b] * F3X4[k, b]
            if TRANSLATION[source] != TRANSLATION[target]:
                q[i, j] *= omega
    q[np.diag_indices_from(q)] = -q.sum(axis=1)
    return q


def main() -> None:
    # Explicit topology/lengths match demo.nwk; no shared parser/pruning code.
    names = ("alpha", "beta", "gamma", "delta", "epsilon", "zeta")
    sequences: dict[str, str] = {}
    name = ""
    for line in (DATA / "demo.fasta").read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]; sequences[name] = ""
        else:
            sequences[name] += line
    sites = len(sequences[names[0]]) // 3
    tips = {}
    for name in names:
        tips[name] = np.zeros((len(CODONS), sites))
        for s in range(sites):
            tips[name][INDEX[sequences[name][3*s:3*s+3]], s] = 1
    alpha_axis, omega_axis, background = axis(.01, 13, 2), axis(.01, 13, 2), axis(.05, 6, 1)
    output = []
    for alpha in alpha_axis:
        # Cache only transition matrices; this is still dense expm, not spectral pruning.
        cache = {}
        for omega in set(omega_axis) | set(background):
            q = generator(GLOBAL_ALPHA * alpha, omega)
            for length in (.03, .04, .05, .06):
                cache[omega, length] = expm(q * length)
        for w1, w2, bg in product(omega_axis, omega_axis, background):
            one = cache[w1, .04] @ ((cache[w1, .05] @ tips["alpha"]) * (cache[w1, .05] @ tips["beta"]))
            two = cache[w2, .04] @ ((cache[w2, .05] @ tips["gamma"]) * (cache[w2, .05] @ tips["delta"]))
            three = cache[bg, .03] @ ((cache[bg, .06] @ tips["epsilon"]) * (cache[bg, .06] @ tips["zeta"]))
            output.append(np.log(PI @ (one * two * three)))
    values = np.asarray(output)
    if values.shape != (128, sites) or not np.isfinite(values).all():
        raise RuntimeError("Invalid independently generated reference")
    np.savetxt(DATA / "scipy-log-likelihoods.tsv", values, delimiter="\t", fmt="%.17g")
    print(f"Generated independent SciPy fixture: {values.shape[0]} categories x {sites} sites")

if __name__ == "__main__":
    main()
