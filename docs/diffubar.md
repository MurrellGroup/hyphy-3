# Native difFUBAR

`hyphy3 diffubar` and `hyphy3-diffubar` implement Bayesian differential codon
selection between two branch groups, optionally with a third background regime.
The analysis uses native hyphy-3 alignments, genetic codes, annotated trees,
MG94 rate matrices, optimizers, and Wheatstone pruning. It does not invoke Julia,
Node, a browser, or an external inference process.

## Run

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target hyphy3 hyphy3-diffubar --parallel 2

build/hyphy3 diffubar \
  --alignment alignment.fasta --tree tagged-tree.nwk \
  --groups G1,G2 --threads 8 --seed 572 \
  --output result.difFUBAR.json
```

A tagged Newick tree is required, either through `--tree` or embedded in a
native NEXUS alignment. The analysis does not infer foreground groups or branch
lengths from an alignment. FASTA/NEXUS loading uses the existing `Alignment`
loader, including its existing format limitations. Use `--help` for all options.

### Branch annotations

Tags are read from `TreeNode::model_tag`, not by substring matching taxon names.
For example:

```text
((a{G1}:0.1,b{G1}:0.1){G1}:0.05,(c{G2}:0.1,d{G2}:0.1){G2}:0.05,e:0.2);
```

A tag labels the **incoming branch** of that node. Internal tags do not
implicitly propagate to descendants. Here, `e` is a background branch.
Without `--groups`, exactly two distinct non-root tags must exist; their lexical
order defines group 1 and group 2. With `--groups G1,G2`, that exact order is used
and every other branch, including branches with other tags, is background.
Tag matching is case-sensitive. The root has no incoming edge and its annotation
is ignored. An untagged root alone does not create a background dimension.
Colored FigTree comments are not interpreted as groups: convert them to the
native `{tag}` representation first.

C++ trees containing native Wheatstone bridges are supported. Each bridge's
three resolution lengths share the regime of its original `bridge_node_id`,
not that of the hub. Fixed mixture weights and bridge lengths are respected in
both global fitting and grid evaluation. The original non-root branch must be
retained in `tree.nodes`; the bridge is evaluated through the native hub API.
The CLI accepts the native Newick loader's inputs; it does not add a separate
network-file importer. JSON records bridge IDs, weights, and lengths for auditing.
A Newick string by itself does not serialize all Wheatstone state.

## Model and reference compatibility

At each site the grid has a shared synonymous rate `a`, foreground ratios
`omega1` and `omega2`, and, only when background branches exist, `omega_background`.
For codons differing at one position, the off-diagonal MG94 rate is

```text
q_ij = a * r_xy * f_position,y                   (synonymous)
q_ij = a * r_xy * f_position,y * omega_group     (nonsynonymous)
```

Other off-diagonal rates are zero. The root distribution is the product of the
three positional frequencies, normalized over sense codons under the selected
native genetic code.

The default global fit has two stages on the **fixed input branch lengths**:
(1) estimate six absolute nucleotide exchangeabilities using the native GTR
likelihood; (2) fix those rates and estimate global codon alpha and beta.
The native bounded optimizers are used in log parameters, with bounds 1e-4 to
1e4. Two coordinate-polishing passes follow the nucleotide simplex fit.
The grid then operates on the equivalent of input branch lengths multiplied
by fitted alpha. Individual input branches are never refitted or overwritten.

The GTR fit restores `GTRMatrix::scale_factor` because that native matrix is
normalized. MG94 is not normalized. Its native AG-reference parameterization
is reconciled with absolute reference exchangeabilities by using `r/r_AG` and
multiplying grid branch time by `global_alpha * r_AG * a`. This avoids an
otherwise silent scale mismatch with Julia/EvoOnline.

### F3x4 is intentionally not CF3x4

The reference analysis uses raw positional nucleotide counts, whereas native
`Alignment::pos_nuc_frequencies` contains corrected CF3x4 estimates. difFUBAR
recomputes raw F3x4 from the alignment without changing the alignment or the
behavior of other analyses. The denominator is taxa times codon sites, including
missing bases, as in the reference `count_F3x4`; rows therefore need not sum to
one when data are missing. Native U-to-T conversion is retained.

Zero empirical entries are allowed without pseudocounts. The native MG94
generator is diagonalized on its closed positive-equilibrium support when the
core's frequency floor would otherwise alter the spectral problem. Frequencies
that give every sense codon zero equilibrium mass are rejected. Global nucleotide
fitting currently requires all four bases represented somewhere in the alignment;
a supplied fixed model can be used for such degenerate datasets.

### Grid and posterior defaults

| Setting | Default |
| --- | --- |
| Foreground density | 6 points below one; 12 values per alpha/foreground axis |
| Background density | 4 points below one; 7 background values |
| Categories | 1,728 without background; 12,096 with background |
| Dirichlet concentration | 0.1 **per category**, not total concentration |
| Gibbs iterations | 2,500 per chain |
| Burn-in | floor(iterations / 5) |
| Chains | 1 |
| Posterior detection threshold | 0.95 |
| Seed | 1 |

Axes use the reference log10 transform with offset 0.05. The mathematically
neutral grid point is stored as exactly 1, preventing floating-point rounding
from classifying neutrality as `omega > 1`. Alpha varies slowest, followed by
omega1, omega2, and background fastest. Category tuples are exported explicitly.

Inference is an uncollapsed categorical-allocation/Dirichlet Gibbs sampler.
Repeated alignment patterns reuse likelihoods, but **each original site receives
an independent latent allocation**. Pattern weights are not treated as one
powered observation. Reported probabilities and marginal means are
Rao-Blackwellized averages of normalized site posteriors over retained theta
draws, not posteriors evaluated at the mean theta. This targets the same posterior
as allocation-count summaries, with different Monte Carlo noise. No arbitrary
conditional-likelihood cutoff is applied; ordinary double-precision underflow
still exists. Impossible patterns cause an error rather than a uniform fallback.

The four principal probabilities are `P(omega1 > omega2)`, `P(omega2 > omega1)`,
`P(omega1 > 1)`, and `P(omega2 > 1)`. Equality mass is reported separately, so the
two directional probabilities are not complements. JSON includes one-based site
indices, posterior means, pattern-level marginal distributions, grid weights,
branch assignments, fitted parameters, and per-chain probability summaries.
It also provides an MLE-style results table. A difFUBAR-specific HTML viewer is
not included and compatibility with existing report viewers is not asserted.

### Numerical and statistical limits

The default two-stage model and grid match the reference structure, not its
optimizer trajectory or random-number stream. Global optimizers, boundary
solutions, and native treatment of ambiguous/stop codons can change results.
The native representation treats a negative codon state as fully missing, not
as a partial IUPAC codon likelihood. The native genetic-code registry determines
available codes. The supplied Julia comparison harness is universal-code only.

Increasing grid density changes both the discretization and the total prior
concentration when per-category concentration is fixed. Do not compare such
runs as though only numerical resolution changed. Increase iterations/chains
when assessing borderline sites; no ESS, R-hat, or formal convergence claim is
made. Per-chain summaries and a boundary flag are provided for inspection.
Seed reproducibility applies within a build/standard-library implementation,
not across Julia, JavaScript, or different C++ standard libraries.

The dominant memory allocation is `patterns * categories` doubles, plus a second
conditional-likelihood matrix in sampling. Branch transition matrices are cached
one alpha slice at a time. `--max-memory-mb` defaults to 1024 and guards dominant
allocations; it is not a hard process-RSS cap. Likelihood work scales approximately
as categories times patterns times branches, with codon-state matrix products.
The implementation does not yet add Julia's advanced subtree-caching schemes.

## Fixed-model interchange and deterministic parity

```sh
build/hyphy3-diffubar \
  --alignment tests/data/diffubar/demo.fasta \
  --tree tests/data/diffubar/demo.nwk \
  --fitted-model tests/data/diffubar/fixed-model.json \
  --foreground-grid 2 --background-grid 1 --threads 1 \
  --save-fitted-model fitted.json --grid-likelihoods grid.tsv \
  --output result.json
```

`--fitted-model` accepts the versioned model document in the fixture or a prior
result JSON containing `fitted_model`. The six rates are absolute, ordered
AC, AG, AT, CG, CT, GT; F3x4 columns are A, C, G, T. The alpha belongs to the
**original input tree**, not an already alpha-rescaled tree. Double-scaling the
tree is an error in the comparison setup. Input/output path collisions are
rejected. A fixed model is caller-supplied scientific input, not automatically
validated as the optimum for a new alignment.

`--grid-likelihoods` writes unnormalized log likelihoods with categories as rows
and original alignment sites as columns. This separates likelihood-engine
comparisons from optimizer and MCMC differences.

EvoOnline provides `models/diffubar/scripts/generate-julia-reference.jl`. After
running it on the same alignment/tree, compare its exports using:

```sh
python scripts/compare_diffubar_julia.py \
  --binary build/hyphy3-diffubar \
  --alignment alignment.fasta --tree original-tagged-tree.nwk \
  --reference-prefix /path/to/julia-reference \
  --foreground-grid 6 --background-grid 4 --groups G1,G2
```

The harness imports Julia's fitted parameters and matches category tuples
explicitly before comparing raw likelihoods. It does **not** establish native
optimizer parity or MCMC agreement. A Julia run that adds background solely
because of an untagged root will intentionally have a different grid dimension.

## Python and C++

After building the optional native Python module:

```python
import hyphy3

settings = hyphy3.DifFUBARSettings()
settings.groups = ["G1", "G2"]
settings.threads = 8
sampler = hyphy3.DifFUBARSamplerSettings()
sampler.seed = 572
sampler.chains = 2
settings.sampler = sampler

result = hyphy3.diffubar(
    hyphy3.Tree.from_newick_file("tagged-tree.nwk"),
    hyphy3.Alignment.load("alignment.fasta"),
    settings,
)
```

The binding returns a dictionary and releases the GIL during computation. An
optional fourth argument is a fitted-model JSON string. The Python wrapper's
alignment loader currently defaults to the universal code; use the C++ API/CLI
for other native genetic codes unless the existing binding is extended.

```cpp
#include "hyphy/analyses/diffubar.hpp"

using namespace hyphy;
analyses::DifFUBARSettings settings;
settings.groups = {"G1", "G2"};
auto analyzer = analyses::DifFUBARAnalyzer::create_and_fit(
    core::Tree::from_newick_file("tagged-tree.nwk"),
    core::Alignment::load("alignment.fasta"), settings);
auto result = analyzer.run();
```

`use_fitted_model` and `compute_grid_log_likelihoods` expose deterministic stages.
Construct a new analyzer after changing topology, taxon identities, or group
selection; its lookup tables are prepared at construction. Thread/sampler
settings can be adjusted between runs. Native objects are owned by value, so
fitting does not mutate the caller's tree or alignment.

## Tests and validation status

```sh
cmake --build build --target hyphy3-diffubar \
  hyphy_diffubar_inference_test hyphy_diffubar_test hyphy_diffubar_json_test
ctest --test-dir build -R '^DifFUBAR' --output-on-failure
```

When nanobind is available, also build `_hyphy3` before running its Python test.

The standalone numerical and native-engine tests were compiled and executed
with warnings treated as errors against byte-verified upstream core headers.
The native test checks group exchange, absolute-rate scaling, one/two-thread
agreement, root/background semantics, zero frequencies, validation failures,
Wheatstone resolution/mixture equivalence, and a two-stage global-fit smoke test.
An independent dense `scipy.linalg.expm` fixture agrees to maximum absolute
log-likelihood error 2.321e-11 over 128 categories by 12 codon sites. Its generator
is `scripts/generate_diffubar_reference.py`; it does not use native likelihoods.
The sampler is checked against exact enumeration and for correct repeated-site
multiplicities.

The full repository build, JSON/CLI/Python integration tests, and direct
Julia/EvoOnline end-to-end numerical comparisons were **not run in the authoring
environment**. Those are explicit review/CI requirements, not claimed successes.
This patch should not be described as having established broad cross-language
result parity until those comparisons are performed.

## Reference sources

- CodonMolecularEvolution.jl, `src/difFUBAR/difFUBAR.jl`:
  https://github.com/MurrellGroup/CodonMolecularEvolution.jl/blob/main/src/difFUBAR/difFUBAR.jl
- EvoOnline difFUBAR parity notes:
  https://github.com/MurrellGroup/EvoOnline/blob/main/models/diffubar/docs/PARITY.md
- EvoOnline grid and raw-frequency model:
  https://github.com/MurrellGroup/EvoOnline/tree/main/models/diffubar/src/model
- Native annotation, bridge, and rate APIs:
  https://github.com/MurrellGroup/hyphy-3/tree/main/include/hyphy/core
