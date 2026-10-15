# difFUBAR test fixtures

`demo.fasta` and `demo.nwk` are unchanged copies of the six-taxon example in
MurrellGroup/EvoOnline, `examples/diffubar-demo.fasta` and
`examples/diffubar-demo.nwk`, read from main during implementation.

Upstream Git blob IDs:

- FASTA: `3d609cb46b658868a91eccfa218c0608f59b4255`
- Newick: `ffd2c3f51b42c59bf0a755d6a233c6360eefe3b7`

This is an example fixture, not a claimed biological validation dataset.

`fixed-model.json` declares synthetic test parameters, not fitted estimates.
`scipy-log-likelihoods.tsv` is an independently generated dense-matrix-exponential
reference for that model and example. It has 128 category rows and 12 site
columns, with foreground density 2/background density 1 and the category order
specified in `docs/diffubar.md`.

Regenerate using `python scripts/generate_diffubar_reference.py` (NumPy/SciPy).
The generator independently constructs the universal-code MG94 generator and
prunes the explicit fixture topology. No native likelihood output is an input.

Source: https://github.com/MurrellGroup/EvoOnline/tree/main/examples
