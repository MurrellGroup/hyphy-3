"""Native Python binding smoke/replay test; run after building _hyphy3."""
from pathlib import Path
import json
import hyphy3


def main() -> None:
    data = Path(__file__).resolve().parents[2] / "tests/data/diffubar"
    settings = hyphy3.DifFUBARSettings()
    settings.foreground_grid = 2
    settings.background_grid = 1
    settings.threads = 1
    sampler = hyphy3.DifFUBARSamplerSettings()
    sampler.iterations = 80
    sampler.seed = 572
    settings.sampler = sampler
    tree = hyphy3.Tree.from_newick_file(str(data / "demo.nwk"))
    alignment = hyphy3.Alignment.load(str(data / "demo.fasta"))
    result = hyphy3.diffubar(tree, alignment, settings, (data / "fixed-model.json").read_text())
    assert result["analysis"]["name"] == "difFUBAR"
    assert result["input"]["groups"] == ["G1", "G2"]
    assert len(result["sites"]) == 12
    assert result["grid"]["category_count"] == 128
    replay = hyphy3.diffubar(tree, alignment, settings, json.dumps(result))
    assert replay["sites"] == result["sites"]
    assert replay["grid"]["mean_weights"] == result["grid"]["mean_weights"]
    print("difFUBAR Python binding tests passed")

if __name__ == "__main__":
    main()
