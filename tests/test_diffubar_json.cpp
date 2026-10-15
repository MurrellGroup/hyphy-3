#include "hyphy/analyses/diffubar_json.hpp"
#include <iostream>

using namespace hyphy::analyses;
using namespace hyphy::core;

int main() {
    try {
        const std::string data = std::string(HYPHY3_SOURCE_DIR) + "/tests/data/diffubar/";
        auto require = [](bool condition, const std::string& message) {
            if (!condition) throw std::runtime_error(message);
        };
        std::ifstream input(data + "fixed-model.json");
        auto fixed_json = nlohmann::json::parse(input);
        auto model = diffubar_fitted_model_from_json(fixed_json);
        const auto serialized = diffubar_fitted_model_to_json(model);
        require(serialized["codon_log_likelihood"].is_null(), "Unknown likelihood is null, not a fabricated value");
        const auto roundtrip = diffubar_fitted_model_from_json(nlohmann::json::parse(serialized.dump()));
        require(roundtrip.exchangeabilities == model.exchangeabilities, "Exchangeability round-trip");
        require((roundtrip.f3x4 - model.f3x4).norm() == 0.0, "Frequency round-trip");
        DifFUBARSettings s;
        s.foreground_grid = 2; s.background_grid = 1; s.threads = 1;
        s.sampler.iterations = 80; s.sampler.seed = 12;
        DifFUBARAnalyzer a(Tree::from_newick_file(data + "demo.nwk"), Alignment::load(data + "demo.fasta"), s);
        a.use_fitted_model(model);
        const auto result = a.run();
        const auto json = nlohmann::json::parse(diffubar_to_json(a, result).dump());
        require(json.at("analysis").at("name") == "difFUBAR", "Analysis name");
        require(json.at("sites").size() == 12, "JSON site count");
        require(json.at("grid").at("category_count") == 128, "JSON category count");
        require(json.at("MLE").at("content").at("0").size() == 12, "MLE table shape");
        require(!json.at("diagnostics").at("convergence_assessed").get<bool>(), "Do not fabricate convergence diagnostics");
        require(json.at("input").at("wheatstone_bridges").empty(), "Ordinary tree metadata");
        require(diffubar_fitted_model_from_json(json).alpha == model.alpha, "Import model from prior result JSON");
        fixed_json["tree_branch_units"] = "already-alpha-rescaled";
        bool rejected = false;
        try { diffubar_fitted_model_from_json(fixed_json); } catch (const std::exception&) { rejected = true; }
        require(rejected, "Reject incompatible branch units");
        std::cout << "difFUBAR JSON tests passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
