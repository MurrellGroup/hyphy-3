#pragma once
#include "hyphy/analyses/diffubar.hpp"
#include "nlohmann/json.hpp"
#include <fstream>
#include <iomanip>

namespace hyphy::analyses {

inline nlohmann::json diffubar_fitted_model_to_json(const DifFUBARFittedModel& model) {
    nlohmann::json frequencies = nlohmann::json::array();
    for (int p = 0; p < 3; ++p)
        frequencies.push_back({model.f3x4(p, 0), model.f3x4(p, 1), model.f3x4(p, 2), model.f3x4(p, 3)});
    const auto finite_or_null = [](double x) -> nlohmann::json {
        return std::isfinite(x) ? nlohmann::json(x) : nlohmann::json(nullptr);
    };
    return {{"schema", "hyphy3-diffubar-fitted-model-v1"},
            {"exchangeability_order", {"AC", "AG", "AT", "CG", "CT", "GT"}},
            {"nucleotide_order", {"A", "C", "G", "T"}},
            {"exchangeabilities", model.exchangeabilities}, {"f3x4", frequencies},
            {"alpha", model.alpha}, {"beta", model.beta},
            {"tree_branch_units", "input-tree"}, {"genetic_code", model.genetic_code},
            {"fit_mode", model.fit_mode}, {"near_parameter_boundary", model.near_parameter_boundary},
            {"nucleotide_log_likelihood", finite_or_null(model.nucleotide_log_likelihood)},
            {"codon_log_likelihood", finite_or_null(model.codon_log_likelihood)}};
}

inline DifFUBARFittedModel diffubar_fitted_model_from_json(const nlohmann::json& document) {
    const auto& j = document.contains("fitted_model") ? document.at("fitted_model") : document;
    if (j.at("schema") != "hyphy3-diffubar-fitted-model-v1" || j.at("tree_branch_units") != "input-tree" ||
        j.at("exchangeability_order") != nlohmann::json({"AC", "AG", "AT", "CG", "CT", "GT"}) ||
        j.at("nucleotide_order") != nlohmann::json({"A", "C", "G", "T"}))
        throw std::invalid_argument("difFUBAR: incompatible fitted-model schema, nucleotide order, or branch units");
    if (!j.at("exchangeabilities").is_array() || j.at("exchangeabilities").size() != 6 ||
        !j.at("f3x4").is_array() || j.at("f3x4").size() != 3)
        throw std::invalid_argument("difFUBAR: fitted model requires 6 exchangeabilities and a 3x4 frequency matrix");
    DifFUBARFittedModel model;
    model.exchangeabilities = j.at("exchangeabilities").get<std::array<double, 6>>();
    for (int p = 0; p < 3; ++p) {
        if (!j.at("f3x4").at(p).is_array() || j.at("f3x4").at(p).size() != 4)
            throw std::invalid_argument("difFUBAR: fitted-model F3x4 must be 3x4");
        for (int n = 0; n < 4; ++n) model.f3x4(p, n) = j.at("f3x4").at(p).at(n).get<double>();
    }
    model.alpha = j.at("alpha").get<double>(); model.beta = j.at("beta").get<double>();
    model.genetic_code = j.at("genetic_code").get<std::string>();
    model.fit_mode = "supplied";
    if (j.contains("nucleotide_log_likelihood") && !j.at("nucleotide_log_likelihood").is_null())
        model.nucleotide_log_likelihood = j.at("nucleotide_log_likelihood").get<double>();
    if (j.contains("codon_log_likelihood") && !j.at("codon_log_likelihood").is_null())
        model.codon_log_likelihood = j.at("codon_log_likelihood").get<double>();
    model.near_parameter_boundary = j.value("near_parameter_boundary", false);
    return model;
}

inline nlohmann::json diffubar_to_json(const DifFUBARAnalyzer& analyzer, const DifFUBARResult& result) {
    const auto& settings = analyzer.settings;
    nlohmann::json j;
    j["analysis"] = {{"name", "difFUBAR"}, {"version", "0.1.0"},
        {"engine", "hyphy-3"}, {"model", "MG94-F3x4: shared alpha, two foreground omegas, optional background omega"},
        {"inference", "Uncollapsed Gibbs with Rao-Blackwellized site posteriors"}};
    j["input"] = {{"number_of_sequences", analyzer.aln.num_taxa}, {"number_of_sites", analyzer.aln.num_codons},
        {"number_of_patterns", analyzer.aln.patterns.size()}, {"tree", analyzer.tree.to_newick()},
        {"groups", result.groups.names}, {"has_background", result.groups.has_background},
        {"genetic_code", result.fitted_model.genetic_code}};
    j["settings"] = {{"foreground_grid", settings.foreground_grid}, {"background_grid", settings.background_grid},
        {"iterations", settings.sampler.iterations}, {"burn_in", settings.sampler.resolved_burn_in()},
        {"chains", settings.sampler.chains}, {"seed", settings.sampler.seed},
        {"concentration_per_category", settings.sampler.concentration},
        {"total_prior_concentration", settings.sampler.concentration * result.grid.categories.size()},
        {"posterior_threshold", settings.posterior_threshold}, {"retained_draws", result.posterior.retained_draws},
        {"branch_length_optimization", false}, {"likelihood_cutoff", 0.0},
        {"frequency_estimator", "raw-F3x4"}, {"threads", settings.threads},
        {"fit_max_iterations", settings.fit_max_iterations}, {"fit_tolerance", settings.fit_tolerance}};
    j["fitted_model"] = diffubar_fitted_model_to_json(result.fitted_model);
    j["timers"] = {{"global_fit", result.time_fit_sec}, {"grid", result.time_grid_sec},
                    {"posterior", result.time_posterior_sec}};
    j["grid"] = {{"alpha", result.grid.alpha}, {"omega", result.grid.omega},
                  {"background_omega", result.grid.has_background ? result.grid.background_omega : std::vector<double>{}},
                  {"category_count", result.grid.categories.size()},
                  {"mean_weights", result.posterior.mean_grid_weights}};
    j["grid"]["categories"] = nlohmann::json::array();
    for (const auto& c : result.grid.categories) {
        nlohmann::json row = {c.alpha, c.omega[0], c.omega[1]};
        if (result.grid.has_background) row.push_back(c.omega[2]);
        j["grid"]["categories"].push_back(std::move(row));
    }
    j["branch_attributes"] = nlohmann::json::array();
    for (std::size_t id = 0; id < analyzer.tree.nodes.size(); ++id) {
        const int group = result.groups.node_group[id];
        if (group < 0) continue;
        const auto& node = analyzer.tree.nodes[id];
        j["branch_attributes"].push_back({{"id", node.id}, {"name", node.name}, {"model_tag", node.model_tag},
            {"group", group == 2 ? "background" : result.groups.names[group]},
            {"input_length", node.branch_length},
            {"grid_baseline_length", node.branch_length * result.fitted_model.alpha * result.fitted_model.exchangeabilities[1]}});
    }
    j["input"]["tree_kind"] = analyzer.tree.bridges.empty() ? "tree" : "wheatstone";
    j["input"]["wheatstone_bridges"] = nlohmann::json::array();
    for (const auto& bridge : analyzer.tree.bridges)
        j["input"]["wheatstone_bridges"].push_back({{"hub_node_id", bridge.hub_node_id},
            {"bridge_node_id", bridge.bridge_node_id}, {"children", bridge.children},
            {"weights", bridge.weights}, {"input_branch_lengths", bridge.branch_lengths}});
    j["site_to_pattern"] = result.site_to_pattern;
    j["pattern_weights"] = result.pattern_weights;
    j["posterior_marginals"] = nlohmann::json::array();
    for (const auto& s : result.posterior.patterns)
        j["posterior_marginals"].push_back({{"alpha", s.alpha_marginal}, {"omega", s.omega_marginal}});
    j["chain_probabilities"] = result.posterior.chain_probabilities;
    j["diagnostics"] = {{"convergence_assessed", false},
        {"note", "Independent chain means are supplied; no ESS or R-hat is computed. Cross-language RNG streams differ."}};
    j["MLE"]["headers"] = {
        {"Codon Sites", "One-based codon position"}, {"P(omega1 > omega2)", "Directional posterior"},
        {"P(omega2 > omega1)", "Directional posterior"}, {"P(omega1 > 1)", "Positive selection in group 1"},
        {"P(omega2 > 1)", "Positive selection in group 2"}, {"mean(alpha)", "Mean synonymous rate"},
        {"mean(omega1)", "Mean group 1 dN/dS"}, {"mean(omega2)", "Mean group 2 dN/dS"}};
    j["MLE"]["content"]["0"] = nlohmann::json::array();
    j["sites"] = nlohmann::json::array();
    for (const auto& s : result.site_results) {
        j["MLE"]["content"]["0"].push_back({s.site, s.prob_omega1_gt_omega2, s.prob_omega2_gt_omega1,
            s.prob_positive[0], s.prob_positive[1], s.mean_alpha, s.mean_omega[0], s.mean_omega[1]});
        nlohmann::json site = {{"site", s.site}, {"pattern", s.pattern}, {"mean_alpha", s.mean_alpha},
            {"mean_omega", {s.mean_omega[0], s.mean_omega[1]}}, {"prob_omega1_gt_omega2", s.prob_omega1_gt_omega2},
            {"prob_omega2_gt_omega1", s.prob_omega2_gt_omega1}, {"prob_omega1_eq_omega2", s.prob_omega1_eq_omega2},
            {"prob_positive", s.prob_positive}, {"detected", s.detected}};
        if (result.grid.has_background) site["mean_background_omega"] = s.mean_omega[2];
        j["sites"].push_back(std::move(site));
    }
    return j;
}

// Category-major, expanded to original sites, matching the Julia harness.
inline void write_diffubar_grid_tsv(const DifFUBARResult& result, const std::string& path) {
    const std::size_t G = result.grid.categories.size();
    if (result.grid_log_likelihoods.size() != result.pattern_weights.size() * G)
        throw std::invalid_argument("difFUBAR: enable store_grid_likelihoods before exporting the likelihood grid");
    std::ofstream output(path);
    if (!output) throw std::runtime_error("difFUBAR: cannot open likelihood export: " + path);
    output << std::setprecision(17);
    for (std::size_t g = 0; g < G; ++g) {
        for (std::size_t s = 0; s < result.site_to_pattern.size(); ++s) {
            if (s) output << '\t';
            output << result.grid_log_likelihoods[result.site_to_pattern[s] * G + g];
        }
        output << '\n';
    }
    if (!output) throw std::runtime_error("difFUBAR: failed while writing likelihood export: " + path);
}

} // namespace hyphy::analyses
