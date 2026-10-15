#include "hyphy/analyses/diffubar.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace hyphy::analyses;
using namespace hyphy::core;

namespace {
void require(bool value, const std::string& message) {
    if (!value) throw std::runtime_error(message);
}
void close(double a, double b, double tolerance, const std::string& message) {
    if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a - b) > tolerance)
        throw std::runtime_error(message + ": " + std::to_string(a) + " vs " + std::to_string(b));
}
template<class F> void rejects(F f, const std::string& message) {
    bool rejected = false;
    try { f(); } catch (const std::exception&) { rejected = true; }
    require(rejected, message);
}
void close_vectors(const std::vector<double>& a, const std::vector<double>& b, double tolerance, const std::string& message) {
    require(a.size() == b.size(), message + " size");
    for (std::size_t i = 0; i < a.size(); ++i) close(a[i], b[i], tolerance, message);
}
DifFUBARSettings small_settings() {
    DifFUBARSettings s;
    s.foreground_grid = 2; s.background_grid = 1;
    s.sampler.iterations = 250; s.sampler.chains = 2; s.sampler.seed = 572;
    s.threads = 1; s.store_grid_likelihoods = true;
    return s;
}
DifFUBARFittedModel fixed_model() {
    DifFUBARFittedModel m;
    m.exchangeabilities = {1.2, 0.9, 0.7, 1.1, 2.1, 0.8};
    m.f3x4 << .30, .20, .20, .30,
               .20, .30, .30, .20,
               .25, .25, .25, .25;
    m.alpha = .8; m.beta = .3; m.genetic_code = "Universal";
    return m;
}
int node_named(const Tree& tree, const std::string& name) {
    for (const auto& node : tree.nodes) if (node.name == name) return node.id;
    throw std::runtime_error("Missing test node " + name);
}
}

int main(int argc, char** argv) {
    try {
        const std::string data = std::string(HYPHY3_SOURCE_DIR) + "/tests/data/diffubar/";
        const auto alignment = Alignment::load(data + "demo.fasta");
        const auto tree = Tree::from_newick_file(data + "demo.nwk");
        const auto config = small_settings();
        DifFUBARAnalyzer analyzer(tree, alignment, config);
        const auto raw_f3x4 = diffubar_detail::count_f3x4(alignment);
        require((raw_f3x4 - alignment.pos_nuc_frequencies).norm() > 1e-5,
                "Fixture must distinguish F3x4 from corrected CF3x4");
        require((analyzer.fitted_model.f3x4 - raw_f3x4).norm() == 0.0, "Use reference raw F3x4");
        analyzer.use_fitted_model(fixed_model());
        require(analyzer.groups.names == std::array<std::string, 2>{"G1", "G2"}, "Lexical group order");
        require(analyzer.groups.has_background, "Background present");
        require(analyzer.groups.branch_counts == std::array<std::size_t, 3>{3, 3, 3}, "Group edge counts");
        const auto grid = create_diffubar_grid(true, 2, 1);
        const auto logs = analyzer.compute_grid_log_likelihoods(grid);
        require(logs.size() == alignment.patterns.size() * grid.categories.size(), "Grid dimensions");
        for (double x : logs) require(std::isfinite(x), "Finite native likelihood");
        auto parallel = analyzer;
        parallel.settings.threads = 2;
        close_vectors(logs, parallel.compute_grid_log_likelihoods(grid), 0.0, "Thread reproducibility");

        auto scaled_model = fixed_model();
        for (double& r : scaled_model.exchangeabilities) r *= 3.0;
        scaled_model.alpha /= 3.0; scaled_model.beta /= 3.0;
        auto scaled = analyzer; scaled.use_fitted_model(scaled_model);
        close_vectors(logs, scaled.compute_grid_log_likelihoods(grid), 1e-10, "Absolute GTR/reference-rate scaling");

        auto reversed_config = config; reversed_config.groups = {"G2", "G1"};
        DifFUBARAnalyzer reversed(tree, alignment, reversed_config);
        reversed.use_fitted_model(fixed_model());
        const auto swapped = reversed.compute_grid_log_likelihoods(grid);
        const std::size_t W = grid.omega.size(), B = grid.background_omega.size(), G = grid.categories.size();
        for (std::size_t g = 0; g < G; ++g) {
            const auto& c = grid.categories[g];
            const auto r = ((c.alpha_index * W + c.omega_index[1]) * W + c.omega_index[0]) * B + c.omega_index[2];
            for (std::size_t p = 0; p < alignment.patterns.size(); ++p)
                close(logs[p * G + g], swapped[p * G + r], 1e-10, "Group exchange symmetry");
        }
        auto root_tagged = tree; root_tagged.nodes[root_tagged.root_id].model_tag = "IgnoredRootTag";
        DifFUBARAnalyzer root_analyzer(root_tagged, alignment, config);
        require(root_analyzer.groups.names == analyzer.groups.names, "Ignore root tag");
        auto no_background = tree;
        for (auto& n : no_background.nodes) if (n.id != no_background.root_id && n.model_tag.empty()) n.model_tag = "G1";
        DifFUBARAnalyzer no_bg(no_background, alignment, config);
        require(!no_bg.groups.has_background, "Untagged root must not create background dimension");
        no_bg.use_fitted_model(fixed_model());
        require(no_bg.compute_grid_log_likelihoods(create_diffubar_grid(false, 2, 1)).size() == alignment.patterns.size() * 64,
                "Three-dimensional no-background likelihood grid");

        const auto result = analyzer.run();
        require(result.site_results.size() == alignment.num_codons, "Site expansion");
        require(result.posterior.retained_draws == 400, "Sampler iterations/burn-in/chains");
        close(std::accumulate(result.posterior.mean_grid_weights.begin(), result.posterior.mean_grid_weights.end(), 0.0),
              1.0, 1e-12, "Grid simplex");
        for (std::size_t i = 0; i < result.site_results.size(); ++i) {
            const auto& site = result.site_results[i];
            require(site.site == i + 1, "One-based codon sites");
            close(site.prob_omega1_gt_omega2 + site.prob_omega2_gt_omega1 + site.prob_omega1_eq_omega2,
                  1.0, 1e-12, "Directional posterior simplex");
        }

        // Empirical zeros: the positive-pi support is a single absorbing codon.
        Alignment constant;
        constant.add_sequence("a", "ATGATG"); constant.add_sequence("b", "ATGATG"); constant.finalize();
        DifFUBARAnalyzer zero_frequency(Tree::from_newick("(a{A}:0.1,b{B}:0.2);"), constant, config);
        auto zero_model = fixed_model(); zero_model.f3x4 = diffubar_detail::count_f3x4(constant);
        zero_frequency.use_fitted_model(zero_model);
        for (double x : zero_frequency.compute_grid_log_likelihoods(create_diffubar_grid(false, 2, 1)))
            close(x, 0.0, 1e-12, "Exact zero-frequency support, no pseudocounts");

        // Native Wheatstone bridge: weights (1,0,0) recover the original topology.
        const auto binary = Tree::from_newick(
            "(((alpha{G1}:0.05,beta{G1}:0.05)V{G1}:0.04,gamma{G2}:0.05)U{G2}:0.02,"
            "delta{G2}:0.05,(epsilon:0.06,zeta:0.06):0.03);");
        auto network = binary;
        const int v = node_named(binary, "V"), u = node_named(binary, "U");
        WheatstoneBridge bridge;
        bridge.hub_node_id = u; bridge.bridge_node_id = v;
        bridge.children = {node_named(binary, "alpha"), node_named(binary, "beta"), node_named(binary, "gamma")};
        bridge.weights = {1, 0, 0}; bridge.branch_lengths = {.04, .07, .09};
        network.nodes[u].children.assign(bridge.children.begin(), bridge.children.end());
        for (int child : bridge.children) network.nodes[child].parent_id = u;
        network.bridges.push_back(bridge); network.compute_post_order(); network.rebuild_bridge_index();
        DifFUBARAnalyzer original(binary, alignment, config), bridged(network, alignment, config);
        original.use_fitted_model(fixed_model()); bridged.use_fitted_model(fixed_model());
        const auto binary_logs = original.compute_grid_log_likelihoods(grid);
        const auto bridge_logs = bridged.compute_grid_log_likelihoods(grid);
        close_vectors(binary_logs, bridge_logs, 1e-10, "Wheatstone resolution-one equivalence");
        require(bridged.groups.node_group[v] == 0 && bridged.groups.node_group[u] == 1,
                "Bridge uses original internal branch group, not hub group");
        std::array<std::vector<double>, 3> resolution_logs;
        for (std::size_t r = 0; r < 3; ++r) {
            auto single = bridged; single.tree.bridges[0].weights = {0, 0, 0}; single.tree.bridges[0].weights[r] = 1;
            resolution_logs[r] = single.compute_grid_log_likelihoods(grid);
        }
        auto mixture = bridged; mixture.tree.bridges[0].weights = {.2, .3, .5};
        const auto mixed_logs = mixture.compute_grid_log_likelihoods(grid);
        for (std::size_t i = 0; i < mixed_logs.size(); ++i) {
            const double shift = std::max({resolution_logs[0][i], resolution_logs[1][i], resolution_logs[2][i]});
            const double expected = shift + std::log(.2 * std::exp(resolution_logs[0][i] - shift) +
                .3 * std::exp(resolution_logs[1][i] - shift) + .5 * std::exp(resolution_logs[2][i] - shift));
            close(mixed_logs[i], expected, 1e-10, "Wheatstone likelihood mixture");
        }

        rejects([&] { auto t = tree; t.nodes[node_named(t, "alpha")].name = "unknown"; DifFUBARAnalyzer a(t, alignment, config); }, "Missing taxon");
        rejects([&] { auto s = config; s.groups = {"G1", "Absent"}; DifFUBARAnalyzer a(tree, alignment, s); }, "Missing group");
        rejects([&] { auto t = tree; for (auto& n : t.nodes) n.branch_length = 0.0; DifFUBARAnalyzer a(t, alignment, config); }, "Zero/missing branch lengths");
        rejects([&] { auto bad = grid; bad.alpha.clear(); analyzer.compute_grid_log_likelihoods(bad); }, "Empty grid axis");
        rejects([&] { auto bad = grid; bad.categories[0].omega_index[0] = 999; analyzer.compute_grid_log_likelihoods(bad); }, "Malformed grid index");
        rejects([&] { auto m = fixed_model(); m.f3x4(0, 0) = -1; auto a = analyzer; a.use_fitted_model(m); }, "Negative frequencies");
        rejects([&] { auto a = analyzer; a.settings.max_memory_mb = 1; a.compute_grid_log_likelihoods(create_diffubar_grid(true)); }, "Memory guard");

        // Reuse the repository's optimizers; fit smoke test, not a Julia-parity claim.
        auto fit_config = config; fit_config.fit_max_iterations = 300; fit_config.fit_tolerance = 1e-4;
        DifFUBARAnalyzer fitted(tree, alignment, fit_config);
        fitted.fit_global_model();
        require(std::isfinite(fitted.fitted_model.codon_log_likelihood), "Finite global codon fit");
        require((fitted.fitted_model.f3x4 - raw_f3x4).norm() == 0.0, "Fitting retains raw F3x4");
        require(fitted.tree.to_newick() == tree.to_newick(), "Fitting preserves input branch lengths and annotations");
        require(fitted.run().site_results.size() == alignment.num_codons, "Fitted-model end-to-end smoke test");

        // Optional raw export for independent scipy.linalg.expm validation.
        if (argc == 2) {
            std::ofstream out(argv[1]); require(bool(out), "Open validation output");
            out << std::setprecision(17);
            for (std::size_t g = 0; g < G; ++g) {
                const auto& c = grid.categories[g];
                out << c.alpha << '\t' << c.omega[0] << '\t' << c.omega[1] << '\t' << c.omega[2];
                for (auto p : alignment.site_to_pattern) out << '\t' << logs[p * G + g];
                out << '\n';
            }
            require(bool(out), "Write validation output");
        }
        // Independent reference values are generated by generate_diffubar_reference.py.
        std::ifstream reference(data + "scipy-log-likelihoods.tsv");
        require(bool(reference), "Open independent SciPy fixture");
        double max_error = 0.0;
        for (std::size_t g = 0; g < G; ++g)
            for (std::size_t s = 0; s < alignment.num_codons; ++s) {
                double expected = 0.0; require(bool(reference >> expected), "Read SciPy fixture");
                const double actual = logs[alignment.site_to_pattern[s] * G + g];
                max_error = std::max(max_error, std::abs(expected - actual));
                close(expected, actual, 1e-8, "Independent dense matrix-exponential parity");
            }
        std::string extra; require(!(reference >> extra), "Unexpected trailing reference values");
        std::cout << "difFUBAR native tests passed; max SciPy log-likelihood error " << std::setprecision(8) << max_error
                  << "; " << G << " categories x " << alignment.num_codons << " codon sites\n";
        std::cout << "Global fit smoke test lnL: " << fitted.fitted_model.codon_log_likelihood << '\n';
    } catch (const std::exception& e) {
        std::cerr << "difFUBAR native test failed: " << e.what() << '\n'; return 1;
    }
    return 0;
}
