#pragma once

#include "hyphy/analyses/detail/diffubar_inference.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/opt/nelder_mead.hpp"

#include <atomic>
#include <chrono>
#include <cctype>
#include <exception>
#include <map>
#include <set>
#include <unordered_set>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

struct DifFUBARSettings {
    DifFUBARSamplerSettings sampler;
    int foreground_grid = 6;
    int background_grid = 4;
    // Empty: infer exactly two distinct non-root tags in lexical order.
    // Explicit: preserve this order; all other branches become background.
    std::vector<std::string> groups;
    double posterior_threshold = 0.95;
    int threads = 0; // 0: runtime default; never changes global OpenMP settings
    std::size_t max_memory_mb = 1024;
    int fit_max_iterations = 2000;
    double fit_tolerance = 1e-6;
    bool store_grid_likelihoods = false;

    void validate() const {
        sampler.validate();
        if (foreground_grid < 1 || foreground_grid > 100 || background_grid < 1 || background_grid > 100)
            throw std::invalid_argument("difFUBAR: grid densities must lie in 1..100");
        if (!groups.empty() && (groups.size() != 2 || groups[0].empty() ||
                                groups[1].empty() || groups[0] == groups[1]))
            throw std::invalid_argument("difFUBAR: supply exactly two distinct, nonempty group tags");
        if (!(posterior_threshold > 0.0 && posterior_threshold <= 1.0))
            throw std::invalid_argument("difFUBAR: posterior threshold must lie in (0,1]");
        if (threads < 0 || max_memory_mb == 0 || fit_max_iterations <= 0 ||
            !(fit_tolerance > 0.0) || !std::isfinite(fit_tolerance))
            throw std::invalid_argument("difFUBAR: invalid thread, memory, or fitting setting");
        diffubar_detail::checked_product(max_memory_mb, 1024 * 1024, "memory budget");
    }
};

struct DifFUBARFittedModel {
    // Absolute exchangeabilities in A,C,G,T order: AC, AG, AT, CG, CT, GT.
    // The grid uses exp[Q(r/r_AG; 1,omega) * alpha_grid *
    //               (input_branch_length * alpha * r_AG)].
    // This explicitly preserves Julia's raw-GTR and alpha-rescaled-tree units.
    std::array<double, 6> exchangeabilities{1, 1, 1, 1, 1, 1};
    Eigen::Matrix<core::Scalar, 3, 4> f3x4 = Eigen::Matrix<core::Scalar, 3, 4>::Constant(0.25);
    double alpha = 1.0;
    double beta = 1.0;
    std::string genetic_code;
    std::string fit_mode = "supplied";
    double nucleotide_log_likelihood = std::numeric_limits<double>::quiet_NaN();
    double codon_log_likelihood = std::numeric_limits<double>::quiet_NaN();
    bool near_parameter_boundary = false;
};

struct DifFUBARBranchGroups {
    std::array<std::string, 2> names;
    std::vector<int> node_group; // 0, 1, 2 (background); root/unused = -1
    std::array<std::size_t, 3> branch_counts{0, 0, 0};
    bool has_background = false;
};

struct DifFUBARSiteResult {
    std::size_t site = 0; // one-based codon site
    std::size_t pattern = 0;
    double mean_alpha = 0.0;
    std::array<double, 3> mean_omega{0, 0, 0};
    double prob_omega1_gt_omega2 = 0.0;
    double prob_omega2_gt_omega1 = 0.0;
    double prob_omega1_eq_omega2 = 0.0;
    std::array<double, 2> prob_positive{0, 0};
    bool detected = false;
};

struct DifFUBARResult {
    DifFUBARGrid grid;
    DifFUBARFittedModel fitted_model;
    DifFUBARBranchGroups groups;
    DifFUBARPosterior posterior;
    std::vector<DifFUBARSiteResult> site_results;
    std::vector<std::size_t> site_to_pattern;
    std::vector<std::size_t> pattern_weights;
    // Optional, pattern-major [p * G + g]; not normalized or rounded.
    std::vector<double> grid_log_likelihoods;
    double time_fit_sec = 0.0;
    double time_grid_sec = 0.0;
    double time_posterior_sec = 0.0;
};

namespace diffubar_detail {

inline core::MG94Parameters relative_rates(const DifFUBARFittedModel& model) {
    core::MG94Parameters p;
    const double reference = model.exchangeabilities[1];
    p.theta_AC = model.exchangeabilities[0] / reference;
    p.theta_AT = model.exchangeabilities[2] / reference;
    p.theta_CG = model.exchangeabilities[3] / reference;
    p.theta_CT = model.exchangeabilities[4] / reference;
    p.theta_GT = model.exchangeabilities[5] / reference;
    p.alpha = p.beta = 1.0;
    return p;
}

inline core::Vector equilibrium(const DifFUBARFittedModel& model, const core::GeneticCode& code) {
    core::Vector pi(code.num_sense_codons);
    for (int c = 0; c < code.num_sense_codons; ++c) {
        const auto& triplet = code.sense_codons[c];
        pi(c) = model.f3x4(0, core::char_to_nuc(triplet[0])) *
                model.f3x4(1, core::char_to_nuc(triplet[1])) *
                model.f3x4(2, core::char_to_nuc(triplet[2]));
    }
    if (!(pi.sum() > 0.0)) throw std::invalid_argument("difFUBAR: no sense-codon equilibrium mass");
    pi /= pi.sum();
    return pi;
}

// Raw F3x4, not Alignment's corrected CF3x4. The denominator includes
// missing bases, matching the reference count_F3x4 convention.
inline Eigen::Matrix<core::Scalar, 3, 4> count_f3x4(const core::Alignment& aln) {
    Eigen::Matrix<core::Scalar, 3, 4> frequencies = Eigen::Matrix<core::Scalar, 3, 4>::Zero();
    const double denominator = static_cast<double>(aln.num_taxa) * static_cast<double>(aln.num_codons);
    for (const auto& sequence : aln.sequences)
        for (std::size_t i = 0; i < sequence.raw_seq.size(); ++i) {
            const int base = core::char_to_nuc(sequence.raw_seq[i]);
            if (base >= 0) frequencies(static_cast<int>(i % 3), base) += 1.0;
        }
    return frequencies / denominator;
}

// Use the native MG94 generator and spectral implementation. At boundary
// frequencies, restrict the spectral problem to the closed positive-pi support
// rather than applying the core's sqrt(pi) floor. Excluded states have zero root
// mass and cannot be entered from that support, so their transition rows do not
// affect the likelihood. This handles empirical zeros without pseudocounts.
class MG94Model {
    core::MG94Matrix native_;
    std::vector<int> support_;
    core::Vector eigenvalues_;
    core::Matrix left_, right_;
    bool restricted_ = false;
public:
    void update(const core::MG94Parameters& parameters,
                const Eigen::Matrix<core::Scalar, 3, 4>& f3x4,
                const core::Vector& pi, const core::GeneticCode& code) {
        native_.update(parameters, f3x4, pi, code);
        restricted_ = pi.minCoeff() < 1e-12;
        if (!restricted_) return;
        support_.clear();
        for (int i = 0; i < pi.size(); ++i) if (pi(i) > 0.0) support_.push_back(i);
        const int n = static_cast<int>(support_.size());
        core::Vector sqrt_pi(n);
        core::Matrix symmetric(n, n);
        for (int i = 0; i < n; ++i) sqrt_pi(i) = std::sqrt(pi(support_[i]));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                symmetric(i, j) = sqrt_pi(i) * native_.Q(support_[i], support_[j]) / sqrt_pi(j);
        symmetric = (0.5 * (symmetric + symmetric.transpose())).eval();
        Eigen::SelfAdjointEigenSolver<core::Matrix> solver(symmetric);
        if (solver.info() != Eigen::Success) throw std::runtime_error("difFUBAR: support eigendecomposition failed");
        eigenvalues_ = solver.eigenvalues().cwiseMin(0.0);
        left_ = sqrt_pi.cwiseInverse().asDiagonal() * solver.eigenvectors();
        right_ = solver.eigenvectors().transpose() * sqrt_pi.asDiagonal();
    }
    core::Matrix transition_matrix(double length) const {
        if (!restricted_) return native_.transition_matrix(length);
        core::Matrix result = core::Matrix::Identity(native_.Q.rows(), native_.Q.cols());
        const core::Matrix block = left_ * (eigenvalues_ * length).array().exp().matrix().asDiagonal() * right_;
        for (std::size_t i = 0; i < support_.size(); ++i)
            for (std::size_t j = 0; j < support_.size(); ++j) result(support_[i], support_[j]) = block(i, j);
        return result;
    }
};

struct PruningWorkspace {
    std::vector<core::Vector> likelihood;
    std::vector<double> log_scale;
    PruningWorkspace(std::size_t nodes, int states)
        : likelihood(nodes, core::Vector::Zero(states)), log_scale(nodes, 0.0) {}
};

// Native Tree handles Wheatstone hubs. Ordinary polytomies are normalized after
// every child combine (rather than after an arbitrarily large product).
template<class MatrixType>
inline double pattern_log_likelihood(
    const core::Tree& tree, const core::SitePattern& pattern,
    const std::vector<std::size_t>& leaf_to_taxon, const core::Vector& pi,
    const std::vector<MatrixType>& transitions,
    const std::vector<std::array<MatrixType, 3>>& bridge_transitions,
    PruningWorkspace& work) {
    const int states = static_cast<int>(pi.size());
    for (const std::int32_t id : tree.post_order) {
        const auto& node = tree.nodes[id];
        auto& value = work.likelihood[id];
        if (node.is_leaf) {
            const int state = pattern.states.at(leaf_to_taxon[id]);
            if (state >= states) throw std::invalid_argument("difFUBAR: alignment contains an invalid state");
            if (state >= 0) { value.setZero(); value(state) = 1.0; }
            else value.setOnes(); // native hyphy-3 gap/ambiguous/stop semantics
            work.log_scale[id] = 0.0;
        } else {
            double scale = 0.0;
            if (tree.bridge_index(id) >= 0) {
                tree.eval_internal_node(id, transitions, work.likelihood, work.log_scale,
                                        value, scale, &bridge_transitions);
                const double sum = value.sum();
                if (!(sum > 0.0)) return -std::numeric_limits<double>::infinity();
                value /= sum;
                scale += std::log(sum);
            } else {
                value.setOnes();
                for (const std::int32_t child : node.children) {
                    value.array() *= (transitions[child] * work.likelihood[child]).array();
                    scale += work.log_scale[child];
                    const double sum = value.sum();
                    if (!(sum > 0.0)) return -std::numeric_limits<double>::infinity();
                    value /= sum;
                    scale += std::log(sum);
                }
            }
            work.log_scale[id] = scale;
        }
    }
    const double root_likelihood = pi.dot(work.likelihood[tree.root_id]);
    return root_likelihood > 0.0 ? std::log(root_likelihood) + work.log_scale[tree.root_id]
                                : -std::numeric_limits<double>::infinity();
}

template<class Model>
inline core::Matrix transition(const Model& model, double length, int states) {
    if (!(length >= 0.0) || !std::isfinite(length))
        throw std::invalid_argument("difFUBAR: invalid scaled branch length");
    // Spectral roundoff at t=0 must not turn an impossible pattern into a
    // possible one. For t>0, only clip numerical negative transition entries.
    if (length == 0.0) return core::Matrix::Identity(states, states);
    core::Matrix P = model.transition_matrix(length);
    if (!P.allFinite() || P.minCoeff() < -1e-9)
        throw std::runtime_error("difFUBAR: invalid transition matrix");
    P = P.cwiseMax(0.0);
    for (int i = 0; i < states; ++i) {
        const double sum = P.row(i).sum();
        if (!(sum > 0.0)) throw std::runtime_error("difFUBAR: zero transition row");
        P.row(i) /= sum;
    }
    return P;
}

} // namespace diffubar_detail

class DifFUBARAnalyzer {
public:
    core::Tree tree;
    core::Alignment aln;
    DifFUBARSettings settings;
    DifFUBARFittedModel fitted_model;
    DifFUBARBranchGroups groups;

    DifFUBARAnalyzer(core::Tree input_tree, core::Alignment alignment, DifFUBARSettings config = {})
        : tree(std::move(input_tree)), aln(std::move(alignment)), settings(std::move(config)) {
        settings.validate();
        validate_data();
        prepare_groups();
        fitted_model.f3x4 = diffubar_detail::count_f3x4(aln);
        fitted_model.genetic_code = code().name;
        fitted_model.fit_mode = "unfitted";
    }

    static DifFUBARAnalyzer create_and_fit(core::Tree tree, core::Alignment aln,
                                           DifFUBARSettings settings = {},
                                           std::function<void(const std::string&, double)> progress = {}) {
        DifFUBARAnalyzer analyzer(std::move(tree), std::move(aln), std::move(settings));
        analyzer.fit_global_model(std::move(progress));
        return analyzer;
    }

    void use_fitted_model(DifFUBARFittedModel model) {
        validate_model(model);
        fitted_model = std::move(model);
        if (fitted_model.genetic_code.empty()) fitted_model.genetic_code = code().name;
        if (fitted_model.fit_mode == "unfitted") fitted_model.fit_mode = "supplied";
    }

    // Two-step fit on fixed input branch lengths, matching difFUBAR's default
    // parameterization. GTRMatrix normalizes Q; its scale_factor restores the
    // six absolute exchangeabilities for this fixed-tree nucleotide fit.
    void fit_global_model(std::function<void(const std::string&, double)> progress = {}) {
        const auto start = std::chrono::steady_clock::now();
        if (progress) progress("Fitting nucleotide exchangeabilities (fixed branches)", 0.0);
        const auto nucleotide_alignment = aln.to_nucleotide_alignment();
        if (!nucleotide_alignment.nuc_frequencies.allFinite() || nucleotide_alignment.nuc_frequencies.minCoeff() <= 0.0)
            throw std::invalid_argument("difFUBAR: global nucleotide fitting needs all four bases represented; supply a fitted model for this degenerate alignment");
        constexpr double bound = 9.210340371976184; // log(10000)
        const std::vector<double> lower6(6, -bound), upper6(6, bound);
        auto nucleotide_objective = [&](const std::vector<double>& x) {
            DifFUBARFittedModel candidate;
            for (std::size_t i = 0; i < 6; ++i) candidate.exchangeabilities[i] = std::exp(x[i]);
            const auto ratios = diffubar_detail::relative_rates(candidate);
            core::GTRParameters p;
            p.theta_AC = ratios.theta_AC; p.theta_AT = ratios.theta_AT;
            p.theta_CG = ratios.theta_CG; p.theta_CT = ratios.theta_CT; p.theta_GT = ratios.theta_GT;
            core::GTRMatrix model;
            model.update(p, nucleotide_alignment.nuc_frequencies);
            return -homogeneous_log_likelihood(nucleotide_alignment, model, core::Vector(model.pi),
                model.scale_factor * candidate.exchangeabilities[1]);
        };
        auto nuc_opt = opt::NelderMeadDynamic::minimize(nucleotide_objective,
            std::vector<double>(6, 0.0), lower6, upper6,
            settings.fit_tolerance, settings.fit_max_iterations);
        // Coordinate polishing uses the existing bounded optimizer, without
        // introducing a separate optimizer implementation or editing branches.
        for (int pass = 0; pass < 2; ++pass) {
            for (std::size_t k = 0; k < 6; ++k) {
                const auto objective = [&](double value) {
                    auto x = nuc_opt.x; x[k] = value; return nucleotide_objective(x);
                };
                const auto optimum = opt::Brent1D::minimize(objective, -bound, nuc_opt.x[k], bound,
                                                           settings.fit_tolerance, 100);
                if (optimum.second < nuc_opt.f) { nuc_opt.x[k] = optimum.first; nuc_opt.f = optimum.second; }
            }
        }
        DifFUBARFittedModel fit;
        fit.f3x4 = diffubar_detail::count_f3x4(aln);
        fit.genetic_code = code().name;
        fit.fit_mode = "two-step-fixed-branches";
        fit.nucleotide_log_likelihood = -nuc_opt.f;
        for (std::size_t i = 0; i < 6; ++i) {
            fit.exchangeabilities[i] = std::exp(nuc_opt.x[i]);
            fit.near_parameter_boundary |= std::abs(nuc_opt.x[i]) > bound - 0.01;
        }
        validate_model(fit);
        if (progress) progress("Fitting global synonymous/nonsynonymous rates", 0.5);
        const core::Vector pi = diffubar_detail::equilibrium(fit, code());
        auto codon_objective = [&](const std::vector<double>& x) {
            auto p = diffubar_detail::relative_rates(fit);
            p.alpha = std::exp(x[0]); p.beta = std::exp(x[1]);
            diffubar_detail::MG94Model model;
            model.update(p, fit.f3x4, pi, code());
            return -homogeneous_log_likelihood(aln, model, pi, fit.exchangeabilities[1]);
        };
        const auto codon_opt = opt::NelderMeadDynamic::minimize(codon_objective, {0.0, 0.0},
            {-bound, -bound}, {bound, bound}, settings.fit_tolerance, settings.fit_max_iterations,
            {{0.0, std::log(0.2)}, {0.0, std::log(2.0)}});
        fit.alpha = std::exp(codon_opt.x[0]); fit.beta = std::exp(codon_opt.x[1]);
        fit.codon_log_likelihood = -codon_opt.f;
        fit.near_parameter_boundary |= std::abs(codon_opt.x[0]) > bound - 0.01 ||
                                       std::abs(codon_opt.x[1]) > bound - 0.01;
        if (!std::isfinite(fit.nucleotide_log_likelihood) || !std::isfinite(fit.codon_log_likelihood))
            throw std::runtime_error("difFUBAR: non-finite global model fit; check the input tree/alignment");
        use_fitted_model(std::move(fit));
        time_fit_sec_ = seconds_since(start);
        if (progress) progress("Global model fit completed", 1.0);
    }

    // Exposed separately for deterministic fixed-model Julia/EvoOnline parity.
    std::vector<double> compute_grid_log_likelihoods(const DifFUBARGrid& grid,
        std::function<void(const std::string&, double)> progress = {}) const {
        settings.validate();
        validate_model(fitted_model);
        if (grid.has_background != groups.has_background || grid.categories.empty() ||
            grid.alpha.empty() || grid.omega.empty() || (grid.has_background && grid.background_omega.empty()))
            throw std::invalid_argument("difFUBAR: grid/background mismatch");
        const std::size_t N = tree.num_nodes(), P = aln.patterns.size(), G = grid.categories.size();
        std::size_t expected = diffubar_detail::checked_product(grid.alpha.size(), grid.omega.size(), "grid");
        expected = diffubar_detail::checked_product(expected, grid.omega.size(), "grid");
        expected = diffubar_detail::checked_product(expected, grid.has_background ? grid.background_omega.size() : 1, "grid");
        if (expected != G) throw std::invalid_argument("difFUBAR: non-Cartesian category grid");
        std::size_t category_index = 0;
        for (std::size_t a = 0; a < grid.alpha.size(); ++a)
            for (std::size_t i = 0; i < grid.omega.size(); ++i)
                for (std::size_t j = 0; j < grid.omega.size(); ++j)
                    for (std::size_t b = 0; b < (grid.has_background ? grid.background_omega.size() : 1); ++b) {
                        const auto& c = grid.categories[category_index++];
                        if (!(c.alpha > 0.0) || !std::isfinite(c.alpha) || c.alpha != grid.alpha[a] || c.alpha_index != a ||
                            c.omega_index != std::array<std::size_t, 3>{i, j, b} ||
                            c.omega != std::array<double, 3>{grid.omega[i], grid.omega[j], grid.has_background ? grid.background_omega[b] : 1.0})
                            throw std::invalid_argument("difFUBAR: invalid category values or ordering");
                        for (double w : c.omega)
                            if (!(w > 0.0) || !std::isfinite(w)) throw std::invalid_argument("difFUBAR: invalid omega grid");
                    }
        const int S = code().num_sense_codons;
        const std::size_t entries = diffubar_detail::checked_product(P, G, "likelihood grid");
        const std::size_t matrix_bytes = diffubar_detail::checked_product(entries, sizeof(double), "likelihood grid bytes");
        const std::size_t limit = settings.max_memory_mb * 1024 * 1024;
        std::size_t cached_matrices = 0;
        for (const auto id : tree.post_order)
            if (id != tree.root_id) cached_matrices += axis(grid, groups.node_group[id]).size();
        for (const auto& bridge : tree.bridges)
            cached_matrices += 3 * axis(grid, groups.node_group[bridge.bridge_node_id]).size();
        const std::size_t cache_bytes = diffubar_detail::checked_product(cached_matrices,
            static_cast<std::size_t>(S) * S * sizeof(double), "transition cache bytes");
        if (matrix_bytes > limit / 2 || cache_bytes > limit - matrix_bytes)
            throw std::length_error("difFUBAR: likelihood/transition storage exceeds --max-memory-mb; reduce grid density or raise the limit");
        std::vector<double> logs(entries, -std::numeric_limits<double>::infinity());
        const auto pi = diffubar_detail::equilibrium(fitted_model, code());
        std::array<std::vector<diffubar_detail::MG94Model>, 2> models;
        for (std::size_t kind = 0; kind < (groups.has_background ? 2u : 1u); ++kind) {
            for (const double omega : (kind == 0 ? grid.omega : grid.background_omega)) {
                auto rates = diffubar_detail::relative_rates(fitted_model);
                rates.beta = omega; // alpha=1; alpha_grid is a branch multiplier
                diffubar_detail::MG94Model model;
                model.update(rates, fitted_model.f3x4, pi, code());
                models[kind].push_back(std::move(model));
            }
        }
        const double baseline_scale = fitted_model.alpha * fitted_model.exchangeabilities[1];
        const core::Matrix identity = core::Matrix::Identity(S, S);
        using MatrixRef = Eigen::Ref<const core::Matrix>;
        const std::size_t block = G / grid.alpha.size();
        if (block * grid.alpha.size() != G) throw std::invalid_argument("difFUBAR: malformed category grid");
        for (std::size_t a = 0; a < grid.alpha.size(); ++a) {
            const double scale = grid.alpha[a] * baseline_scale;
            std::vector<std::vector<core::Matrix>> cache(N);
            std::vector<std::vector<std::array<core::Matrix, 3>>> bridge_cache(tree.bridges.size());
            for (const auto id : tree.post_order) {
                if (id == tree.root_id) continue;
                const int kind = groups.node_group[id] == 2 ? 1 : 0;
                for (const auto& model : models[kind])
                    cache[id].push_back(diffubar_detail::transition(model, tree.nodes[id].branch_length * scale, S));
            }
            for (std::size_t b = 0; b < tree.bridges.size(); ++b) {
                const auto& bridge = tree.bridges[b];
                const int kind = groups.node_group[bridge.bridge_node_id] == 2 ? 1 : 0;
                for (const auto& model : models[kind]) {
                    std::array<core::Matrix, 3> matrices;
                    for (std::size_t r = 0; r < 3; ++r)
                        matrices[r] = diffubar_detail::transition(model, bridge.branch_lengths[r] * scale, S);
                    bridge_cache[b].push_back(std::move(matrices));
                }
            }
            std::exception_ptr failure;
            std::atomic<bool> failed{false};
#ifdef _OPENMP
            const int threads = settings.threads > 0 ? settings.threads : omp_get_max_threads();
#pragma omp parallel num_threads(threads)
#endif
            {
                diffubar_detail::PruningWorkspace work(N, S);
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
                for (std::size_t offset = 0; offset < block; ++offset) {
                    if (failed.load(std::memory_order_relaxed)) continue;
                    try {
                        const std::size_t g = a * block + offset;
                        const auto& category = grid.categories[g];
                        if (category.alpha_index != a) throw std::invalid_argument("difFUBAR: unexpected grid ordering");
                        std::vector<MatrixRef> transitions;
                        transitions.reserve(N);
                        for (std::size_t id = 0; id < N; ++id) {
                            if (cache[id].empty()) transitions.emplace_back(identity);
                            else transitions.emplace_back(cache[id].at(category.omega_index[groups.node_group[id]]));
                        }
                        std::vector<std::array<MatrixRef, 3>> bridge_transitions;
                        bridge_transitions.reserve(tree.bridges.size());
                        for (std::size_t b = 0; b < tree.bridges.size(); ++b) {
                            const int group = groups.node_group[tree.bridges[b].bridge_node_id];
                            const auto& matrices = bridge_cache[b].at(category.omega_index[group]);
                            bridge_transitions.push_back({MatrixRef(matrices[0]), MatrixRef(matrices[1]), MatrixRef(matrices[2])});
                        }
                        for (std::size_t p = 0; p < P; ++p)
                            logs[p * G + g] = diffubar_detail::pattern_log_likelihood(tree, aln.patterns[p],
                                leaf_to_taxon_, pi, transitions, bridge_transitions, work);
                    } catch (...) {
#ifdef _OPENMP
#pragma omp critical(hyphy_diffubar_failure)
#endif
                        { if (!failure) failure = std::current_exception(); }
                        failed.store(true, std::memory_order_relaxed);
                    }
                }
            }
            if (failure) std::rethrow_exception(failure);
            if (progress) progress("Computing difFUBAR conditional likelihoods", static_cast<double>(a + 1) / grid.alpha.size());
        }
        return logs;
    }

    DifFUBARResult run(std::function<void(const std::string&, double)> progress = {}) const {
        if (fitted_model.fit_mode == "unfitted")
            throw std::logic_error("difFUBAR: call create_and_fit, fit_global_model, or use_fitted_model before run");
        DifFUBARResult result;
        result.grid = create_diffubar_grid(groups.has_background, settings.foreground_grid, settings.background_grid);
        result.groups = groups; result.fitted_model = fitted_model; result.time_fit_sec = time_fit_sec_;
        result.site_to_pattern = aln.site_to_pattern;
        for (const auto& pattern : aln.patterns) result.pattern_weights.push_back(pattern.weight);
        auto start = std::chrono::steady_clock::now();
        auto logs = compute_grid_log_likelihoods(result.grid, [&](const std::string& message, double fraction) {
            if (progress) progress(message, 0.7 * fraction);
        });
        result.time_grid_sec = seconds_since(start);
        start = std::chrono::steady_clock::now();
        result.posterior = sample_diffubar(result.grid, logs, result.pattern_weights, settings.sampler,
            [&](std::size_t completed, std::size_t total) {
                if (progress) progress("Sampling difFUBAR posterior", 0.7 + 0.3 * static_cast<double>(completed) / total);
            });
        result.time_posterior_sec = seconds_since(start);
        for (std::size_t site = 0; site < aln.num_codons; ++site) {
            const std::size_t pattern = aln.site_to_pattern[site];
            const auto& s = result.posterior.patterns[pattern];
            const bool detected = std::max({s.prob_omega1_gt_omega2, s.prob_omega2_gt_omega1,
                                            s.prob_positive[0], s.prob_positive[1]}) >= settings.posterior_threshold;
            result.site_results.push_back({site + 1, pattern, s.mean_alpha, s.mean_omega,
                s.prob_omega1_gt_omega2, s.prob_omega2_gt_omega1, s.prob_omega1_eq_omega2, s.prob_positive, detected});
        }
        if (settings.store_grid_likelihoods) result.grid_log_likelihoods = std::move(logs);
        return result;
    }

private:
    std::vector<std::size_t> leaf_to_taxon_;
    std::vector<std::int32_t> modeled_edges_;
    double time_fit_sec_ = 0.0;
    const core::GeneticCode& code() const { return *(aln.code ? aln.code : core::GeneticCode::universal()); }
    static double seconds_since(std::chrono::steady_clock::time_point start) {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }
    static const std::vector<double>& axis(const DifFUBARGrid& grid, int group) {
        return group == 2 ? grid.background_omega : grid.omega;
    }
    static std::string lower(std::string text) {
        std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) { return std::tolower(c); });
        return text;
    }
    void validate_model(const DifFUBARFittedModel& model) const {
        if (!model.genetic_code.empty() && model.genetic_code != code().name)
            throw std::invalid_argument("difFUBAR: fitted model and alignment use different genetic codes");
        for (double r : model.exchangeabilities)
            if (!(r > 0.0) || !std::isfinite(r)) throw std::invalid_argument("difFUBAR: exchangeabilities must be finite and positive");
        if (!(model.alpha > 0.0 && model.beta > 0.0) || !std::isfinite(model.alpha) || !std::isfinite(model.beta))
            throw std::invalid_argument("difFUBAR: global alpha/beta must be finite and positive");
        for (int pos = 0; pos < 3; ++pos) {
            if (!model.f3x4.row(pos).allFinite() || model.f3x4.row(pos).minCoeff() < 0.0 ||
                !(model.f3x4.row(pos).sum() > 0.0 && model.f3x4.row(pos).sum() <= 1.0 + 1e-8))
                throw std::invalid_argument("difFUBAR: each F3x4 row must be finite, nonnegative, and have total mass in (0,1]");
        }
    }
    void validate_data() {
        const std::size_t N = tree.num_nodes();
        if (N == 0 || tree.root_id < 0 || static_cast<std::size_t>(tree.root_id) >= N ||
            tree.post_order.empty() || tree.post_order.back() != tree.root_id)
            throw std::invalid_argument("difFUBAR: empty or invalid phylogenetic tree");
        if (aln.num_codons == 0 || aln.num_taxa < 2 || aln.patterns.empty() ||
            aln.site_to_pattern.size() != aln.num_codons || aln.num_nucleotides != 3 * aln.num_codons)
            throw std::invalid_argument("difFUBAR: a nonempty, in-frame codon alignment with at least two taxa is required");
        if (aln.taxon_names.size() != aln.num_taxa || aln.sequences.size() != aln.num_taxa)
            throw std::invalid_argument("difFUBAR: inconsistent alignment taxon names");
        for (const auto& sequence : aln.sequences)
            if (sequence.raw_seq.size() != aln.num_nucleotides)
                throw std::invalid_argument("difFUBAR: inconsistent raw sequence lengths");
        for (std::size_t id = 0; id < N; ++id)
            if (tree.nodes[id].id != static_cast<std::int32_t>(id))
                throw std::invalid_argument("difFUBAR: tree node IDs must index tree.nodes");
        std::vector<bool> visited(N, false), used_taxa(aln.num_taxa, false);
        leaf_to_taxon_.assign(N, std::numeric_limits<std::size_t>::max());
        std::set<std::int32_t> edges;
        for (const auto id : tree.post_order) {
            if (id < 0 || static_cast<std::size_t>(id) >= N || visited[id])
                throw std::invalid_argument("difFUBAR: invalid post-order traversal");
            const auto& node = tree.nodes[id];
            for (const auto child : node.children)
                if (child < 0 || static_cast<std::size_t>(child) >= N || !visited[child])
                    throw std::invalid_argument("difFUBAR: traversal is not in post-order");
            visited[id] = true;
            if (id != tree.root_id) edges.insert(id);
            if (node.is_leaf) {
                std::size_t found = aln.num_taxa;
                const auto it = aln.taxon_to_index.find(node.name);
                if (it != aln.taxon_to_index.end()) found = it->second;
                else {
                    for (std::size_t taxon = 0; taxon < aln.num_taxa; ++taxon)
                        if (lower(node.name) == lower(aln.taxon_names[taxon])) {
                            if (found != aln.num_taxa) throw std::invalid_argument("difFUBAR: ambiguous case-insensitive taxon match");
                            found = taxon;
                        }
                }
                if (found >= aln.num_taxa || used_taxa[found])
                    throw std::invalid_argument("difFUBAR: missing or duplicate alignment taxon for tree tip '" + node.name + "'");
                leaf_to_taxon_[id] = found; used_taxa[found] = true;
            }
        }
        if (std::find(used_taxa.begin(), used_taxa.end(), false) != used_taxa.end())
            throw std::invalid_argument("difFUBAR: alignment and tree must contain the same taxa");
        tree.rebuild_bridge_index();
        for (const auto& bridge : tree.bridges) {
            if (bridge.bridge_node_id < 0 || static_cast<std::size_t>(bridge.bridge_node_id) >= N ||
                bridge.bridge_node_id == tree.root_id || bridge.bridge_node_id == bridge.hub_node_id || bridge.hub_node_id < 0 ||
                static_cast<std::size_t>(bridge.hub_node_id) >= N || !visited[bridge.hub_node_id])
                throw std::invalid_argument("difFUBAR: Wheatstone bridge needs a valid original non-root branch for its group");
            edges.insert(bridge.bridge_node_id);
            double weight_sum = 0.0;
            for (std::size_t r = 0; r < 3; ++r) {
                if (!std::isfinite(bridge.weights[r]) || bridge.weights[r] < 0.0 ||
                    !std::isfinite(bridge.branch_lengths[r]) || bridge.branch_lengths[r] < 0.0)
                    throw std::invalid_argument("difFUBAR: invalid Wheatstone bridge parameters");
                weight_sum += bridge.weights[r];
                if (bridge.children[r] < 0 || static_cast<std::size_t>(bridge.children[r]) >= N || !visited[bridge.children[r]])
                    throw std::invalid_argument("difFUBAR: invalid Wheatstone child");
            }
            if (std::abs(weight_sum - 1.0) > 1e-8)
                throw std::invalid_argument("difFUBAR: Wheatstone weights must sum to one");
        }
        modeled_edges_.assign(edges.begin(), edges.end());
        bool positive_length = false;
        for (const auto id : modeled_edges_) {
            const double length = tree.nodes[id].branch_length;
            if (!(length >= 0.0) || !std::isfinite(length))
                throw std::invalid_argument("difFUBAR: branch lengths must be finite and nonnegative");
            positive_length |= length > 0.0;
        }
        for (const auto& b : tree.bridges) for (double length : b.branch_lengths) positive_length |= length > 0.0;
        if (!positive_length) throw std::invalid_argument("difFUBAR: all branch lengths are zero/missing");
        std::vector<std::size_t> weights(aln.patterns.size(), 0);
        for (const auto p : aln.site_to_pattern) {
            if (p >= weights.size()) throw std::invalid_argument("difFUBAR: invalid site-to-pattern mapping");
            ++weights[p];
        }
        for (std::size_t p = 0; p < weights.size(); ++p)
            if (weights[p] == 0 || weights[p] != aln.patterns[p].weight || aln.patterns[p].states.size() != aln.num_taxa)
                throw std::invalid_argument("difFUBAR: inconsistent compressed pattern weights/states");
    }
    void prepare_groups() {
        std::set<std::string> tags;
        for (const auto id : modeled_edges_)
            if (!tree.nodes[id].model_tag.empty()) tags.insert(tree.nodes[id].model_tag);
        std::vector<std::string> selected = settings.groups;
        if (selected.empty()) {
            if (tags.size() != 2)
                throw std::invalid_argument("difFUBAR: annotate two foreground groups with {tags}, or select two existing tags using --groups");
            selected.assign(tags.begin(), tags.end());
        }
        groups.names = {selected[0], selected[1]};
        groups.node_group.assign(tree.num_nodes(), -1);
        for (const auto id : modeled_edges_) {
            const auto& tag = tree.nodes[id].model_tag;
            const int group = tag == selected[0] ? 0 : tag == selected[1] ? 1 : 2;
            groups.node_group[id] = group;
            ++groups.branch_counts[group];
        }
        if (groups.branch_counts[0] == 0 || groups.branch_counts[1] == 0)
            throw std::invalid_argument("difFUBAR: each foreground tag must label at least one modeled non-root branch");
        groups.has_background = groups.branch_counts[2] > 0;
    }
    template<class Model>
    double homogeneous_log_likelihood(const core::Alignment& data, const Model& model,
                                      const core::Vector& pi, double scale) const {
        const int S = static_cast<int>(pi.size());
        std::vector<core::Matrix> transitions(tree.num_nodes(), core::Matrix::Identity(S, S));
        for (const auto id : tree.post_order)
            if (id != tree.root_id)
                transitions[id] = diffubar_detail::transition(model, tree.nodes[id].branch_length * scale, S);
        std::vector<std::array<core::Matrix, 3>> bridges(tree.bridges.size());
        for (std::size_t b = 0; b < tree.bridges.size(); ++b)
            for (std::size_t r = 0; r < 3; ++r)
                bridges[b][r] = diffubar_detail::transition(model, tree.bridges[b].branch_lengths[r] * scale, S);
        diffubar_detail::PruningWorkspace workspace(tree.num_nodes(), S);
        double log_likelihood = 0.0;
        for (const auto& pattern : data.patterns) {
            const double log_l = diffubar_detail::pattern_log_likelihood(tree, pattern, leaf_to_taxon_, pi,
                                                                         transitions, bridges, workspace);
            if (!std::isfinite(log_l)) return -std::numeric_limits<double>::infinity();
            log_likelihood += static_cast<double>(pattern.weight) * log_l;
        }
        return log_likelihood;
    }
};

} // namespace hyphy::analyses
