#pragma once

// Numerical difFUBAR machinery. Deliberately independent of Eigen and the
// phylogenetic engine so the Bayesian calculation can be tested in isolation.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace hyphy::analyses {

struct DifFUBARCategory {
    double alpha = 1.0;
    std::array<double, 3> omega{1.0, 1.0, 1.0};
    std::size_t alpha_index = 0;
    std::array<std::size_t, 3> omega_index{0, 0, 0};
};

struct DifFUBARGrid {
    std::vector<double> alpha;
    std::vector<double> omega;
    std::vector<double> background_omega;
    std::vector<DifFUBARCategory> categories;
    bool has_background = false;
};

struct DifFUBARSamplerSettings {
    std::size_t iterations = 2500;
    // -1 means floor(iterations / 5), including when iterations is changed.
    std::int64_t burn_in = -1;
    std::size_t chains = 1;
    double concentration = 0.1; // per category, NOT total concentration
    std::uint64_t seed = 1;

    std::size_t resolved_burn_in() const {
        if (burn_in < -1) throw std::invalid_argument("difFUBAR: burn-in must be -1 or nonnegative");
        return burn_in == -1 ? iterations / 5 : static_cast<std::size_t>(burn_in);
    }
    void validate() const {
        if (iterations == 0 || chains == 0 || resolved_burn_in() >= iterations)
            throw std::invalid_argument("difFUBAR: need positive chains and iterations > burn-in");
        if (!(concentration > 0.0) || !std::isfinite(concentration))
            throw std::invalid_argument("difFUBAR: concentration must be finite and positive");
        if (chains > std::numeric_limits<std::size_t>::max() / iterations)
            throw std::invalid_argument("difFUBAR: total iteration count overflows size_t");
    }
};

struct DifFUBARPosteriorSummary {
    double mean_alpha = 0.0;
    std::array<double, 3> mean_omega{0.0, 0.0, 0.0};
    double prob_omega1_gt_omega2 = 0.0;
    double prob_omega2_gt_omega1 = 0.0;
    double prob_omega1_eq_omega2 = 0.0;
    std::array<double, 2> prob_positive{0.0, 0.0};
    std::vector<double> alpha_marginal;
    std::array<std::vector<double>, 3> omega_marginal;
};

struct DifFUBARPosterior {
    std::vector<DifFUBARPosteriorSummary> patterns;
    std::vector<double> mean_grid_weights;
    // Per-chain means of the four reported probabilities, for comparisons
    // between independent chains. These are not an ESS/R-hat estimate.
    std::vector<std::vector<std::array<double, 4>>> chain_probabilities;
    std::size_t retained_draws = 0;
};

namespace diffubar_detail {

inline std::size_t checked_product(std::size_t a, std::size_t b, const char* what) {
    if (b != 0 && a > std::numeric_limits<std::size_t>::max() / b)
        throw std::length_error(std::string("difFUBAR: size overflow: ") + what);
    return a * b;
}

inline std::vector<double> transformed_grid(double lower, double upper, int below_one) {
    if (!(lower > 0.0 && lower < 1.0 && upper > 1.0) ||
        !std::isfinite(upper) || below_one <= 0 || below_one > 100)
        throw std::invalid_argument("difFUBAR: invalid grid bounds/density (allowed density 1..100)");
    const double start = std::log10(lower + 0.05);
    const double step = (std::log10(1.05) - start) / below_one;
    const auto count = static_cast<std::size_t>(
        std::floor((std::log10(upper + 0.05) - start) / step + 1e-12)) + 1;
    std::vector<double> values(count);
    for (std::size_t i = 0; i < count; ++i)
        values[i] = std::pow(10.0, start + step * static_cast<double>(i)) - 0.05;
    // The category at the neutral boundary is exactly neutral, independent of
    // libm rounding; it must not leak into P(omega > 1).
    values.at(static_cast<std::size_t>(below_one)) = 1.0;
    return values;
}

inline void initialize_summary(DifFUBARPosteriorSummary& s, const DifFUBARGrid& grid) {
    s.alpha_marginal.assign(grid.alpha.size(), 0.0);
    s.omega_marginal[0].assign(grid.omega.size(), 0.0);
    s.omega_marginal[1].assign(grid.omega.size(), 0.0);
    if (grid.has_background) s.omega_marginal[2].assign(grid.background_omega.size(), 0.0);
}

inline void accumulate(DifFUBARPosteriorSummary& s, const DifFUBARCategory& c,
                       double probability, bool background) {
    s.mean_alpha += probability * c.alpha;
    s.alpha_marginal[c.alpha_index] += probability;
    for (std::size_t j = 0; j < (background ? 3u : 2u); ++j) {
        s.mean_omega[j] += probability * c.omega[j];
        s.omega_marginal[j][c.omega_index[j]] += probability;
    }
    if (c.omega[0] > c.omega[1]) s.prob_omega1_gt_omega2 += probability;
    else if (c.omega[1] > c.omega[0]) s.prob_omega2_gt_omega1 += probability;
    else s.prob_omega1_eq_omega2 += probability;
    for (std::size_t j = 0; j < 2; ++j)
        if (c.omega[j] > 1.0) s.prob_positive[j] += probability;
}

inline void divide_summary(DifFUBARPosteriorSummary& s, double n) {
    s.mean_alpha /= n;
    for (double& x : s.mean_omega) x /= n;
    s.prob_omega1_gt_omega2 /= n;
    s.prob_omega2_gt_omega1 /= n;
    s.prob_omega1_eq_omega2 /= n;
    for (double& x : s.prob_positive) x /= n;
    for (double& x : s.alpha_marginal) x /= n;
    for (auto& marginal : s.omega_marginal) for (double& x : marginal) x /= n;
}

} // namespace diffubar_detail

inline DifFUBARGrid create_diffubar_grid(bool background, int foreground_grid = 6,
                                        int background_grid = 4,
                                        std::size_t max_categories = 1000000) {
    DifFUBARGrid grid;
    grid.has_background = background;
    grid.alpha = diffubar_detail::transformed_grid(0.01, 13.0, foreground_grid);
    grid.omega = grid.alpha;
    grid.background_omega = diffubar_detail::transformed_grid(0.05, 6.0, background_grid);
    auto count = diffubar_detail::checked_product(grid.alpha.size(), grid.omega.size(), "grid");
    count = diffubar_detail::checked_product(count, grid.omega.size(), "grid");
    count = diffubar_detail::checked_product(count, background ? grid.background_omega.size() : 1, "grid");
    if (count > max_categories) throw std::length_error("difFUBAR: grid exceeds category limit");
    grid.categories.reserve(count);
    // Same order as EvoOnline: background varies fastest, alpha slowest.
    for (std::size_t a = 0; a < grid.alpha.size(); ++a)
        for (std::size_t i = 0; i < grid.omega.size(); ++i)
            for (std::size_t j = 0; j < grid.omega.size(); ++j)
                for (std::size_t b = 0; b < (background ? grid.background_omega.size() : 1); ++b)
                    grid.categories.push_back({grid.alpha[a],
                        {grid.omega[i], grid.omega[j], background ? grid.background_omega[b] : 1.0},
                        a, {i, j, b}});
    return grid;
}

// log_likelihoods is pattern-major [pattern * G + category]. Multiplicities
// represent independent observations, not powers of the conditional likelihood
// of one shared allocation. All latent allocations are sampled independently.
// Rao-Blackwellization averages p(z_site | theta, D_site) over retained theta
// draws, never p(z_site | E[theta], D_site).
inline DifFUBARPosterior sample_diffubar(
    const DifFUBARGrid& grid,
    const std::vector<double>& log_likelihoods,
    const std::vector<std::size_t>& multiplicities,
    DifFUBARSamplerSettings settings = {},
    const std::function<void(std::size_t, std::size_t)>& progress = {}) {
    settings.validate();
    const std::size_t G = grid.categories.size(), P = multiplicities.size();
    if (G == 0 || P == 0 || log_likelihoods.size() != diffubar_detail::checked_product(P, G, "likelihoods"))
        throw std::invalid_argument("difFUBAR: empty or inconsistent conditional-likelihood matrix");
    std::size_t sites = 0;
    for (std::size_t n : multiplicities) {
        if (n == 0 || n > std::numeric_limits<std::size_t>::max() - sites)
            throw std::invalid_argument("difFUBAR: invalid pattern multiplicity");
        sites += n;
    }
    if (sites > (std::uint64_t{1} << 53))
        throw std::length_error("difFUBAR: number of sites exceeds exact double integer range");
    for (const auto& c : grid.categories) {
        if (!(c.alpha > 0.0) || !std::isfinite(c.alpha) || c.alpha_index >= grid.alpha.size())
            throw std::invalid_argument("difFUBAR: invalid alpha category");
        for (std::size_t j = 0; j < (grid.has_background ? 3u : 2u); ++j) {
            const auto& axis = j == 2 ? grid.background_omega : grid.omega;
            if (!(c.omega[j] > 0.0) || !std::isfinite(c.omega[j]) || c.omega_index[j] >= axis.size())
                throw std::invalid_argument("difFUBAR: invalid omega category");
        }
    }
    std::vector<double> conditionals(log_likelihoods.size());
    for (std::size_t p = 0; p < P; ++p) {
        double max_log = -std::numeric_limits<double>::infinity();
        for (std::size_t g = 0; g < G; ++g) {
            const double x = log_likelihoods[p * G + g];
            if (std::isnan(x) || x == std::numeric_limits<double>::infinity())
                throw std::invalid_argument("difFUBAR: non-finite conditional likelihood");
            max_log = std::max(max_log, x);
        }
        if (!std::isfinite(max_log))
            throw std::runtime_error("difFUBAR: every grid category is impossible for a site pattern");
        for (std::size_t g = 0; g < G; ++g)
            conditionals[p * G + g] = std::exp(log_likelihoods[p * G + g] - max_log);
    }
    DifFUBARPosterior result;
    result.patterns.resize(P);
    for (auto& s : result.patterns) diffubar_detail::initialize_summary(s, grid);
    result.mean_grid_weights.assign(G, 0.0);
    result.chain_probabilities.resize(settings.chains,
        std::vector<std::array<double, 4>>(P, {0.0, 0.0, 0.0, 0.0}));
    const std::size_t burn = settings.resolved_burn_in();
    result.retained_draws = (settings.iterations - burn) * settings.chains;
    std::vector<std::size_t> counts(G), next_counts(G);
    std::vector<double> theta(G), weights(G), cdf(G);
    for (std::size_t chain = 0; chain < settings.chains; ++chain) {
        std::mt19937_64 rng(settings.seed + std::uint64_t{0x9e3779b97f4a7c15} * chain);
        std::fill(counts.begin(), counts.end(), 0);
        for (std::size_t iter = 0; iter < settings.iterations; ++iter) {
            double theta_sum = 0.0;
            for (std::size_t g = 0; g < G; ++g) {
                std::gamma_distribution<double> gamma(settings.concentration + static_cast<double>(counts[g]), 1.0);
                theta[g] = gamma(rng);
                theta_sum += theta[g];
            }
            if (!(theta_sum > 0.0) || !std::isfinite(theta_sum))
                throw std::runtime_error("difFUBAR: Dirichlet draw underflow/overflow; use a numerically resolvable concentration");
            for (double& w : theta) w /= theta_sum;
            const bool retain = iter >= burn;
            if (retain)
                for (std::size_t g = 0; g < G; ++g) result.mean_grid_weights[g] += theta[g];
            std::fill(next_counts.begin(), next_counts.end(), 0);
            for (std::size_t p = 0; p < P; ++p) {
                double denominator = 0.0;
                for (std::size_t g = 0; g < G; ++g) {
                    weights[g] = theta[g] * conditionals[p * G + g];
                    denominator += weights[g];
                    cdf[g] = denominator;
                }
                if (!(denominator > 0.0) || !std::isfinite(denominator))
                    throw std::runtime_error("difFUBAR: posterior normalization underflow; no supported category in this draw");
                if (retain) {
                    auto& per_chain = result.chain_probabilities[chain][p];
                    for (std::size_t g = 0; g < G; ++g) {
                        const double probability = weights[g] / denominator;
                        const auto& c = grid.categories[g];
                        diffubar_detail::accumulate(result.patterns[p], c, probability, grid.has_background);
                        if (c.omega[0] > c.omega[1]) per_chain[0] += probability;
                        if (c.omega[1] > c.omega[0]) per_chain[1] += probability;
                        if (c.omega[0] > 1.0) per_chain[2] += probability;
                        if (c.omega[1] > 1.0) per_chain[3] += probability;
                    }
                }
                // Build one CDF per pattern, but sample one independent latent
                // category for EACH occurrence in the uncompressed alignment.
                for (std::size_t n = 0; n < multiplicities[p]; ++n) {
                    const double u = std::min(std::generate_canonical<double, 53>(rng) * denominator,
                                              std::nextafter(denominator, 0.0));
                    const auto it = std::upper_bound(cdf.begin(), cdf.end(), u);
                    const std::size_t g = static_cast<std::size_t>(it - cdf.begin());
                    if (g >= G) throw std::runtime_error("difFUBAR: invalid categorical draw");
                    ++next_counts[g];
                }
            }
            counts.swap(next_counts);
            if (progress && (iter % 25 == 0 || iter + 1 == settings.iterations))
                progress(chain * settings.iterations + iter + 1, settings.chains * settings.iterations);
        }
        for (auto& probabilities : result.chain_probabilities[chain])
            for (double& x : probabilities) x /= static_cast<double>(settings.iterations - burn);
    }
    const double retained = static_cast<double>(result.retained_draws);
    for (auto& s : result.patterns) diffubar_detail::divide_summary(s, retained);
    for (double& x : result.mean_grid_weights) x /= retained;
    return result;
}

} // namespace hyphy::analyses
