#include "hyphy/analyses/detail/diffubar_inference.hpp"
#include <iostream>
#include <sstream>

using namespace hyphy::analyses;

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void near(double actual, double expected, double tolerance, const char* message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::ostringstream out;
        out << message << ": actual=" << actual << ", expected=" << expected;
        throw std::runtime_error(out.str());
    }
}
template<class F> void rejects(F&& f) {
    bool rejected = false;
    try { f(); } catch (const std::exception&) { rejected = true; }
    require(rejected, "invalid input was not rejected");
}
DifFUBARGrid two_categories() {
    DifFUBARGrid grid;
    grid.alpha = {1.0};
    grid.omega = {0.5, 2.0};
    grid.categories = {{1.0, {0.5, 2.0, 1.0}, 0, {0, 1, 0}},
                       {1.0, {2.0, 0.5, 1.0}, 0, {1, 0, 0}}};
    return grid;
}
// Independent, exact integration over theta for a two-category, S-site model:
// p(z,D) = prod_s C[s,z_s] * B(n_0+a,n_1+a)/B(a,a).
std::vector<double> enumerate(const std::vector<std::array<double, 2>>& C, double a) {
    const std::size_t S = C.size();
    std::vector<double> probabilities(S, 0.0);
    double denominator = 0.0;
    for (std::size_t mask = 0; mask < (std::size_t{1} << S); ++mask) {
        std::size_t n1 = 0;
        double log_weight = 0.0;
        for (std::size_t s = 0; s < S; ++s) {
            const std::size_t z = (mask >> s) & 1;
            n1 += z;
            log_weight += std::log(C[s][z]);
        }
        log_weight += std::lgamma(a + static_cast<double>(n1)) +
            std::lgamma(a + static_cast<double>(S - n1)) -
            std::lgamma(2.0 * a + static_cast<double>(S));
        const double weight = std::exp(log_weight);
        denominator += weight;
        for (std::size_t s = 0; s < S; ++s)
            if ((mask >> s) & 1) probabilities[s] += weight;
    }
    for (double& x : probabilities) x /= denominator;
    return probabilities;
}
}

int main() {
    try {
        const auto default_grid = create_diffubar_grid(false);
        const auto background_grid = create_diffubar_grid(true);
        require(default_grid.categories.size() == 1728, "default foreground grid size");
        require(background_grid.categories.size() == 12096, "default background grid size");
        require(default_grid.alpha.size() == 12 && background_grid.background_omega.size() == 7,
                "default axis lengths");
        near(default_grid.alpha.front(), 0.01, 1e-15, "lower endpoint");
        near(default_grid.alpha[6], 1.0, 0.0, "exact neutral boundary");
        near(background_grid.background_omega[4], 1.0, 0.0, "background neutral boundary");
        require(background_grid.categories[1].omega_index == std::array<std::size_t, 3>{0, 0, 1},
                "background is fastest-varying coordinate");
        require(default_grid.categories[12].omega_index == std::array<std::size_t, 3>{1, 0, 0},
                "foreground category order");
        require(default_grid.categories[144].alpha_index == 1, "alpha is slowest coordinate");
        rejects([] { create_diffubar_grid(false, 0); });
        rejects([] { create_diffubar_grid(true, 100); });
        rejects([] { create_diffubar_grid(true, 6, -1); });

        auto grid = two_categories();
        const std::vector<double> logs{std::log(0.1), 0.0, 0.0, std::log(0.3)};
        DifFUBARSamplerSettings settings;
        settings.iterations = 70000;
        settings.burn_in = 10000;
        settings.chains = 2;
        settings.seed = 572;
        settings.concentration = 0.3;
        const auto posterior = sample_diffubar(grid, logs, {2, 1}, settings);
        const auto exact = enumerate({{0.1, 1.0}, {0.1, 1.0}, {1.0, 0.3}}, settings.concentration);
        near(posterior.patterns[0].prob_omega1_gt_omega2, exact[0], 0.009, "compressed pattern posterior");
        near(posterior.patterns[1].prob_omega1_gt_omega2, exact[2], 0.009, "other pattern posterior");
        // A one-allocation-per-pattern bug gives a materially different answer.
        const auto wrong = enumerate({{0.1, 1.0}, {1.0, 0.3}}, settings.concentration);
        require(std::abs(exact[2] - wrong[1]) > 0.08, "fixture does not detect lost multiplicities");
        require(posterior.retained_draws == 120000, "burn-in/chain retained count");
        near(std::accumulate(posterior.mean_grid_weights.begin(), posterior.mean_grid_weights.end(), 0.0),
             1.0, 1e-11, "theta simplex");
        for (const auto& summary : posterior.patterns) {
            near(summary.prob_omega1_gt_omega2 + summary.prob_omega2_gt_omega1 + summary.prob_omega1_eq_omega2,
                 1.0, 1e-11, "comparison probabilities sum to one");
            near(summary.prob_positive[0], summary.prob_omega1_gt_omega2, 1e-13, "positive-selection marginal");
            near(summary.mean_omega[0], 0.5 + 1.5 * summary.prob_positive[0], 1e-10, "omega posterior mean");
            near(std::accumulate(summary.omega_marginal[0].begin(), summary.omega_marginal[0].end(), 0.0),
                 1.0, 1e-11, "omega marginal simplex");
        }
        const auto expanded = sample_diffubar(grid,
            {logs[0], logs[1], logs[0], logs[1], logs[2], logs[3]}, {1, 1, 1}, settings);
        near(posterior.patterns[0].prob_positive[0], expanded.patterns[0].prob_positive[0], 1e-14,
             "compression must preserve seeded independent allocations");
        near(expanded.patterns[0].prob_positive[0], expanded.patterns[1].prob_positive[0], 0.0,
             "identical patterns have identical Rao-Blackwellized posteriors");

        settings.iterations = 2000; settings.burn_in = -1; settings.chains = 1;
        auto repeat1 = sample_diffubar(grid, logs, {2, 1}, settings);
        auto repeat2 = sample_diffubar(grid, logs, {2, 1}, settings);
        require(repeat1.mean_grid_weights == repeat2.mean_grid_weights, "seed reproducibility");
        require(repeat1.retained_draws == 1600, "automatic burn-in follows iteration count");
        for (auto& category : grid.categories) {
            std::swap(category.omega[0], category.omega[1]);
            std::swap(category.omega_index[0], category.omega_index[1]);
        }
        auto swapped = sample_diffubar(grid, logs, {2, 1}, settings);
        near(repeat1.patterns[0].prob_omega1_gt_omega2, swapped.patterns[0].prob_omega2_gt_omega1,
             0.0, "group swap symmetry");
        grid = two_categories();
        auto shifted_logs = logs;
        for (double& x : shifted_logs) x -= 10000.0;
        auto shifted = sample_diffubar(grid, shifted_logs, {2, 1}, settings);
        near(shifted.patterns[0].prob_positive[0], repeat1.patterns[0].prob_positive[0],
             1e-10, "log likelihood offset invariance/underflow protection");
        DifFUBARGrid neutral;
        neutral.alpha = {1.0}; neutral.omega = {1.0}; neutral.categories = {DifFUBARCategory{}};
        auto tied = sample_diffubar(neutral, {0.0}, {1}, settings);
        near(tied.patterns[0].prob_omega1_eq_omega2, 1.0, 0.0, "neutral equality mass");
        near(tied.patterns[0].prob_positive[0], 0.0, 0.0, "neutral is not positive selection");
        auto invalid = settings;
        invalid.concentration = 0.0;
        rejects([&] { sample_diffubar(grid, logs, {2, 1}, invalid); });
        invalid = settings; invalid.burn_in = 2000;
        rejects([&] { sample_diffubar(grid, logs, {2, 1}, invalid); });
        invalid = settings; invalid.burn_in = -2;
        rejects([&] { sample_diffubar(grid, logs, {2, 1}, invalid); });
        rejects([&] { sample_diffubar(grid, logs, {0, 1}, settings); });
        rejects([&] { sample_diffubar(grid, logs, {1}, settings); });
        rejects([&] { sample_diffubar(grid, {NAN, 0.0}, {1}, settings); });
        rejects([&] { sample_diffubar(grid, {INFINITY, 0.0}, {1}, settings); });
        rejects([&] { sample_diffubar(grid, {-INFINITY, -INFINITY}, {1}, settings); });
        auto partially_impossible = sample_diffubar(grid, {-INFINITY, 0.0}, {1}, settings);
        near(partially_impossible.patterns[0].prob_positive[0], 1.0, 0.0, "impossible category has zero mass");
        std::cout << "difFUBAR inference tests passed: grids, exact posterior enumeration, pattern weights,\n"
                     "marginals, group symmetry, deterministic seeds, neutral boundary, numerical scaling, validation.\n";
        std::cout << "Exact P(omega1>omega2): " << exact[0] << ", " << exact[2] << '\n';
        std::cout << "Gibbs P(omega1>omega2): " << posterior.patterns[0].prob_omega1_gt_omega2
                  << ", " << posterior.patterns[1].prob_omega1_gt_omega2 << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
