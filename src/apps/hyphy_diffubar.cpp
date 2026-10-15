#include "hyphy/analyses/diffubar_json.hpp"
#include <filesystem>
#include <iostream>
#include <sstream>

using namespace hyphy::analyses;
using namespace hyphy::core;

namespace {
void usage(const char* program) {
    std::cout << "Usage: " << program << " --alignment FILE [--tree FILE] [OPTIONS]\n"
        "difFUBAR: Bayesian differential codon selection between two tagged branch groups.\n\n"
        "  --alignment, -a FILE     Codon FASTA/NEXUS alignment\n"
        "  --tree, -t FILE          Tagged Newick; otherwise use the alignment's embedded tree\n"
        "  --groups A,B            Select/order two exact tags; others are background\n"
        "                          Default: two distinct non-root tags in lexical order\n"
        "  --code NAME             Native genetic code name (default: Universal)\n"
        "  --foreground-grid N     Points below 1; default 6 gives 12 values per foreground axis\n"
        "  --background-grid N     Default 4 gives 7 background omega values\n"
        "  --iterations N          Gibbs iterations per chain (default: 2500)\n"
        "  --burn-in N             Default: floor(iterations/5)\n"
        "  --chains N              Independent chains (default: 1)\n"
        "  --seed N                Unsigned 64-bit RNG seed (default: 1)\n"
        "  --concentration X       Dirichlet concentration PER CATEGORY (default: 0.1)\n"
        "  --posterior X           Reporting threshold (default: 0.95)\n"
        "  --threads, -T N         Grid worker threads (0: runtime default)\n"
        "  --max-memory-mb N       Dominant likelihood/cache storage budget (default: 1024)\n"
        "  --fitted-model FILE     Skip fitting; read fixed-model JSON or a prior result JSON\n"
        "  --save-fitted-model FILE  Export fitted parameters for reproducible comparisons\n"
        "  --grid-likelihoods FILE  Export category-by-site log likelihoods as TSV\n"
        "  --output, -o FILE       JSON output (default: ALIGNMENT.difFUBAR.json)\n"
        "  --quiet                 Suppress status messages\n"
        "  --help, -h              Show help\n\n"
        "Tags label incoming branches, e.g. ((A{g1}:0.1,B{g1}:0.1){g1}:0.1,\n"
        "                                (C{g2}:0.1,D{g2}:0.1){g2}:0.1);\n"
        "Input branch lengths are retained; the global fitted alpha rescales them for the grid.\n";
}
std::uint64_t integer(const std::string& value, const std::string& option) {
    if (value.empty() || value.front() == '-' || value.front() == '+')
        throw std::invalid_argument(option + " requires an unsigned integer");
    std::size_t consumed = 0;
    const auto number = std::stoull(value, &consumed);
    if (consumed != value.size()) throw std::invalid_argument("Invalid integer for " + option);
    return number;
}
int small_integer(const std::string& value, const std::string& option) {
    const auto number = integer(value, option);
    if (number > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Integer out of range for " + option);
    return static_cast<int>(number);
}
double real(const std::string& value, const std::string& option) {
    std::size_t consumed = 0;
    const auto number = std::stod(value, &consumed);
    if (consumed != value.size() || !std::isfinite(number)) throw std::invalid_argument("Invalid number for " + option);
    return number;
}
nlohmann::json read_json(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot read " + path);
    return nlohmann::json::parse(input);
}
void write_json(const nlohmann::json& j, const std::string& path) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Cannot write " + path);
    output << j.dump(2) << '\n';
    if (!output) throw std::runtime_error("Failed while writing " + path);
}
std::filesystem::path normalized_path(const std::string& path) {
    return std::filesystem::weakly_canonical(std::filesystem::absolute(path));
}
}

int run_diffubar(int argc, char* argv[]) {
    try {
        DifFUBARSettings settings;
        std::string alignment_path, tree_path, output_path, model_path, saved_model_path, grid_path;
        std::string code_name = "Universal";
        bool quiet = false;
        for (int i = 1; i < argc; ++i) {
            const std::string option = argv[i];
            if (option == "--help" || option == "-h") { usage(argv[0]); return 0; }
            if (option == "--quiet") { quiet = true; continue; }
            if (i + 1 == argc) throw std::invalid_argument("Missing value for " + option);
            const std::string value = argv[++i];
            if (option == "--alignment" || option == "-a") alignment_path = value;
            else if (option == "--tree" || option == "-t") tree_path = value;
            else if (option == "--output" || option == "-o") output_path = value;
            else if (option == "--code") code_name = value;
            else if (option == "--fitted-model") model_path = value;
            else if (option == "--save-fitted-model") saved_model_path = value;
            else if (option == "--grid-likelihoods") { grid_path = value; settings.store_grid_likelihoods = true; }
            else if (option == "--groups") {
                const auto comma = value.find(',');
                if (comma == std::string::npos || value.find(',', comma + 1) != std::string::npos)
                    throw std::invalid_argument("--groups requires two comma-separated exact tag names");
                settings.groups = {value.substr(0, comma), value.substr(comma + 1)};
            } else if (option == "--foreground-grid") settings.foreground_grid = small_integer(value, option);
            else if (option == "--background-grid") settings.background_grid = small_integer(value, option);
            else if (option == "--iterations") settings.sampler.iterations = integer(value, option);
            else if (option == "--chains") settings.sampler.chains = integer(value, option);
            else if (option == "--seed") settings.sampler.seed = integer(value, option);
            else if (option == "--burn-in") {
                const auto burn = integer(value, option);
                if (burn > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
                    throw std::invalid_argument("--burn-in is too large");
                settings.sampler.burn_in = static_cast<std::int64_t>(burn);
            } else if (option == "--concentration") settings.sampler.concentration = real(value, option);
            else if (option == "--posterior") settings.posterior_threshold = real(value, option);
            else if (option == "--threads" || option == "-T") settings.threads = small_integer(value, option);
            else if (option == "--max-memory-mb") settings.max_memory_mb = integer(value, option);
            else throw std::invalid_argument("Unknown option: " + option);
        }
        if (alignment_path.empty()) throw std::invalid_argument("--alignment is required");
        settings.validate();
        if (output_path.empty()) output_path = alignment_path + ".difFUBAR.json";
        // Do not overwrite inputs or collide result/model/grid output files.
        std::set<std::filesystem::path> outputs;
        for (const auto& path : {output_path, saved_model_path, grid_path}) {
            if (path.empty()) continue;
            const auto normalized = normalized_path(path);
            if (!outputs.insert(normalized).second) throw std::invalid_argument("Output paths must be distinct");
            for (const auto& input : {alignment_path, tree_path, model_path})
                if (!input.empty() && normalized == normalized_path(input))
                    throw std::invalid_argument("An output path would overwrite an input file");
        }
        const auto code = GeneticCode::from_name(code_name);
        auto alignment = Alignment::load(alignment_path, code);
        Tree tree;
        if (!tree_path.empty()) tree = Tree::from_newick_file(tree_path);
        else if (!alignment.embedded_tree_newick.empty()) tree = Tree::from_newick(alignment.embedded_tree_newick);
        else throw std::invalid_argument("A tagged --tree or embedded NEXUS tree is required; foreground groups cannot be inferred from an alignment alone");
        DifFUBARAnalyzer analyzer(std::move(tree), std::move(alignment), settings);
        std::string last_stage;
        const auto progress = [&](const std::string& stage, double) {
            if (!quiet && stage != last_stage) { std::cerr << "[difFUBAR] " << stage << '\n'; last_stage = stage; }
        };
        if (!model_path.empty()) analyzer.use_fitted_model(diffubar_fitted_model_from_json(read_json(model_path)));
        else analyzer.fit_global_model(progress);
        const auto result = analyzer.run(progress);
        auto json = diffubar_to_json(analyzer, result);
        json["input"]["alignment_file"] = alignment_path;
        json["input"]["tree_file"] = tree_path;
        if (!saved_model_path.empty()) write_json(diffubar_fitted_model_to_json(result.fitted_model), saved_model_path);
        if (!grid_path.empty()) write_diffubar_grid_tsv(result, grid_path);
        write_json(json, output_path);
        if (!quiet) std::cerr << "[difFUBAR] Wrote " << result.site_results.size() << " codon sites to " << output_path << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "difFUBAR: " << error.what() << '\n';
        return 1;
    }
}

#ifndef HYPHY3_COMBINED_DRIVER
int main(int argc, char* argv[]) { return run_diffubar(argc, argv); }
#endif
