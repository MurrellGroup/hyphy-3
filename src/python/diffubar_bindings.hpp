#pragma once
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include "hyphy/analyses/diffubar_json.hpp"

inline void bind_diffubar(nanobind::module_& module) {
    namespace nb = nanobind;
    using namespace nb::literals;
    using namespace hyphy::analyses;
    nb::class_<DifFUBARSamplerSettings>(module, "DifFUBARSamplerSettings")
        .def(nb::init<>())
        .def_rw("iterations", &DifFUBARSamplerSettings::iterations)
        .def_rw("burn_in", &DifFUBARSamplerSettings::burn_in)
        .def_rw("chains", &DifFUBARSamplerSettings::chains)
        .def_rw("seed", &DifFUBARSamplerSettings::seed)
        .def_rw("concentration", &DifFUBARSamplerSettings::concentration);
    nb::class_<DifFUBARSettings>(module, "DifFUBARSettings")
        .def(nb::init<>())
        .def_rw("sampler", &DifFUBARSettings::sampler)
        .def_rw("foreground_grid", &DifFUBARSettings::foreground_grid)
        .def_rw("background_grid", &DifFUBARSettings::background_grid)
        .def_rw("groups", &DifFUBARSettings::groups)
        .def_rw("posterior_threshold", &DifFUBARSettings::posterior_threshold)
        .def_rw("threads", &DifFUBARSettings::threads)
        .def_rw("max_memory_mb", &DifFUBARSettings::max_memory_mb)
        .def_rw("fit_max_iterations", &DifFUBARSettings::fit_max_iterations)
        .def_rw("fit_tolerance", &DifFUBARSettings::fit_tolerance);
    module.def("diffubar", [](hyphy::core::Tree tree, hyphy::core::Alignment alignment,
                              DifFUBARSettings settings, const std::string& fitted_model_json) -> nb::object {
        std::string serialized;
        {
            nb::gil_scoped_release release;
            DifFUBARAnalyzer analyzer(std::move(tree), std::move(alignment), std::move(settings));
            if (fitted_model_json.empty()) analyzer.fit_global_model();
            else analyzer.use_fitted_model(diffubar_fitted_model_from_json(nlohmann::json::parse(fitted_model_json)));
            serialized = diffubar_to_json(analyzer, analyzer.run()).dump();
        }
        return nb::module_::import_("json").attr("loads")(serialized);
    }, "tree"_a, "alignment"_a, "settings"_a = DifFUBARSettings{}, "fitted_model_json"_a = "",
       "Run native difFUBAR and return a result dictionary. An optional fitted-model JSON string bypasses global fitting.");
}
