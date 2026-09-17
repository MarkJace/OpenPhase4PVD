#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "DrivingForce.h"
#include "PhaseField.h"

namespace py = pybind11;
using openphase::DrivingForce;
using openphase::PhaseField;

namespace {

std::pair<int, int> planar_shape(const PhaseField& phase) {
    return {static_cast<int>(phase.Fields.sizeY()),
            static_cast<int>(phase.Fields.sizeX())};
}

double phase_fraction_at(const PhaseField& phase, int i, int j,
                         std::size_t phase_index) {
    double value = 0.0;
    const auto& node = phase.Fields(i, j, 0);
    for (auto field = node.cbegin(); field != node.cend(); ++field) {
        if (phase.FieldsProperties[field->index].Phase == phase_index) {
            value += field->value;
        }
    }
    return value;
}

struct FluxParameters {
    double ux;
    double uy;
    double peak_driving;
    double partition_power;
};

FluxParameters validate_flux(double flux_x, double flux_y,
                             double flux_to_driving,
                             double partition_power) {
    if (flux_to_driving < 0.0) {
        throw std::invalid_argument("flux_to_driving must be non-negative");
    }
    if (partition_power <= 0.0) {
        throw std::invalid_argument("partition_power must be positive");
    }
    const double magnitude = std::hypot(flux_x, flux_y);
    if (magnitude == 0.0) {
        return {0.0, 0.0, 0.0, partition_power};
    }
    return {flux_x / magnitude, flux_y / magnitude,
            magnitude * flux_to_driving, partition_power};
}

double exposure_at(const PhaseField& phase, int i, int j,
                   std::size_t vapor_phase, const FluxParameters& flux) {
    const auto [ny, nx] = planar_shape(phase);
    const int im = (i + nx - 1) % nx;
    const int ip = (i + 1) % nx;
    const int jm = std::max(j - 1, 0);
    const int jp = std::min(j + 1, ny - 1);

    // Only the gradient direction is needed.  Omitting dx here avoids a
    // numerically pointless common factor in the normalization.
    const double gx = 0.5 * (phase_fraction_at(phase, ip, j, vapor_phase) -
                             phase_fraction_at(phase, im, j, vapor_phase));
    const double gy_denom = static_cast<double>(std::max(jp - jm, 1));
    const double gy = (phase_fraction_at(phase, i, jp, vapor_phase) -
                       phase_fraction_at(phase, i, jm, vapor_phase)) / gy_denom;
    const double norm = std::hypot(gx, gy);
    if (norm < 1.0e-14 || flux.peak_driving == 0.0) {
        return 0.0;
    }

    // grad(phi_v) points from film into vapor.  The incident flux vector
    // points towards the film, hence exposed surfaces have -Jhat.n > 0.
    return std::max(-(flux.ux * gx + flux.uy * gy) / norm, 0.0);
}

double solid_partition_denominator(const PhaseField& phase,
                                   const decltype(phase.Fields(0, 0, 0))& node,
                                   std::size_t solid_phase,
                                   double power) {
    double denominator = 0.0;
    for (auto field = node.cbegin(); field != node.cend(); ++field) {
        if (phase.FieldsProperties[field->index].Phase == solid_phase) {
            denominator += std::pow(std::max(field->value, 0.0), power);
        }
    }
    return denominator;
}

py::array_t<double> evaluate_deposition_impl(
        const PhaseField& phase, std::size_t vapor_phase,
        double flux_x, double flux_y, double flux_to_driving,
        double partition_power) {
    const auto flux = validate_flux(flux_x, flux_y, flux_to_driving,
                                    partition_power);
    const auto [ny, nx] = planar_shape(phase);
    py::array_t<double> result({ny, nx});
    auto out = result.mutable_unchecked<2>();
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            out(j, i) = flux.peak_driving *
                        exposure_at(phase, i, j, vapor_phase, flux);
        }
    }
    return result;
}

}  // namespace

PYBIND11_MODULE(openphase_deposition, module) {
    module.doc() = "Directional PVD driving-force coupling for OpenPhase";

    module.def(
        "add_deposition_driving_force",
        [](DrivingForce& force, const PhaseField& phase,
           std::size_t vapor_phase, std::size_t solid_phase,
           double flux_x, double flux_y, double flux_to_driving,
           double partition_power) {
            if (vapor_phase == solid_phase) {
                throw std::invalid_argument(
                    "vapor_phase and solid_phase must be different");
            }
            const auto flux = validate_flux(flux_x, flux_y, flux_to_driving,
                                            partition_power);
            const auto [ny, nx] = planar_shape(phase);
            py::array_t<double> result({ny, nx});
            auto out = result.mutable_unchecked<2>();

            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    out(j, i) = 0.0;
                    const auto& node = phase.Fields(i, j, 0);
                    if (!node.interface()) {
                        continue;
                    }
                    const double exposure = exposure_at(
                        phase, i, j, vapor_phase, flux);
                    const double local_drive = flux.peak_driving * exposure;
                    out(j, i) = local_drive;
                    if (local_drive == 0.0) {
                        continue;
                    }

                    const double denominator = solid_partition_denominator(
                        phase, node, solid_phase, flux.partition_power);
                    if (denominator <= 1.0e-30) {
                        continue;
                    }

                    for (auto a = node.cbegin(); a != node.cend() - 1; ++a) {
                        for (auto b = a + 1; b != node.cend(); ++b) {
                            const auto phase_a =
                                phase.FieldsProperties[a->index].Phase;
                            const auto phase_b =
                                phase.FieldsProperties[b->index].Phase;
                            if (!((phase_a == vapor_phase && phase_b == solid_phase) ||
                                  (phase_a == solid_phase && phase_b == vapor_phase))) {
                                continue;
                            }
                            const double solid_value =
                                phase_a == solid_phase ? a->value : b->value;
                            const double share = std::pow(
                                std::max(solid_value, 0.0), flux.partition_power) /
                                denominator;
                            const double dg = local_drive * share;

                            // OpenPhase add_raw(A,B,dG) increases A relative
                            // to B for positive dG.  Deposition must therefore
                            // favor the solid member of the vapor-solid pair.
                            if (phase_a == solid_phase) {
                                force.Force(i, j, 0).add_raw(a->index, b->index, dg);
                            } else {
                                force.Force(i, j, 0).add_raw(a->index, b->index, -dg);
                            }
                        }
                    }
                }
            }
            return result;
        },
        py::arg("force"), py::arg("phase"), py::arg("vapor_phase"),
        py::arg("solid_phase"), py::arg("flux_x"), py::arg("flux_y"),
        py::arg("flux_to_driving"), py::arg("partition_power") = 4.0);

    module.def(
        "evaluate_deposition",
        &evaluate_deposition_impl,
        py::arg("phase"), py::arg("vapor_phase"), py::arg("flux_x"),
        py::arg("flux_y"), py::arg("flux_to_driving"),
        py::arg("partition_power") = 4.0);

    module.def("phase_fraction",
        [](const PhaseField& phase, std::size_t phase_index) {
            const auto [ny, nx] = planar_shape(phase);
            py::array_t<double> result({ny, nx});
            auto out = result.mutable_unchecked<2>();
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    out(j, i) = phase_fraction_at(phase, i, j, phase_index);
                }
            }
            return result;
        }, py::arg("phase"), py::arg("phase_index"));

    module.def("dominant_grain",
        [](const PhaseField& phase, std::size_t solid_phase) {
            const auto [ny, nx] = planar_shape(phase);
            py::array_t<int> result({ny, nx});
            auto out = result.mutable_unchecked<2>();
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    int best = -1;
                    double best_value = 0.0;
                    const auto& node = phase.Fields(i, j, 0);
                    for (auto field = node.cbegin(); field != node.cend(); ++field) {
                        if (phase.FieldsProperties[field->index].Phase == solid_phase &&
                            field->value > best_value) {
                            best = static_cast<int>(field->index);
                            best_value = field->value;
                        }
                    }
                    out(j, i) = best;
                }
            }
            return result;
        }, py::arg("phase"), py::arg("solid_phase"));

    module.def("maximum_partition_error",
        [](const PhaseField& phase) {
            const auto [ny, nx] = planar_shape(phase);
            double maximum = 0.0;
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    double sum = 0.0;
                    const auto& node = phase.Fields(i, j, 0);
                    for (auto field = node.cbegin(); field != node.cend(); ++field) {
                        sum += field->value;
                    }
                    maximum = std::max(maximum, std::abs(sum - 1.0));
                }
            }
            return maximum;
        }, py::arg("phase"));
}
