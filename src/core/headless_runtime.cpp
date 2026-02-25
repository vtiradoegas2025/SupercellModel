/**
 * @file headless_runtime.cpp
 * @brief Core runtime implementation for the tornado model.
 *
 * Provides simulation orchestration and subsystem integration
 * for dynamics, numerics, physics, and runtime execution paths.
 * This file belongs to the primary src/core execution layer.
 */

#include "headless_runtime.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "advection.hpp"
#include "field_contract.hpp"
#include "field_validation.hpp"
#include "runtime_config.hpp"
#include "simulation.hpp"

/**
 * @brief Executes the simulation loop in headless mode with periodic exports.
 * @param options Headless runtime options from CLI/config integration.
 * @return Zero on success, non-zero on runtime or export failure.
 */
int run_headless_simulation(const HeadlessRunOptions& options)
{
    const int export_ms = options.export_ms;
    const int duration_s = options.duration_s;
    const int write_every_s = options.write_every_s;
    const std::string& outdir = options.outdir;
    const int thetaIndex = 0;
    const bool verbose_export_debug = log_debug_enabled() || (std::getenv("TORNADO_DEBUG_EXPORTS") != nullptr);
    std::vector<tmv::ValidationReport> export_validation_reports;
    std::string validation_error_message;

    auto validate_core_fields = [&](const std::string& context,
                                    int step_index,
                                    bool include_percentiles,
                                    const std::filesystem::path* report_path,
                                    bool persist_report) -> bool
    {
        tmv::ValidationReport report;
        report.context = context;
        report.step_index = step_index;
        report.guard_mode = tmv::to_string(global_validation_policy.mode);
        report.guard_fail_on = tmv::to_string(global_validation_policy.fail_on);
        report.guard_scope = tmv::to_string(global_validation_policy.strict_scope);

        const std::vector<std::pair<std::string, Field3D*>> fields = 
        {
            {"u", &u},
            {"v", &v_theta},
            {"w", &w},
            {"rho", &rho},
            {"p", &p},
            {"theta", &theta},
            {"qv", &qv},
            {"qc", &qc},
            {"qr", &qr},
            {"qh", &qh},
            {"qg", &qg},
            {"radar", &radar_reflectivity},
            {"tracer", &tracer},
            {"qi", &qi},
            {"qs", &qs},
            {"vorticity_r", &vorticity_r},
            {"vorticity_theta", &vorticity_theta},
            {"vorticity_z", &vorticity_z},
            {"stretching_term", &stretching_term},
            {"tilting_term", &tilting_term},
            {"baroclinic_term", &baroclinic_term},
            {"p_prime", &p_prime},
            {"dynamic_pressure", &dynamic_pressure},
            {"buoyancy_pressure", &buoyancy_pressure},
            {"angular_momentum", &angular_momentum},
            {"angular_momentum_tendency", &angular_momentum_tendency},
        };

        for (const auto& binding : fields)
        {
            const tmv::FieldContract* contract = tmv::find_field_contract(binding.first);
            if (contract == nullptr)
            {
                report.missing_required.push_back(binding.first);
                continue;
            }
            if (contract->status == tmv::FieldImplementationStatus::NotImplemented)
            {
                if (contract->requirement == tmv::FieldRequirementTier::RequiredNow)
                {
                    report.missing_required.push_back(binding.first);
                }
                continue;
            }

            tmv::FieldValidationReport field_report;
            field_report.field_id = contract->id;
            field_report.status = tmv::to_string(contract->status);
            field_report.requirement = tmv::to_string(contract->requirement);
            field_report.result = tmv::validate_field3d_inplace(
                *binding.second,
                *contract,
                global_validation_policy,
                include_percentiles);

            if (field_report.result.failed)
            {
                report.failed = true;
            }
            report.fields.push_back(std::move(field_report));
        }

        if (persist_report)
        {
            const auto not_implemented = tmv::contracts_with_status(tmv::FieldImplementationStatus::NotImplemented);
            report.known_not_implemented.reserve(not_implemented.size());
            for (const tmv::FieldContract* contract : not_implemented)
            {
                report.known_not_implemented.push_back(contract->id);
                if (contract->requirement == tmv::FieldRequirementTier::RequiredNow)
                {
                    report.missing_not_implemented.push_back(contract->id);
                }
            }
        }

        if (!report.missing_required.empty() && global_validation_policy.mode == tmv::GuardMode::Strict)
        {
            report.failed = true;
        }

        if (report_path != nullptr)
        {
            std::string write_error;
            if (!tmv::write_validation_report_json(report, *report_path, write_error))
            {
                validation_error_message = write_error;
                return false;
            }
        }

        if (persist_report)
        {
            export_validation_reports.push_back(report);
        }

        if (report.failed && global_validation_policy.mode == tmv::GuardMode::Strict)
        {
            validation_error_message = "strict guard failure in context='" + context +
                "' step=" + std::to_string(step_index);
            std::cerr << "[VALIDATION] strict guard failure in context='" << context
                      << "' step=" << step_index << std::endl;
            return false;
        }

        return true;
    };

    auto write_validation_summary = [&](const std::filesystem::path& summary_path) -> bool
    {
        tmv::ValidationReport summary;
        summary.context = "run_summary";
        summary.step_index = static_cast<int>(export_validation_reports.size());
        summary.guard_mode = tmv::to_string(global_validation_policy.mode);
        summary.guard_fail_on = tmv::to_string(global_validation_policy.fail_on);
        summary.guard_scope = tmv::to_string(global_validation_policy.strict_scope);
        summary.failed = false;

        std::size_t failed_reports = 0;
        std::size_t total_violations = 0;
        std::size_t total_nonfinite = 0;
        std::size_t total_bounds = 0;
        std::size_t total_sanitized_nonfinite = 0;
        std::size_t total_sanitized_bounds = 0;

        for (const auto& report : export_validation_reports)
        {
            if (report.failed)
            {
                ++failed_reports;
                summary.failed = true;
            }

            for (const auto& field : report.fields)
            {
                total_nonfinite += field.result.stats.nan_count + field.result.stats.inf_count;
                total_bounds += field.result.stats.below_min_count + field.result.stats.above_max_count;
                total_sanitized_nonfinite += field.result.stats.sanitized_nonfinite_count;
                total_sanitized_bounds += field.result.stats.sanitized_bounds_count;
                total_violations += field.result.violations.size();
            }
        }

        std::error_code ec;
        std::filesystem::create_directories(summary_path.parent_path(), ec);
        if (ec)
        {
            validation_error_message = "failed to create summary directory '" +
                summary_path.parent_path().string() + "': " + ec.message();
            return false;
        }

        std::ofstream out(summary_path);
        if (!out)
        {
            validation_error_message = "failed to open summary report: " + summary_path.string();
            return false;
        }

        out << "{\n";
        out << "  \"context\": \"run_summary\",\n";
        out << "  \"guard_mode\": \"" << json_escape_local(summary.guard_mode) << "\",\n";
        out << "  \"guard_fail_on\": \"" << json_escape_local(summary.guard_fail_on) << "\",\n";
        out << "  \"guard_scope\": \"" << json_escape_local(summary.guard_scope) << "\",\n";
        out << "  \"report_count\": " << export_validation_reports.size() << ",\n";
        out << "  \"failed_report_count\": " << failed_reports << ",\n";
        out << "  \"total_violations\": " << total_violations << ",\n";
        out << "  \"total_nonfinite\": " << total_nonfinite << ",\n";
        out << "  \"total_bounds\": " << total_bounds << ",\n";
        out << "  \"total_sanitized_nonfinite\": " << total_sanitized_nonfinite << ",\n";
        out << "  \"total_sanitized_bounds\": " << total_sanitized_bounds << ",\n";
        out << "  \"failed\": " << (summary.failed ? "true" : "false") << "\n";
        out << "}\n";

        if (!out.good())
        {
            validation_error_message = "failed to write summary report: " + summary_path.string();
            return false;
        }

        return true;
    };
#ifdef EXPORT_NPY
    auto write_npy_2d = [&](const std::vector<float>& buf, const std::string& filename)
    {
        const std::size_t expected_size = static_cast<std::size_t>(NR) * static_cast<std::size_t>(NZ);
        if (buf.size() != expected_size)
        {
            return;
        }

        std::string header_dict = "{'descr': '<f4', 'fortran_order': False, 'shape': (" +
            std::to_string(NZ) + ", " + std::to_string(NR) + "), }";
        std::size_t header_len = header_dict.size() + 1;

        const std::size_t preamble = 6 + 2 + 2;
        const std::size_t total = preamble + header_len;
        const std::size_t padding = (16 - (total % 16)) % 16;

        header_len += padding;
        std::ofstream out(filename, std::ios::binary);

        if (!out){return;}

        out.write("\x93NUMPY", 6);
        out.put(static_cast<char>(1));
        out.put(static_cast<char>(0));

        const uint16_t hl = static_cast<uint16_t>(header_len);
        char lenb[2];

        lenb[0] = static_cast<char>(hl & 0xFF);
        lenb[1] = static_cast<char>((hl >> 8) & 0xFF);

        out.write(lenb, 2);
        out.write(header_dict.c_str(), static_cast<std::streamsize>(header_dict.size()));

        for (std::size_t i = 0; i < header_len - (header_dict.size() + 1); ++i)
        {
            out.put(' ');
        }
        out.put('\n');
        out.write(reinterpret_cast<const char*>(buf.data()),
                  static_cast<std::streamsize>(buf.size() * sizeof(float)));
        out.close();
    };

    auto save_field_slice_npy = [&](const Field3D& field, int theta, const std::string& filename)
    {
        std::vector<float> buf;
        buf.resize(static_cast<size_t>(NR) * static_cast<size_t>(NZ));
        size_t idx = 0;
        
        static bool debug_save = true;
        if (verbose_export_debug && debug_save && theta == 0) 
        {
            float field_min = 1e10, field_max = -1e10;
            int nan_count = 0;

            for (int i = 0; i < NR && i < 10; ++i) 
            {
                for (int k = 0; k < NZ && k < 10; ++k) 
                {
                    float val = field[i][theta][k];
                    if (std::isnan(val)) nan_count++;
                    if (val < field_min) field_min = val;
                    if (val > field_max) field_max = val;
                }
            }
            std::cout << "[SAVE DEBUG] Saving field slice theta=" << theta << ":" << std::endl;
            std::cout << "  Sample range: [" << field_min << ", " << field_max << "]" << std::endl;
            std::cout << "  NaN count (sample): " << nan_count << std::endl;
            debug_save = false;
        }
        
        for (int k = 0; k < NZ; ++k)
        {
            for (int i = 0; i < NR; ++i)
            {
                float val = static_cast<float>(field[i][theta][k]);
                buf[idx++] = val;
            }
        }
        write_npy_2d(buf, filename);
    };

    auto write_all_fields = [&](int export_index) -> bool
    {
        if (verbose_export_debug && export_index == 0) 
        {
            std::cout << "\n[EXPORT DEBUG] Writing timestep " << export_index << " (t=" << simulation_time << "s)" << std::endl;

            float theta_min = 1e10, theta_max = -1e10;
            float u_min = 1e10, u_max = -1e10;
            int nan_count = 0;

            for (int i = 0; i < NR && i < 10; ++i) 
            {
                for (int j = 0; j < NTH && j < 3; ++j) 
                {
                    for (int k = 0; k < NZ && k < 3; ++k) 
                    {
                        if (std::isnan(theta[i][j][k])) nan_count++;
                        if (theta[i][j][k] < theta_min) theta_min = theta[i][j][k];
                        if (theta[i][j][k] > theta_max) theta_max = theta[i][j][k];
                        if (u[i][j][k] < u_min) u_min = u[i][j][k];
                        if (u[i][j][k] > u_max) u_max = u[i][j][k];
                    }
                }
            }
            std::cout << "  Theta sample: [" << theta_min << ", " << theta_max << "] K" << std::endl;
            std::cout << "  Wind (u) sample: [" << u_min << ", " << u_max << "] m/s" << std::endl;
            std::cout << "  NaN count (sample): " << nan_count << std::endl;

            if (theta_min < 0 || theta_max > 500) 
            {
                std::cerr << "  ⚠️  ERROR: Theta values are wrong before export!" << std::endl;
            }
            std::cout << std::endl;
        }
        
        std::ostringstream stepdir;
        stepdir << outdir << "/step_" << std::setfill('0') << std::setw(6) << export_index;
        const std::filesystem::path step_path(stepdir.str());
        std::error_code remove_ec;
        std::filesystem::remove_all(step_path, remove_ec);
        std::filesystem::create_directories(step_path);

        std::filesystem::path report_path = step_path / "validation_report.json";
        if (!global_validation_report_path.empty())
        {
            std::filesystem::path report_root(global_validation_report_path);
            std::ostringstream report_name;
            report_name << "step_" << std::setfill('0') << std::setw(6) << export_index
                        << "_validation_report.json";
            report_path = report_root / report_name.str();
        }

        if (!validate_core_fields("pre_export", export_index, true, &report_path, true))
        {
            return false;
        }

        const std::size_t slice_size = static_cast<std::size_t>(NR) * static_cast<std::size_t>(NZ);
        float reflectivity_dbz_min = -30.0f;
        float reflectivity_dbz_max = 120.0f;

        if (const tmv::FieldContract* reflectivity_contract = tmv::find_field_contract("reflectivity_dbz"))
        {
            if (reflectivity_contract->default_bounds.has_min)
            {
                reflectivity_dbz_min =
                    static_cast<float>(reflectivity_contract->default_bounds.min_value);
            }
            if (reflectivity_contract->default_bounds.has_max)
            {
                reflectivity_dbz_max =
                    static_cast<float>(reflectivity_contract->default_bounds.max_value);
            }
        }
        std::vector<float> theta_azimuth_mean(slice_size, 0.0f);
        for (int i = 0; i < NR; ++i)
        {
            for (int k = 0; k < NZ; ++k)
            {
                double sum = 0.0;
                for (int j = 0; j < NTH; ++j)
                {
                    sum += static_cast<double>(theta[i][j][k]);
                }
                theta_azimuth_mean[static_cast<std::size_t>(i) * static_cast<std::size_t>(NZ) + static_cast<std::size_t>(k)] =
                    static_cast<float>(sum / static_cast<double>(NTH));
            }
        }

        std::vector<float> temperature_slice(slice_size, 0.0f);
        std::vector<float> theta_prime_slice(slice_size, 0.0f);
        std::vector<float> theta_v_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> theta_e_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> dewpoint_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> relative_humidity_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> saturation_mixing_ratio_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> total_condensate_slice(slice_size, 0.0f);
        std::vector<float> reflectivity_dbz_slice(slice_size, -30.0f);
        std::vector<float> vorticity_magnitude_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> divergence_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> buoyancy_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> horizontal_vorticity_streamwise_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> horizontal_vorticity_crosswise_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> pressure_gradient_force_x_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> pressure_gradient_force_y_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> pressure_gradient_force_z_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> storm_relative_winds_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> helicity_density_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> okubo_weiss_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> theta_w_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> zdr_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> kdp_slice(slice_size, std::numeric_limits<float>::quiet_NaN());
        std::vector<float> rhohv_slice(slice_size, std::numeric_limits<float>::quiet_NaN());

        const float kappa = static_cast<float>(R_d / cp);
        const float p0f = static_cast<float>(p0);
        const float t_freezing_k = 273.15f;
        const float latent_heat_v = 2.5e6f;
        const float cp_f = static_cast<float>(cp);
        const float g_f = static_cast<float>(g);

        std::vector<float> storm_u_level_mean(static_cast<std::size_t>(NZ), 0.0f);
        std::vector<float> storm_v_level_mean(static_cast<std::size_t>(NZ), 0.0f);
        for (int k = 0; k < NZ; ++k)
        {
            double u_sum = 0.0;
            double v_sum = 0.0;
            std::size_t count = 0;

            for (int th = 0; th < NTH; ++th)
            {
                for (int i = 0; i < NR; ++i)
                {
                    const float u_sample = static_cast<float>(u[i][th][k]);
                    const float v_sample = static_cast<float>(v_theta[i][th][k]);
                    if (std::isfinite(u_sample) && std::isfinite(v_sample))
                    {
                        u_sum += static_cast<double>(u_sample);
                        v_sum += static_cast<double>(v_sample);
                        ++count;
                    }
                }
            }
            if (count > 0)
            {
                const double denom = static_cast<double>(count);
                storm_u_level_mean[static_cast<std::size_t>(k)] =
                    static_cast<float>(u_sum / denom);
                storm_v_level_mean[static_cast<std::size_t>(k)] =
                    static_cast<float>(v_sum / denom);
            }
            else
            {
                storm_u_level_mean[static_cast<std::size_t>(k)] = std::numeric_limits<float>::quiet_NaN();
                storm_v_level_mean[static_cast<std::size_t>(k)] = std::numeric_limits<float>::quiet_NaN();
            }
        }

        auto validate_derived_export_slice = [&](const char* field_id, std::vector<float>& values, int theta_index) -> bool
        {
            const tmv::FieldContract* contract = tmv::find_field_contract(field_id);
            if (contract == nullptr || contract->status != tmv::FieldImplementationStatus::ExportedNow)
            {
                return true;
            }

            tmv::FieldValidationResult result = tmv::validate_buffer_inplace(
                values.data(),
                values.size(),
                *contract,
                global_validation_policy,
                false);

            if (result.failed && global_validation_policy.mode == tmv::GuardMode::Strict)
            {
                validation_error_message =
                    "strict guard failure for derived export field='" +
                    std::string(field_id) +
                    "' step=" + std::to_string(export_index) +
                    " theta=" + std::to_string(theta_index);
                std::cerr << "[VALIDATION] strict guard failure for derived export field='"
                          << field_id << "' step=" << export_index
                          << " theta=" << theta_index << std::endl;
                return false;
            }

            return true;
        };

        const std::pair<const char*, const char*> exported_core_fields[] = {
            {"u", "u"},
            {"v", "v"},
            {"w", "w"},
            {"rho", "rho"},
            {"p", "p"},
            {"theta", "theta"},
            {"qv", "qv"},
            {"qc", "qc"},
            {"qr", "qr"},
            {"qi", "qi"},
            {"qs", "qs"},
            {"qh", "qh"},
            {"qg", "qg"},
            {"radar", "radar"},
            {"tracer", "tracer"},
            {"vorticity_r", "vorticity_r"},
            {"vorticity_theta", "vorticity_theta"},
            {"vorticity_z", "vorticity_z"},
            {"stretching_term", "stretching_term"},
            {"tilting_term", "tilting_term"},
            {"baroclinic_term", "baroclinic_term"},
            {"p_prime", "p_prime"},
            {"dynamic_pressure", "dynamic_pressure"},
            {"buoyancy_pressure", "buoyancy_pressure"},
            {"angular_momentum", "angular_momentum"},
            {"angular_momentum_tendency", "angular_momentum_tendency"},
        };

        const std::pair<const char*, const char*> exported_derived_fields[] = 
        {
            {"temperature", "temperature"},
            {"theta_prime", "theta_prime"},
            {"theta_v", "theta_v"},
            {"theta_e", "theta_e"},
            {"dewpoint", "dewpoint"},
            {"relative_humidity", "relative_humidity"},
            {"saturation_mixing_ratio", "saturation_mixing_ratio"},
            {"total_condensate", "total_condensate"},
            {"reflectivity_dbz", "reflectivity_dbz"},
            {"vorticity_magnitude", "vorticity_magnitude"},
            {"divergence", "divergence"},
            {"buoyancy", "buoyancy"},
            {"horizontal_vorticity_streamwise", "horizontal_vorticity_streamwise"},
            {"horizontal_vorticity_crosswise", "horizontal_vorticity_crosswise"},
            {"pressure_gradient_force_x", "pressure_gradient_force_x"},
            {"pressure_gradient_force_y", "pressure_gradient_force_y"},
            {"pressure_gradient_force_z", "pressure_gradient_force_z"},
            {"storm_relative_winds", "storm_relative_winds"},
            {"helicity_density", "helicity_density"},
            {"okubo_weiss", "okubo_weiss"},
            {"theta_w", "theta_w"},
            {"zdr", "zdr"},
            {"kdp", "kdp"},
            {"rhohv", "rhohv"},
        };

        auto write_step_manifest = [&](const std::filesystem::path& manifest_path) -> bool
        {
            std::ofstream out(manifest_path);
            if (!out)
            {
                return false;
            }

            out << "{\n";
            out << "  \"step_index\": " << export_index << ",\n";
            out << "  \"simulation_time_s\": " << simulation_time << ",\n";
            out << "  \"step_dir\": \"" << json_escape_local(step_path.filename().string()) << "\",\n";
            out << "  \"grid\": {\n";
            out << "    \"nr\": " << NR << ",\n";
            out << "    \"nth\": " << NTH << ",\n";
            out << "    \"nz\": " << NZ << ",\n";
            out << "    \"dr_m\": " << dr << ",\n";
            out << "    \"dtheta_rad\": " << dtheta << ",\n";
            out << "    \"dz_m\": " << dz << "\n";
            out << "  },\n";
            out << "  \"theta_index\": {\n";
            out << "    \"min\": 0,\n";
            out << "    \"max\": " << (NTH > 0 ? NTH - 1 : 0) << ",\n";
            out << "    \"count\": " << NTH << ",\n";
            out << "    \"file_prefix\": \"th{theta}_\"\n";
            out << "  },\n";
            out << "  \"soundings\": {\n";
            out << "    \"enabled\": " << (global_sounding_enabled ? "true" : "false") << ",\n";
            out << "    \"scheme\": \"" << json_escape_local(global_runtime_sounding_config.scheme_id) << "\",\n";
            out << "    \"file_path\": \"" << json_escape_local(global_runtime_sounding_config.file_path) << "\",\n";
            out << "    \"interpolation_method\": \""
                << json_escape_local(sounding_interpolation_method_name(global_runtime_sounding_config.interpolation_method))
                << "\"\n";
            out << "  },\n";
            out << "  \"fields\": [\n";

            bool first_field = true;

            auto write_field_entry = [&](const char* field_id,const char* suffix, const char* category, bool derived)
            {
                if (!first_field)
                {
                    out << ",\n";
                }
                first_field = false;

                out << "    {\n";
                out << "      \"field_id\": \"" << json_escape_local(field_id) << "\",\n";
                out << "      \"category\": \"" << json_escape_local(category) << "\",\n";
                out << "      \"derived\": " << (derived ? "true" : "false") << ",\n";
                out << "      \"file_suffix\": \"" << json_escape_local(suffix) << "\",\n";
                out << "      \"file_pattern\": \"th{theta}_" << json_escape_local(suffix) << ".npy\"";

                if (const tmv::FieldContract* contract = tmv::find_field_contract(field_id))
                {
                    out << ",\n";
                    out << "      \"units\": \"" << json_escape_local(contract->units) << "\",\n";
                    out << "      \"description\": \"" << json_escape_local(contract->description) << "\"";

                    if (contract->default_bounds.has_min || contract->default_bounds.has_max)
                    {
                        out << ",\n";
                        out << "      \"bounds\": {";
                        bool first_bound = true;
                        if (contract->default_bounds.has_min)
                        {
                            out << "\"min\": " << contract->default_bounds.min_value;
                            first_bound = false;
                        }
                        if (contract->default_bounds.has_max)
                        {
                            if (!first_bound)
                            {
                                out << ", ";
                            }
                            out << "\"max\": " << contract->default_bounds.max_value;
                        }
                        out << "}";
                    }
                }

                out << "\n";
                out << "    }";
            };

            for (const auto& field : exported_core_fields)
            {
                write_field_entry(field.first, field.second, "core", false);
            }
            for (const auto& field : exported_derived_fields)
            {
                write_field_entry(field.first, field.second, "derived", true);
            }

            out << "\n";
            out << "  ]\n";
            out << "}\n";
            return true;
        };

        for (int th = 0; th < NTH; ++th)
        {
            std::string base_path = stepdir.str() + "/th" + std::to_string(th);

            save_field_slice_npy(u, th, base_path + "_u.npy");
            save_field_slice_npy(v_theta, th, base_path + "_v.npy");
            save_field_slice_npy(w, th, base_path + "_w.npy");
            save_field_slice_npy(rho, th, base_path + "_rho.npy");
            save_field_slice_npy(p, th, base_path + "_p.npy");
            save_field_slice_npy(theta, th, base_path + "_theta.npy");
            save_field_slice_npy(qv, th, base_path + "_qv.npy");
            save_field_slice_npy(qc, th, base_path + "_qc.npy");
            save_field_slice_npy(qr, th, base_path + "_qr.npy");
            save_field_slice_npy(qi, th, base_path + "_qi.npy");
            save_field_slice_npy(qs, th, base_path + "_qs.npy");
            save_field_slice_npy(qh, th, base_path + "_qh.npy");
            save_field_slice_npy(qg, th, base_path + "_qg.npy");
            save_field_slice_npy(radar_reflectivity, th, base_path + "_radar.npy");
            save_field_slice_npy(tracer, th, base_path + "_tracer.npy");
            save_field_slice_npy(vorticity_r, th, base_path + "_vorticity_r.npy");
            save_field_slice_npy(vorticity_theta, th, base_path + "_vorticity_theta.npy");
            save_field_slice_npy(vorticity_z, th, base_path + "_vorticity_z.npy");
            save_field_slice_npy(stretching_term, th, base_path + "_stretching_term.npy");
            save_field_slice_npy(tilting_term, th, base_path + "_tilting_term.npy");
            save_field_slice_npy(baroclinic_term, th, base_path + "_baroclinic_term.npy");
            save_field_slice_npy(p_prime, th, base_path + "_p_prime.npy");
            save_field_slice_npy(dynamic_pressure, th, base_path + "_dynamic_pressure.npy");
            save_field_slice_npy(buoyancy_pressure, th, base_path + "_buoyancy_pressure.npy");
            save_field_slice_npy(angular_momentum, th, base_path + "_angular_momentum.npy");
            save_field_slice_npy(angular_momentum_tendency, th, base_path + "_angular_momentum_tendency.npy");

            std::size_t idx = 0;
            for (int k = 0; k < NZ; ++k)
            {
                for (int i = 0; i < NR; ++i)
                {
                    const float p_pa = static_cast<float>(p[i][th][k]);
                    const float theta_value = static_cast<float>(theta[i][th][k]);
                    const float qv_value = static_cast<float>(qv[i][th][k]);
                    if (std::isfinite(p_pa) && std::isfinite(theta_value) && p_pa > 0.0f)
                    {
                        temperature_slice[idx] = theta_value * std::pow(p_pa / p0f, kappa);
                    }
                    else
                    {
                        temperature_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                    }

                    const float temperature_k = temperature_slice[idx];
                    if (std::isfinite(temperature_k) && std::isfinite(p_pa) && p_pa > 0.0f)
                    {
                        const float p_hpa = p_pa * 0.01f;
                        float es_hpa = 6.112f * std::exp(
                            17.67f * (temperature_k - t_freezing_k) / (temperature_k - 29.65f));
                        if (std::isfinite(es_hpa) && p_hpa > 1.0e-3f)
                        {
                            es_hpa = std::max(0.0f, std::min(es_hpa, 0.99f * p_hpa));
                            const float denom = std::max(1.0e-3f, p_hpa - es_hpa);
                            float qsat = 0.622f * es_hpa / denom;

                            qsat = std::max(1.0e-8f, std::min(0.10f, qsat));
                            saturation_mixing_ratio_slice[idx] = qsat;

                            if (std::isfinite(qv_value) && std::isfinite(qsat) && qsat > 0.0f)
                            {
                                const float rh = std::max(0.0f, std::min(200.0f, 100.0f * (qv_value / qsat)));
                                relative_humidity_slice[idx] = rh;
                                const float rh_frac = std::max(1.0e-6f, rh * 0.01f);
                                const float temperature_c = temperature_k - t_freezing_k;
                                const float gamma = std::log(rh_frac) +
                                    (17.625f * temperature_c) / (243.04f + temperature_c);
                                const float gamma_denom = 17.625f - gamma;
                                if (std::abs(gamma_denom) > 1.0e-6f)
                                {
                                    const float td_c = 243.04f * gamma / gamma_denom;
                                    dewpoint_slice[idx] = td_c + t_freezing_k;
                                }
                            }
                        }
                    }

                    const std::size_t mean_idx = static_cast<std::size_t>(i) * static_cast<std::size_t>(NZ) + static_cast<std::size_t>(k);
                    theta_prime_slice[idx] = theta_value - theta_azimuth_mean[mean_idx];

                    const float qc_value = std::max(0.0f, static_cast<float>(qc[i][th][k]));
                    const float qr_value = std::max(0.0f, static_cast<float>(qr[i][th][k]));
                    const float qi_value = std::max(0.0f, static_cast<float>(qi[i][th][k]));
                    const float qs_value = std::max(0.0f, static_cast<float>(qs[i][th][k]));
                    const float qg_value = std::max(0.0f, static_cast<float>(qg[i][th][k]));
                    const float qh_value = std::max(0.0f, static_cast<float>(qh[i][th][k]));

                    total_condensate_slice[idx] = qc_value + qr_value + qi_value + qs_value + qg_value + qh_value;

                    const float qv_nonnegative = std::isfinite(qv_value) ? std::max(0.0f, qv_value) : std::numeric_limits<float>::quiet_NaN();

                    if (std::isfinite(theta_value) &&std::isfinite(qv_nonnegative) && std::isfinite(total_condensate_slice[idx]))
                    {
                        theta_v_slice[idx] = theta_value * (1.0f + (0.61f * qv_nonnegative) - total_condensate_slice[idx]);
                    }
                    else
                    {
                        theta_v_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                    }

                    if (std::isfinite(theta_value) &&
                        std::isfinite(qv_nonnegative) &&
                        std::isfinite(temperature_k) &&
                        temperature_k > 0.0f)
                    {
                        const float exponent = (latent_heat_v * qv_nonnegative) / (cp_f * temperature_k);
                        theta_e_slice[idx] = theta_value * std::exp(exponent);
                    }
                    else
                    {
                        theta_e_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                    }

                    const float z_linear = static_cast<float>(radar_reflectivity[i][th][k]);
                    if (!std::isfinite(z_linear))
                    {
                        reflectivity_dbz_slice[idx] = reflectivity_dbz_min;
                    }
                    else if (z_linear <= 0.0f)
                    {
                        reflectivity_dbz_slice[idx] = reflectivity_dbz_min;
                    }
                    else
                    {
                        float dbz_value = 10.0f * std::log10(z_linear);
                        if (!std::isfinite(dbz_value))
                        {
                            dbz_value = reflectivity_dbz_min;
                        }
                        reflectivity_dbz_slice[idx] =
                            std::max(reflectivity_dbz_min, std::min(reflectivity_dbz_max, dbz_value));
                    }

                    const float vort_r = static_cast<float>(vorticity_r[i][th][k]);
                    const float vort_theta = static_cast<float>(vorticity_theta[i][th][k]);
                    const float vort_z = static_cast<float>(vorticity_z[i][th][k]);
                    if (std::isfinite(vort_r) && std::isfinite(vort_theta) && std::isfinite(vort_z))
                    {
                        vorticity_magnitude_slice[idx] = std::sqrt(
                            (vort_r * vort_r) + (vort_theta * vort_theta) + (vort_z * vort_z));
                    }
                    else
                    {
                        vorticity_magnitude_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                    }

                    const int th_plus = (th + 1) % NTH;
                    const int th_minus = (th + NTH - 1) % NTH;
                    const float u_center = static_cast<float>(u[i][th][k]);
                    const float v_center = static_cast<float>(v_theta[i][th][k]);
                    const float w_center = static_cast<float>(w[i][th][k]);
                    const float u_plus = static_cast<float>(u[(i + 1 < NR) ? i + 1 : i][th][k]);
                    const float u_minus = static_cast<float>(u[(i > 0) ? i - 1 : i][th][k]);
                    const float u_th_plus = static_cast<float>(u[i][th_plus][k]);
                    const float u_th_minus = static_cast<float>(u[i][th_minus][k]);
                    const float v_plus = static_cast<float>(v_theta[i][th_plus][k]);
                    const float v_minus = static_cast<float>(v_theta[i][th_minus][k]);
                    const float v_r_plus = static_cast<float>(v_theta[(i + 1 < NR) ? i + 1 : i][th][k]);
                    const float v_r_minus = static_cast<float>(v_theta[(i > 0) ? i - 1 : i][th][k]);
                    const float w_plus = static_cast<float>(w[i][th][(k + 1 < NZ) ? k + 1 : k]);
                    const float w_minus = static_cast<float>(w[i][th][(k > 0) ? k - 1 : k]);

                    float radial_term = std::numeric_limits<float>::quiet_NaN();
                    if (std::isfinite(u_center) && std::isfinite(u_plus) && std::isfinite(u_minus) && dr > 0.0)
                    {
                        const float dr_f = static_cast<float>(dr);
                        if (i == 0 && NR > 1)
                        {
                            radial_term = 2.0f * (u_plus - u_center) / dr_f;
                        }
                        else
                        {
                            const float r_center = static_cast<float>(i) * dr_f;
                            if (r_center > 0.0f)
                            {
                                const float r_plus = static_cast<float>(std::min(i + 1, NR - 1)) * dr_f;
                                const float r_minus = static_cast<float>((i > 0) ? i - 1 : 0) * dr_f;
                                const float denominator = (i == NR - 1) ? dr_f : (2.0f * dr_f);
                                const float d_ru_dr = ((r_plus * u_plus) - (r_minus * u_minus)) / denominator;
                                radial_term = d_ru_dr / r_center;
                            }
                        }
                    }

                    float azimuthal_term = 0.0f;
                    if (i > 0 && dtheta > 0.0 && dr > 0.0 &&
                        std::isfinite(v_plus) && std::isfinite(v_minus))
                    {
                        const float r_center = static_cast<float>(i) * static_cast<float>(dr);
                        if (r_center > 0.0f)
                        {
                            azimuthal_term =
                                (v_plus - v_minus) / (2.0f * static_cast<float>(dtheta) * r_center);
                        }
                    }

                    float vertical_term = std::numeric_limits<float>::quiet_NaN();
                    if (std::isfinite(w_plus) && std::isfinite(w_minus) && dz > 0.0)
                    {
                        const float dz_f = static_cast<float>(dz);
                        const float denominator = (k == 0 || k == NZ - 1) ? dz_f : (2.0f * dz_f);
                        vertical_term = (w_plus - w_minus) / denominator;
                    }

                    if (std::isfinite(radial_term) && std::isfinite(azimuthal_term) && std::isfinite(vertical_term))
                    {
                        divergence_slice[idx] = radial_term + azimuthal_term + vertical_term;
                    }
                    else
                    {
                        divergence_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                    }

                    const float theta_mean = theta_azimuth_mean[mean_idx];
                    if (std::isfinite(theta_prime_slice[idx]) &&
                        std::isfinite(theta_mean) &&
                        std::abs(theta_mean) > 1.0e-6f)
                    {
                        buoyancy_slice[idx] = g_f * (theta_prime_slice[idx] / theta_mean);
                    }
                    else
                    {
                        buoyancy_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                    }

                    if (std::isfinite(vort_r) &&
                        std::isfinite(vort_theta) &&
                        std::isfinite(u_center) &&
                        std::isfinite(v_center))
                    {
                        const float speed_h = std::sqrt((u_center * u_center) + (v_center * v_center));
                        if (speed_h > 1.0e-6f)
                        {
                            horizontal_vorticity_streamwise_slice[idx] =
                                ((vort_r * u_center) + (vort_theta * v_center)) / speed_h;
                            horizontal_vorticity_crosswise_slice[idx] =
                                ((-vort_r * v_center) + (vort_theta * u_center)) / speed_h;
                        }
                        else
                        {
                            horizontal_vorticity_streamwise_slice[idx] = 0.0f;
                            horizontal_vorticity_crosswise_slice[idx] = 0.0f;
                        }
                    }
                    else
                    {
                        horizontal_vorticity_streamwise_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                        horizontal_vorticity_crosswise_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                    }

                    const float storm_u_mean = storm_u_level_mean[static_cast<std::size_t>(k)];
                    const float storm_v_mean = storm_v_level_mean[static_cast<std::size_t>(k)];
                    if (std::isfinite(u_center) &&
                        std::isfinite(v_center) &&
                        std::isfinite(storm_u_mean) &&
                        std::isfinite(storm_v_mean))
                    {
                        const float u_sr = u_center - storm_u_mean;
                        const float v_sr = v_center - storm_v_mean;
                        storm_relative_winds_slice[idx] = std::sqrt((u_sr * u_sr) + (v_sr * v_sr));
                    }
                    else
                    {
                        storm_relative_winds_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                    }

                    if (std::isfinite(u_center) &&
                        std::isfinite(v_center) &&
                        std::isfinite(w_center) &&
                        std::isfinite(vort_r) &&
                        std::isfinite(vort_theta) &&
                        std::isfinite(vort_z))
                    {
                        helicity_density_slice[idx] =
                            (u_center * vort_r) +
                            (v_center * vort_theta) +
                            (w_center * vort_z);
                    }
                    else
                    {
                        helicity_density_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                    }

                    if (dr > 0.0 && dtheta > 0.0)
                    {
                        const float dr_f = static_cast<float>(dr);
                        const float radial_denom = (i == 0 || i == NR - 1) ? dr_f : (2.0f * dr_f);
                        float du_dr = std::numeric_limits<float>::quiet_NaN();
                        float dv_dr = std::numeric_limits<float>::quiet_NaN();
                        if (std::isfinite(u_plus) && std::isfinite(u_minus))
                        {
                            du_dr = (u_plus - u_minus) / radial_denom;
                        }
                        if (std::isfinite(v_r_plus) && std::isfinite(v_r_minus))
                        {
                            dv_dr = (v_r_plus - v_r_minus) / radial_denom;
                        }

                        float du_dy = 0.0f;
                        float dv_dy = 0.0f;
                        bool azimuthal_derivatives_valid = false;
                        const float r_center = static_cast<float>(i) * dr_f;
                        if (r_center > 0.0f &&
                            std::isfinite(u_th_plus) &&
                            std::isfinite(u_th_minus) &&
                            std::isfinite(v_plus) &&
                            std::isfinite(v_minus))
                        {
                            const float inv_arc = 1.0f / (2.0f * static_cast<float>(dtheta) * r_center);
                            du_dy = (u_th_plus - u_th_minus) * inv_arc;
                            dv_dy = (v_plus - v_minus) * inv_arc;
                            azimuthal_derivatives_valid = true;
                        }
                        else if (r_center <= 0.0f)
                        {
                            du_dy = 0.0f;
                            dv_dy = 0.0f;
                            azimuthal_derivatives_valid = true;
                        }

                        if (std::isfinite(du_dr) &&
                            std::isfinite(dv_dr) &&
                            azimuthal_derivatives_valid &&
                            std::isfinite(vort_z))
                        {
                            const float normal_strain = du_dr - dv_dy;
                            const float shear_strain = dv_dr + du_dy;
                            okubo_weiss_slice[idx] =
                                (normal_strain * normal_strain) +
                                (shear_strain * shear_strain) -
                                (vort_z * vort_z);
                        }
                        else
                        {
                            okubo_weiss_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                        }
                    }
                    else
                    {
                        okubo_weiss_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                    }

                    const float p_center = static_cast<float>(p[i][th][k]);
                    const float p_r_plus = static_cast<float>(p[(i + 1 < NR) ? i + 1 : i][th][k]);
                    const float p_r_minus = static_cast<float>(p[(i > 0) ? i - 1 : i][th][k]);
                    const float p_th_plus = static_cast<float>(p[i][th_plus][k]);
                    const float p_th_minus = static_cast<float>(p[i][th_minus][k]);
                    const float p_z_plus = static_cast<float>(p[i][th][(k + 1 < NZ) ? k + 1 : k]);
                    const float p_z_minus = static_cast<float>(p[i][th][(k > 0) ? k - 1 : k]);
                    const float rho_center = static_cast<float>(rho[i][th][k]);

                    if (std::isfinite(p_center) &&
                        std::isfinite(p_r_plus) &&
                        std::isfinite(p_r_minus) &&
                        std::isfinite(p_th_plus) &&
                        std::isfinite(p_th_minus) &&
                        std::isfinite(p_z_plus) &&
                        std::isfinite(p_z_minus) &&
                        std::isfinite(rho_center) &&
                        rho_center > 1.0e-6f &&
                        dr > 0.0 &&
                        dz > 0.0)
                    {
                        const float dr_f = static_cast<float>(dr);
                        const float dz_f = static_cast<float>(dz);
                        const float radial_denom = (i == 0 || i == NR - 1) ? dr_f : (2.0f * dr_f);
                        const float vertical_denom = (k == 0 || k == NZ - 1) ? dz_f : (2.0f * dz_f);
                        const float dp_dr = (p_r_plus - p_r_minus) / radial_denom;
                        const float dp_dz = (p_z_plus - p_z_minus) / vertical_denom;

                        float dp_dy = 0.0f;
                        const float r_center = static_cast<float>(i) * dr_f;
                        if (r_center > 0.0f && dtheta > 0.0)
                        {
                            dp_dy = (p_th_plus - p_th_minus) /
                                (2.0f * static_cast<float>(dtheta) * r_center);
                        }

                        pressure_gradient_force_x_slice[idx] = -dp_dr / rho_center;
                        pressure_gradient_force_y_slice[idx] = -dp_dy / rho_center;
                        pressure_gradient_force_z_slice[idx] = -dp_dz / rho_center;
                    }
                    else
                    {
                        pressure_gradient_force_x_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                        pressure_gradient_force_y_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                        pressure_gradient_force_z_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                    }

                    if (std::isfinite(temperature_k) &&
                        std::isfinite(relative_humidity_slice[idx]) &&
                        std::isfinite(p_pa) &&
                        p_pa > 0.0f)
                    {
                        const float rh_clamped = std::max(1.0e-3f, std::min(100.0f, relative_humidity_slice[idx]));
                        const float t_c = temperature_k - t_freezing_k;
                        const float tw_c =
                            (t_c * std::atan(0.151977f * std::sqrt(rh_clamped + 8.313659f))) +
                            std::atan(t_c + rh_clamped) -
                            std::atan(rh_clamped - 1.676331f) +
                            (0.00391838f * std::pow(rh_clamped, 1.5f) * std::atan(0.023101f * rh_clamped)) -
                            4.686035f;
                        const float tw_k = tw_c + t_freezing_k;
                        if (std::isfinite(tw_k) && tw_k > 0.0f)
                        {
                            theta_w_slice[idx] = tw_k * std::pow(p0f / p_pa, kappa);
                        }
                        else
                        {
                            theta_w_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                        }
                    }
                    else
                    {
                        theta_w_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                    }

                    if (std::isfinite(qr_value))
                    {
                        const float qrain_nonnegative = std::max(0.0f, qr_value);
                        if (qrain_nonnegative <= 1.0e-10f)
                        {
                            zdr_slice[idx] = 0.0f;
                        }
                        else
                        {
                            const float rain_rate_proxy = qrain_nonnegative * 3600.0f;
                            float axis_ratio = 0.95f;
                            if (rain_rate_proxy >= 10.0f)
                            {
                                axis_ratio = 0.75f;
                            }
                            else if (rain_rate_proxy >= 1.0f)
                            {
                                axis_ratio = 0.85f;
                            }
                            zdr_slice[idx] = -40.0f * std::log10(axis_ratio);
                        }

                        const float mixed_ice = qi_value + qs_value + qg_value + qh_value;
                        kdp_slice[idx] = std::max(0.0f, (1500.0f * qrain_nonnegative) + (200.0f * mixed_ice));
                        const float hydro_sum = qrain_nonnegative + mixed_ice + 1.0e-8f;
                        const float mixed_fraction = mixed_ice / hydro_sum;
                        float rhohv = 1.0f - (0.35f * mixed_fraction) -
                            (0.05f * std::min(1.0f, qrain_nonnegative / 0.005f));
                        rhohv = std::max(0.5f, std::min(1.0f, rhohv));
                        rhohv_slice[idx] = rhohv;
                    }
                    else
                    {
                        zdr_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                        kdp_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                        rhohv_slice[idx] = std::numeric_limits<float>::quiet_NaN();
                    }

                    ++idx;
                }
            }

            const std::pair<const char*, std::vector<float>*> derived_bindings[] = {
                {"temperature", &temperature_slice},
                {"theta_prime", &theta_prime_slice},
                {"theta_v", &theta_v_slice},
                {"theta_e", &theta_e_slice},
                {"dewpoint", &dewpoint_slice},
                {"relative_humidity", &relative_humidity_slice},
                {"saturation_mixing_ratio", &saturation_mixing_ratio_slice},
                {"total_condensate", &total_condensate_slice},
                {"reflectivity_dbz", &reflectivity_dbz_slice},
                {"vorticity_magnitude", &vorticity_magnitude_slice},
                {"divergence", &divergence_slice},
                {"buoyancy", &buoyancy_slice},
                {"horizontal_vorticity_streamwise", &horizontal_vorticity_streamwise_slice},
                {"horizontal_vorticity_crosswise", &horizontal_vorticity_crosswise_slice},
                {"pressure_gradient_force_x", &pressure_gradient_force_x_slice},
                {"pressure_gradient_force_y", &pressure_gradient_force_y_slice},
                {"pressure_gradient_force_z", &pressure_gradient_force_z_slice},
                {"storm_relative_winds", &storm_relative_winds_slice},
                {"helicity_density", &helicity_density_slice},
                {"okubo_weiss", &okubo_weiss_slice},
                {"theta_w", &theta_w_slice},
                {"zdr", &zdr_slice},
                {"kdp", &kdp_slice},
                {"rhohv", &rhohv_slice},
            };
            for (const auto& binding : derived_bindings)
            {
                if (!validate_derived_export_slice(binding.first, *binding.second, th))
                {
                    return false;
                }
            }

            write_npy_2d(temperature_slice, base_path + "_temperature.npy");
            write_npy_2d(theta_prime_slice, base_path + "_theta_prime.npy");
            write_npy_2d(theta_v_slice, base_path + "_theta_v.npy");
            write_npy_2d(theta_e_slice, base_path + "_theta_e.npy");
            write_npy_2d(dewpoint_slice, base_path + "_dewpoint.npy");
            write_npy_2d(relative_humidity_slice, base_path + "_relative_humidity.npy");
            write_npy_2d(saturation_mixing_ratio_slice, base_path + "_saturation_mixing_ratio.npy");
            write_npy_2d(total_condensate_slice, base_path + "_total_condensate.npy");
            write_npy_2d(reflectivity_dbz_slice, base_path + "_reflectivity_dbz.npy");
            write_npy_2d(vorticity_magnitude_slice, base_path + "_vorticity_magnitude.npy");
            write_npy_2d(divergence_slice, base_path + "_divergence.npy");
            write_npy_2d(buoyancy_slice, base_path + "_buoyancy.npy");
            write_npy_2d(horizontal_vorticity_streamwise_slice, base_path + "_horizontal_vorticity_streamwise.npy");
            write_npy_2d(horizontal_vorticity_crosswise_slice, base_path + "_horizontal_vorticity_crosswise.npy");
            write_npy_2d(pressure_gradient_force_x_slice, base_path + "_pressure_gradient_force_x.npy");
            write_npy_2d(pressure_gradient_force_y_slice, base_path + "_pressure_gradient_force_y.npy");
            write_npy_2d(pressure_gradient_force_z_slice, base_path + "_pressure_gradient_force_z.npy");
            write_npy_2d(storm_relative_winds_slice, base_path + "_storm_relative_winds.npy");
            write_npy_2d(helicity_density_slice, base_path + "_helicity_density.npy");
            write_npy_2d(okubo_weiss_slice, base_path + "_okubo_weiss.npy");
            write_npy_2d(theta_w_slice, base_path + "_theta_w.npy");
            write_npy_2d(zdr_slice, base_path + "_zdr.npy");
            write_npy_2d(kdp_slice, base_path + "_kdp.npy");
            write_npy_2d(rhohv_slice, base_path + "_rhohv.npy");
        }

        if (!write_step_manifest(step_path / "manifest.json"))
        {
            std::cerr << "[EXPORT] failed to write step manifest: "
                      << (step_path / "manifest.json") << std::endl;
            return false;
        }
        return true;
    };

#endif

    if (write_every_s > 0)
    {
        std::filesystem::create_directories(outdir);

        static const std::regex step_pattern(R"(^step_[0-9]{6}$)");
        for (const auto& entry : std::filesystem::directory_iterator(outdir))
        {
            if (!entry.is_directory())
            {
                continue;
            }
            const std::string name = entry.path().filename().string();
            if (std::regex_match(name, step_pattern))
            {
                std::error_code ec;
                std::filesystem::remove_all(entry.path(), ec);
            }
        }
    }

    auto lastGuiExport = std::chrono::steady_clock::now();
    int steps = 0;
    int export_index = 0;
    ::simulation_time = 0.0;
    double next_field_export_time_s = (write_every_s > 0) ? static_cast<double>(write_every_s) : -1.0;

    using PerfClock = std::chrono::steady_clock;
    /**
     * @brief Aggregated runtime timing totals for headless simulation profiling.
     */
    struct PerfTotals
    {
        uint64_t steps = 0;
        double simulated_s = 0.0;
        double total_step_s = 0.0;
        double initial_export_s = 0.0;
        double radiation_s = 0.0;
        double boundary_layer_s = 0.0;
        double chaos_noise_s = 0.0;
        double chaos_tendency_s = 0.0;
        double dynamics_s = 0.0;
        double export_s = 0.0;
    } perf_totals;

    auto timed_call = [&](double& accumulator, auto&& fn)
    {
        if (!global_perf_timing_enabled)
        {
            fn();
            return;
        }
        const auto t0 = PerfClock::now();
        fn();
        const auto t1 = PerfClock::now();
        accumulator += std::chrono::duration<double>(t1 - t0).count();
    };

    if (global_perf_timing_enabled)
    {
        reset_advection_perf_stats();
    }
    
    if (verbose_export_debug)
    {
        float theta_min = 1e10, theta_max = -1e10;
        float u_min = 1e10, u_max = -1e10;
        int nan_count = 0;
        for (int i = 0; i < NR && i < 10; ++i) {
            for (int j = 0; j < NTH && j < 5; ++j) {
                for (int k = 0; k < NZ && k < 5; ++k) {
                    if (std::isnan(theta[i][j][k])) nan_count++;
                    if (theta[i][j][k] < theta_min) theta_min = theta[i][j][k];
                    if (theta[i][j][k] > theta_max) theta_max = theta[i][j][k];
                    if (u[i][j][k] < u_min) u_min = u[i][j][k];
                    if (u[i][j][k] > u_max) u_max = u[i][j][k];
                }
            }
        }
        std::cout << "\n[TIME STEP DEBUG] Before time stepping (t=0):" << std::endl;
        std::cout << "  Theta sample: min=" << theta_min << "K, max=" << theta_max << "K" << std::endl;
        std::cout << "  Wind (u) sample: min=" << u_min << "m/s, max=" << u_max << "m/s" << std::endl;
        std::cout << "  NaN count (sample): " << nan_count << std::endl;
        if (theta_min < 0 || theta_max > 500) {
            std::cerr << "  ⚠️  ERROR: Theta already corrupted before time stepping!" << std::endl;
        }
        std::cout << "  Headless duration/export cadence use simulation seconds" << std::endl;
        std::cout << std::endl;
    }

#ifdef EXPORT_NPY
    if (write_every_s > 0)
    {
        bool initial_export_ok = true;
        timed_call(perf_totals.initial_export_s, [&] {
            initial_export_ok = write_all_fields(export_index++);
        });
        if (!initial_export_ok)
        {
            std::cerr << "[VALIDATION] " << validation_error_message << std::endl;
            return 2;
        }
    }
#endif
    
    while (true)
    {
        const double runtime_dt = choose_runtime_timestep();
        if (std::isfinite(runtime_dt) && runtime_dt > 0.0)
        {
            dt = runtime_dt;
        }

        const auto step_t0 = PerfClock::now();

        timed_call(perf_totals.radiation_s, [&] { step_radiation(simulation_time); });
        if (!validate_core_fields("after_radiation", steps, false, nullptr, false))
        {
            std::cerr << "[VALIDATION] " << validation_error_message << std::endl;
            return 2;
        }
        timed_call(perf_totals.boundary_layer_s, [&] { step_boundary_layer(simulation_time); });
        if (!validate_core_fields("after_boundary_layer", steps, false, nullptr, false))
        {
            std::cerr << "[VALIDATION] " << validation_error_message << std::endl;
            return 2;
        }
        timed_call(perf_totals.chaos_noise_s, [&] { step_chaos_noise(dt); });
        if (!validate_core_fields("after_chaos_noise", steps, false, nullptr, false))
        {
            std::cerr << "[VALIDATION] " << validation_error_message << std::endl;
            return 2;
        }
        timed_call(perf_totals.chaos_tendency_s, [&] { apply_chaos_tendencies(); });
        if (!validate_core_fields("after_chaos_tendency", steps, false, nullptr, false))
        {
            std::cerr << "[VALIDATION] " << validation_error_message << std::endl;
            return 2;
        }
        timed_call(perf_totals.dynamics_s, [&] { step_dynamics(simulation_time); });
        if (!validate_core_fields("after_dynamics", steps, false, nullptr, false))
        {
            std::cerr << "[VALIDATION] " << validation_error_message << std::endl;
            return 2;
        }
        simulation_time += dt;
        perf_totals.simulated_s += dt;
        
        if (verbose_export_debug && steps % 100 == 0 && steps < 1000) 
        {
            float theta_min = 1e10, theta_max = -1e10;
            float u_min = 1e10, u_max = -1e10;
            int nan_count = 0, inf_count = 0;
            for (int i = 0; i < NR && i < 10; ++i) 
            {
                for (int j = 0; j < NTH && j < 5; ++j) 
                {
                    for (int k = 0; k < NZ && k < 5; ++k) 
                    {
                        if (std::isnan(theta[i][j][k])) nan_count++;
                        if (std::isinf(theta[i][j][k])) inf_count++;
                        if (theta[i][j][k] < theta_min) theta_min = theta[i][j][k];
                        if (theta[i][j][k] > theta_max) theta_max = theta[i][j][k];
                        if (u[i][j][k] < u_min) u_min = u[i][j][k];
                        if (u[i][j][k] > u_max) u_max = u[i][j][k];
                    }
                }
            }
            std::cout << "[TIME STEP DEBUG] Step " << steps << " (t=" << simulation_time << "s):" << std::endl;
            std::cout << "  Theta sample: min=" << theta_min << "K, max=" << theta_max << "K" << std::endl;
            std::cout << "  Wind (u) sample: min=" << u_min << "m/s, max=" << u_max << "m/s" << std::endl;
            std::cout << "  NaN/Inf count (sample): " << nan_count << "/" << inf_count << std::endl;

            if (theta_min < 0 || theta_max > 500 || std::abs(u_min) > 150 || std::abs(u_max) > 150) 
            {
                std::cerr << "  ⚠️  WARNING: Values going wrong at step " << steps << "!" << std::endl;
            }
            std::cout << std::endl;
        }
#ifdef EXPORT_NPY
        if (export_ms > 0)
        {
            auto now = std::chrono::steady_clock::now();
            if (std::chrono::duration_cast<std::chrono::milliseconds>(now - lastGuiExport).count() >= export_ms)
            {
                lastGuiExport = now;
                std::string tmp = std::string("data/.tracer_slice_th0.npy.tmp");
                std::string fin = std::string("data/tracer_slice_th0.npy");
                save_field_slice_npy(tracer, thetaIndex, tmp);
                std::rename(tmp.c_str(), fin.c_str());
            }
        }
        if (write_every_s > 0)
        {
            while (simulation_time + 1.0e-9 >= next_field_export_time_s)
            {
                bool export_ok = true;
                timed_call(perf_totals.export_s, [&] {
                    export_ok = write_all_fields(export_index++);
                });
                if (!export_ok)
                {
                    std::cerr << "[VALIDATION] " << validation_error_message << std::endl;
                    return 2;
                }
                next_field_export_time_s += static_cast<double>(write_every_s);
            }
        }
#endif
        ++steps;
        if (global_perf_timing_enabled)
        {
            const auto step_t1 = PerfClock::now();
            perf_totals.total_step_s += std::chrono::duration<double>(step_t1 - step_t0).count();
            ++perf_totals.steps;

            if (global_perf_report_every_steps > 0 &&
                (perf_totals.steps % static_cast<uint64_t>(global_perf_report_every_steps) == 0))
            {
                const double step_ms = (perf_totals.total_step_s / std::max<uint64_t>(1, perf_totals.steps)) * 1000.0;
                const double sim_seconds_per_wall_second =
                    perf_totals.simulated_s / std::max(1e-9, perf_totals.total_step_s);
                std::cout << "[PERF] steps=" << perf_totals.steps
                          << ", avg_step_ms=" << step_ms
                          << ", sim_s_per_wall_s=" << sim_seconds_per_wall_second
                          << std::endl;
            }
        }
        if (steps % 1000 == 0) { }
        if (duration_s >= 0 && simulation_time + 1.0e-9 >= static_cast<double>(duration_s))
        {
            break;
        }
    }

    if (write_every_s > 0 && !export_validation_reports.empty())
    {
        std::filesystem::path summary_path = std::filesystem::path(outdir) / "validation_summary.json";
        if (!global_validation_report_path.empty())
        {
            summary_path = std::filesystem::path(global_validation_report_path) / "validation_summary.json";
        }

        if (!write_validation_summary(summary_path))
        {
            std::cerr << "[VALIDATION] " << validation_error_message << std::endl;
            return 2;
        }
    }

    if (global_perf_timing_enabled && perf_totals.steps > 0)
    {
        const double step_total = std::max(perf_totals.total_step_s, 1e-9);
        const double total_profiled = step_total + perf_totals.initial_export_s;
        const auto pct_step = [&](double component) { return 100.0 * component / step_total; };
        std::cout << "\n[PERF SUMMARY] steps=" << perf_totals.steps
                  << ", step_wall_s=" << step_total
                  << ", avg_step_ms=" << (1000.0 * step_total / perf_totals.steps)
                  << ", sim_s_per_wall_s="
                  << (perf_totals.simulated_s / step_total)
                  << std::endl;
        std::cout << "  radiation_s=" << perf_totals.radiation_s << " (" << pct_step(perf_totals.radiation_s) << "%)" << std::endl;
        std::cout << "  boundary_layer_s=" << perf_totals.boundary_layer_s << " (" << pct_step(perf_totals.boundary_layer_s) << "%)" << std::endl;
        std::cout << "  chaos_noise_s=" << perf_totals.chaos_noise_s << " (" << pct_step(perf_totals.chaos_noise_s) << "%)" << std::endl;
        std::cout << "  chaos_tendency_s=" << perf_totals.chaos_tendency_s << " (" << pct_step(perf_totals.chaos_tendency_s) << "%)" << std::endl;
        std::cout << "  dynamics_s=" << perf_totals.dynamics_s << " (" << pct_step(perf_totals.dynamics_s) << "%)" << std::endl;
        std::cout << "  export_s=" << perf_totals.export_s << " (" << pct_step(perf_totals.export_s) << "%)" << std::endl;
        
        if (perf_totals.initial_export_s > 0.0)
        {
            std::cout << "  initial_export_s=" << perf_totals.initial_export_s << std::endl;
            std::cout << "  total_profiled_s=" << total_profiled << std::endl;
        }
        log_advection_perf_summary();
        std::cout << std::endl;
    }

    return 0;
}
