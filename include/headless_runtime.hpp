#pragma once

#include <string>

/**
 * @file headless_runtime.hpp
 * @brief Command-line runtime options for non-GUI simulation execution.
 *
 * Defines the option bundle consumed by the headless runner and
 * exposes the entry point used by tests, scripts, and batch jobs.
 * This header is intentionally minimal to keep startup dependencies low.
 */

struct HeadlessRunOptions
{
    int export_ms = 0;
    int duration_s = -1;
    int write_every_s = 0;
    std::string outdir = "data/exports";
};

/**
 * @brief Runs the simulation in headless mode.
 * @param options Runtime options for export cadence and output paths.
 * @return Process-style status code where zero indicates success.
 */
int run_headless_simulation(const HeadlessRunOptions& options);
