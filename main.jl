# -*- coding: utf-8 -*-
# =============================================================================
# main.jl — Entry point for the Gene Transfer Simulation
# =============================================================================
#
# Summary:
#   Parses command-line arguments (falling back to config.json defaults),
#   initialises genotype pools, runs the stochastic simulation for T
#   generations, and writes summary statistics to a CSV output file.
#
#   Three top-level functions are defined:
#     - parse_commandline() : merges config.json defaults with CLI overrides
#     - simulation()        : runs the main stochastic loop
#     - AnalyzeData()       : computes and saves summary statistics
#     - main()              : orchestrates the full pipeline
#
# Usage:
#   julia main.jl --T 200000 --L 40 --J 0.01 --N 1E6 --rR_random 3E-5 --rH_random 3E-5
#
# Merge multiple output files:
#   awk '!/^(L,)/' results/*.txt > combined_output.txt
#
# Author: Wenping Cui
# =============================================================================

include("src/simulation_lightning.jl")
using .simulation_light
using ArgParse, BenchmarkTools, Combinatorics, DataStructures, JSON, Printf, SpecialFunctions, StatsBase

# Path to the default-parameter configuration file (same directory as this script)
const CONFIG_FILE = joinpath(@__DIR__, "config.json")


# -----------------------------------------------------------------------------
# Configuration loading
# -----------------------------------------------------------------------------

"""
    load_config(path)

Load parameter defaults from a JSON file.  Returns an empty Dict if the
file is not found, allowing the hard-coded ArgParse defaults to take over.
"""
function load_config(path::String=CONFIG_FILE)
    isfile(path) || return Dict{String,Any}()
    return JSON.parsefile(path)
end


"""
    cfg_get(config, section, key, fallback)

Look up `config[section][key]`, returning `fallback` if either the section
or the key is absent.  Keeps default retrieval concise in `parse_commandline`.
"""
function cfg_get(config, section, key, fallback)
    section_dict = get(config, section, Dict())
    return get(section_dict, key, fallback)
end


# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------

"""
    parse_commandline()

Build an ArgParseSettings table whose defaults are sourced from config.json.
Any value explicitly supplied on the command line overrides the config file.

Returns an ArgParse result Dict with all simulation parameters.
"""
function parse_commandline()
    config = load_config()

    s = ArgParseSettings()
    @add_arg_table s begin
        "--L"
            help     = "total number of genes"
            arg_type = Int64
            default  = Int64(cfg_get(config, "population", "L", 30))
        "--d"
            help     = "number of genes in each genotype"
            arg_type = Int64
            default  = Int64(cfg_get(config, "population", "d", 2))
        "--r"
            help     = "host and phage recombination rate (legacy shorthand)"
            arg_type = Float64
            default  = Float64(cfg_get(config, "recombination_rates", "r", 1e-3))
        "--J"
            help     = "interaction strength (host-phage coupling)"
            arg_type = Float64
            default  = Float64(cfg_get(config, "genetics", "J", 5.0))
        "--Jprop"
            help     = "ratio between phage and host interaction strength"
            arg_type = Float64
            default  = Float64(cfg_get(config, "genetics", "Jprop", 1.0))
        "--N"
            help     = "population size (hosts; phages scaled by Nprop_phagedivhost)"
            arg_type = Float64
            default  = Float64(cfg_get(config, "population", "N", 1.0e6))
        "--T"
            help     = "number of simulation generations"
            arg_type = Float64
            default  = Float64(cfg_get(config, "simulation", "T", 5.0e5))
        "--taskID"
            help     = "task ID (for batch job arrays)"
            arg_type = Int64
            default  = Int64(cfg_get(config, "run", "taskID", 1))
        "--threadID"
            help     = "thread ID (string, for output file naming)"
            arg_type = String
            default  = string(cfg_get(config, "run", "threadID", "1"))
        "--trialID"
            help     = "trial / replicate ID"
            arg_type = Int64
            default  = Int64(cfg_get(config, "run", "trialID", 1))
        "--rR_hosthost"
            help     = "host-host (sexual) recombination rate"
            arg_type = Float64
            default  = Float64(cfg_get(config, "recombination_rates", "rR_hosthost", 0.0))
        "--rR_random"
            help     = "host random-mutation recombination rate"
            arg_type = Float64
            default  = Float64(cfg_get(config, "recombination_rates", "rR_random", 0.0))
        "--rH_random"
            help     = "phage random-mutation recombination rate"
            arg_type = Float64
            default  = Float64(cfg_get(config, "recombination_rates", "rH_random", 0.0))
        "--rH_phagephage"
            help     = "phage-phage recombination rate"
            arg_type = Float64
            default  = Float64(cfg_get(config, "recombination_rates", "rH_phagephage", 0.0))
        "--rH_phagehost"
            help     = "phage-host horizontal gene transfer (2-body) rate"
            arg_type = Float64
            default  = Float64(cfg_get(config, "recombination_rates", "rH_phagehost", 0.0))
        "--rH_3body"
            help     = "phage coinfection (3-body) horizontal gene transfer rate"
            arg_type = Float64
            default  = Float64(cfg_get(config, "recombination_rates", "rH_3body", 0.0))
        "--Nprop_phagedivhost"
            help     = "ratio of phage to host abundance (NY = N * Nprop_phagedivhost)"
            arg_type = Float64
            default  = Float64(cfg_get(config, "population", "Nprop_phagedivhost", 1.0))
        "--r_prop"
            help     = "global multiplier applied to all phage recombination rates"
            arg_type = Float64
            default  = Float64(cfg_get(config, "recombination_rates", "r_prop", 1.0))
        "--output_details"
            help     = "if true, record full time-series and compute summary statistics"
            arg_type = Bool
            default  = Bool(cfg_get(config, "io", "output_details", true))
        "--population_constraint"
            help     = "population constraint: 'hard', 'soft', or 'no'"
            arg_type = String
            default  = string(cfg_get(config, "population", "population_constraint", "soft"))
        "--save_dir"
            help     = "directory for output files"
            arg_type = String
            default  = string(cfg_get(config, "io", "save_dir", "results_simulation_test"))
        "--soft_threshold"
            help     = "soft constraint threshold (Nmax = soft_threshold * N)"
            arg_type = Float64
            default  = Float64(cfg_get(config, "population", "soft_threshold", 5.0))
    end
    return parse_args(s)
end


# -----------------------------------------------------------------------------
# Simulation loop
# -----------------------------------------------------------------------------

"""
    simulation(dict_genenotype2index, G, Freq_host, Freq_phage, args, T, relaxation_time;
               stop=false, output_details=true, del_old=true)

Run the stochastic simulation for up to `T` generations.

Each generation applies (if the corresponding rate is non-zero):
  1. Host recombination (random or pairwise)
  2. Phage recombination / horizontal gene transfer
  3. Fitness evaluation and Poisson-sampled selection

The first `relaxation_time` generations are treated as a burn-in and are
not stored in the time series.

# Arguments
- `stop`           : if `true`, halt early when a gene goes to fixation or near-extinction
- `output_details` : if `true`, return full T×K time-series matrices;
                     otherwise return only the stopping time `t`
- `del_old`        : if `true`, remove recombining individuals from the pool before replacement

# Returns
- `(TimeSeries_genotype_host, TimeSeries_genotype_phage)` when `output_details=true`
- `t` (Int) when `output_details=false`
"""
function simulation(dict_genenotype2index, G, Freq_host, Freq_phage, args, T, relaxation_time;
                    stop=false, output_details=true, del_old=true)
    L = args["L"]

    if output_details
        TimeSeries_genotype_host  = zeros(Int32, T - relaxation_time, size(G, 1))
        TimeSeries_genotype_phage = zeros(Int32, T - relaxation_time, size(G, 1))
    end

    for t in 1:T
        # --- Host recombination ---
        NX = sum(Freq_host)
        if args["rR_random"] != 0
            G_host_new, FreqG_outflow_host, FreqG_inflow_host =
                recombination_random(dict_genenotype2index, G, Freq_host, L,
                                     args["rR_random"], NX, del_old=del_old)
        end
        if args["rR_hosthost"] != 0
            G_host_new, FreqG_outflow_host, FreqG_inflow_host =
                recombination(dict_genenotype2index, G, Freq_host, L,
                              args["rR_hosthost"], NX, del_old=del_old)
        end

        # --- Phage recombination / HGT ---
        NY = sum(Freq_phage)
        if args["rH_random"] != 0
            G_phage_new =
                recombination_random(dict_genenotype2index, G, Freq_phage, L,
                                     args["rH_random"], NY, del_old=del_old)
        end
        if args["rH_phagephage"] != 0
            G_phage_new, FreqG_outflow_phage, FreqG_inflow_phage =
                recombination(dict_genenotype2index, G, Freq_phage, L,
                              args["rH_phagephage"], NY, del_old=del_old)
        end
        if args["rH_phagehost"] != 0
            G_phage_new, FreqG_outflow_phage, FreqG_inflow_phage =
                HorizontalGeneTransfer_2body(dict_genenotype2index, G, Freq_host, Freq_phage,
                                             L, args["rH_phagehost"], NY, del_old=del_old)
        end
        if args["rH_3body"] != 0
            G_phage_new, G_phage_sampled =
                HorizontalGeneTransfer_3body_ImportantSampling(dict_genenotype2index, G,
                                                                Freq_host, Freq_phage,
                                                                L, args["rH_3body"], NY)
        end

        # --- Selection ---
        if args["population_constraint"] == "hard"
            Fit_host, Fit_phage =
                cal_fitness_hard_constraint(Freq_host, Freq_phage, args["J"], args["Jprop"] * args["J"])
            Freq_host,  meanfit_host  = select_hard_constraint(Freq_host,  Fit_host,  args["NX"])
            Freq_phage, meanfit_phage = select_hard_constraint(Freq_phage, Fit_phage, args["NY"])

        elseif args["population_constraint"] == "soft"
            Fit_host, Fit_phage =
                cal_fitness_no_constraint(Freq_host, Freq_phage, args["J"], args["J"] * args["Jprop"],
                                          args["NX"], args["NY"])
            Freq_host,  meanfit_host  = select_soft_constraint(Freq_host,  Fit_host,  args["soft_threshold"] * args["NX"])
            Freq_phage, meanfit_phage = select_soft_constraint(Freq_phage, Fit_phage, args["soft_threshold"] * args["NY"])

        elseif args["population_constraint"] == "no"
            Fit_host, Fit_phage =
                cal_fitness_no_constraint(Freq_host, Freq_phage, args["J"], args["J"] * args["Jprop"],
                                          args["NX"], args["NY"])
            Freq_host,  meanfit_host  = select_no_constraint(Freq_host,  Fit_host)
            Freq_phage, meanfit_phage = select_no_constraint(Freq_phage, Fit_phage)
        end

        # --- Record time series (post-relaxation) ---
        if output_details && t - relaxation_time > 0
            TimeSeries_genotype_host[t - relaxation_time, :]  = Freq_host
            TimeSeries_genotype_phage[t - relaxation_time, :] = Freq_phage
        end

        # --- Optional early-stopping check ---
        margin_host  = get_margin_gene(G, Freq_host,  L)
        margin_phage = get_margin_gene(G, Freq_phage, L)
        if stop && (sum(margin_host .== 0.) > 0 || sum(margin_phage .> 0) < 3)
            if output_details
                return TimeSeries_genotype_host[1:t-1, :], TimeSeries_genotype_phage[1:t-1, :]
            else
                return t
            end
        end
    end

    if output_details
        return TimeSeries_genotype_host, TimeSeries_genotype_phage
    else
        return T
    end
end


# -----------------------------------------------------------------------------
# Post-simulation analysis and output
# -----------------------------------------------------------------------------

"""
    compute_stats(args, TimeSeries_genotype_host, TimeSeries_genotype_phage,
                  G, L, relaxation_time) → (column_names, save_data)

Compute all summary statistics from time-series matrices and return them as
`(column_names, save_data)` vectors ready for `write2file`.  Separated from
I/O so the result can be inspected or tested without touching the filesystem.

Frees the large genotype matrices before computing gene-level statistics to
reduce peak memory usage.
"""
function compute_stats(args, TimeSeries_genotype_host, TimeSeries_genotype_phage,
                       G, L, relaxation_time)
    T = size(TimeSeries_genotype_host, 1) + relaxation_time

    # Cache row-sums once to avoid recomputing for both mean and std
    host_sums  = sum(TimeSeries_genotype_host,  dims=2)
    phage_sums = sum(TimeSeries_genotype_phage, dims=2)
    mean_Host_N,  std_Host_N  = mean(host_sums),  std(host_sums)
    mean_Phage_N, std_Phage_N = mean(phage_sums), std(phage_sums)

    estab_size      = find_establishment_size(TimeSeries_genotype_host, TimeSeries_genotype_phage)
    estab_size_mean = isempty(estab_size) ? 0 : mean(estab_size)

    persTimeHost, gapTimeHost   = genotypePers(TimeSeries_genotype_host)
    persTimePhage, gapTimePhage = genotypePers(TimeSeries_genotype_phage)

    persTimeHost_mean  = isempty(persTimeHost)  ? 0 : mean(persTimeHost)
    persTimePhage_mean = isempty(persTimePhage) ? 0 : mean(persTimePhage)
    gapTimeHost_mean   = isempty(gapTimeHost)   ? 0 : mean(gapTimeHost)
    gapTimePhage_mean  = isempty(gapTimePhage)  ? 0 : mean(gapTimePhage)

    theta_GH             = cal_theta(TimeSeries_genotype_host)
    TimeSeries_gene_host = get_margin_gene_array(G, TimeSeries_genotype_host, L)
    TimeSeries_genotype_host = nothing   # release before allocating gene-phage array
    theta_gH             = cal_theta(TimeSeries_gene_host)

    theta_GP              = cal_theta(TimeSeries_genotype_phage)
    TimeSeries_gene_phage = get_margin_gene_array(G, TimeSeries_genotype_phage, L)
    TimeSeries_genotype_phage = nothing  # release before allocating gene-phage array
    theta_gP              = cal_theta(TimeSeries_gene_phage)

    column_names = ["L", "d", "J", "Jprop", "NX", "NY",
                    "rR_random", "rH_random", "rR_hosthost", "rH_phagephage",
                    "rH_phagehost", "rH_3body",
                    "simulation_time", "GenePersisentTime", "mean_establish_size",
                    "persTimeHost", "persTimePhage", "gapTimeHost", "gapTimePhage",
                    "theta_GH", "theta_GP", "theta_gH", "theta_gP",
                    "constraint", "soft_threshold",
                    "mean_host_N", "std_host_N", "mean_phage_N", "std_phage_N"]

    save_data = [args["L"], args["d"], args["J"], args["Jprop"], args["NX"], args["NY"],
                 args["rR_random"], args["rH_random"], args["rR_hosthost"], args["rH_phagephage"],
                 args["rH_phagehost"], args["rH_3body"], args["T"],
                 T, estab_size_mean,
                 persTimeHost_mean, persTimePhage_mean, gapTimeHost_mean, gapTimePhage_mean,
                 theta_GH, theta_GP, theta_gH, theta_gP,
                 args["population_constraint"], args["soft_threshold"],
                 mean_Host_N, std_Host_N, mean_Phage_N, std_Phage_N]

    return column_names, save_data
end


"""
    AnalyzeData(args, TimeSeries_genotype_host, TimeSeries_genotype_phage,
                G, L, save_dir, relaxation_time)

Compute summary statistics via `compute_stats` and append them to a CSV file.
Output filename encodes key parameters for easy identification in batch runs.
"""
function AnalyzeData(args, TimeSeries_genotype_host, TimeSeries_genotype_phage,
                     G, L, save_dir, relaxation_time)
    column_names, save_data = compute_stats(args, TimeSeries_genotype_host,
                                            TimeSeries_genotype_phage, G, L, relaxation_time)
    file_name = "L_$(args["L"])_N_$(@sprintf("%.2E", args["NX"]))_J_$(args["J"])" *
                "_outputTask$(args["taskID"])Thread$(args["threadID"]).txt"
    write2file(joinpath(save_dir, file_name), save_data, column_names)
end


# -----------------------------------------------------------------------------
# Main pipeline
# -----------------------------------------------------------------------------

"""
    main(pars; relaxation_time=5_000)

Orchestrate the full simulation pipeline:
  1. Unpack and store all parameters in an OrderedDict.
  2. Initialise genotypes with `initial_genotypes`.
  3. Run `simulation` for `T` generations.
  4. If `output_details=true`, call `AnalyzeData` to save summary statistics;
     otherwise save only the stopping time.

`pars` is a 20-element tuple in the order defined by the entry-point block below.
"""
function main(pars; relaxation_time=5_000)
    (L, d, rR_random, rH_random, rR_hosthost, rH_phagephage, rH_phagehost, rH_3body,
     N, J, Jprop, Nprop_phagedivhost, T, taskID, threadID, trialID,
     output_details, population_constraint, soft_threshold, save_dir) = pars

    args = OrderedDict{String,Any}(
        "L"                     => Int64(L),
        "d"                     => d,
        "rR_random"             => rR_random,
        "rR_hosthost"           => rR_hosthost,
        "rH_random"             => rH_random,
        "rH_phagephage"         => rH_phagephage,
        "rH_phagehost"          => rH_phagehost,
        "rH_3body"              => rH_3body,
        "NX"                    => Int64(N),
        "NY"                    => N * Nprop_phagedivhost,
        "J"                     => J,
        "Jprop"                 => Jprop,
        "T"                     => Int64(T),
        "taskID"                => taskID,
        "threadID"              => threadID,
        "trialID"               => Int64(trialID),
        "population_constraint" => population_constraint,
        "soft_threshold"        => soft_threshold,
    )

    G, dict_genenotype2index, dict_index2genenotype, Freq_host, Freq_phage =
        initial_genotypes(args["L"], args["NX"], args["NY"], args["d"])

    if output_details
        TimeSeries_genotype_host, TimeSeries_genotype_phage =
            simulation(dict_genenotype2index, G, Freq_host, Freq_phage, args,
                       Int64(T), relaxation_time; stop=true, output_details=true, del_old=true)
        t = size(TimeSeries_genotype_host, 1) + relaxation_time
        AnalyzeData(args, TimeSeries_genotype_host, TimeSeries_genotype_phage,
                    G, Int64(L), save_dir, relaxation_time)
    else
        t = simulation(dict_genenotype2index, G, Freq_host, Freq_phage, args,
                       Int64(T), relaxation_time; stop=true, output_details=false, del_old=true)
        column_names = ["L", "d", "J", "Jprop", "NX", "NY",
                        "rR_random", "rH_random", "rR_hosthost", "rH_phagephage",
                        "rH_phagehost", "rH_3body",
                        "GenePersisentTime", "constraint", "soft_threshold"]
        save_data    = [args["L"], args["d"], args["J"], args["Jprop"], args["NX"], args["NY"],
                        args["rR_random"], args["rH_random"], args["rR_hosthost"], args["rH_phagephage"],
                        args["rH_phagehost"], args["rH_3body"],
                        t, args["population_constraint"], args["soft_threshold"]]
        file_name = "L_$(args["L"])_N_$(@sprintf("%.2E", args["NX"]))_J_$(args["J"])" *
                    "_outputTask$(args["taskID"])Thread$(args["threadID"])_geneduration.txt"
        write2file(joinpath(save_dir, file_name), save_data, column_names)
    end

    println(args)
    println("Stopped at generation: ", t)
end


# =============================================================================
# Entry point
# =============================================================================

parsed_args = parse_commandline()

L                  = parsed_args["L"]
N                  = parsed_args["N"]
J                  = parsed_args["J"]
T                  = floor(Int64, parsed_args["T"])
rR_random          = parsed_args["rR_random"]
rH_random          = parsed_args["rH_random"]
rR_hosthost        = parsed_args["rR_hosthost"]
rH_phagephage      = parsed_args["rH_phagephage"]
rH_phagehost       = parsed_args["rH_phagehost"]
rH_3body           = parsed_args["rH_3body"]
Jprop              = parsed_args["Jprop"]
Nprop_phagedivhost = parsed_args["Nprop_phagedivhost"]
d                  = parsed_args["d"]
r_prop             = parsed_args["r_prop"]

# Apply the global phage-rate multiplier
rH_random, rH_phagephage, rH_phagehost, rH_3body =
    r_prop .* [rH_random, rH_phagephage, rH_phagehost, rH_3body]

isdir(parsed_args["save_dir"]) || mkdir(parsed_args["save_dir"])

taskID         = parsed_args["taskID"]
threadID       = parse(Int64, parsed_args["threadID"])
trialID        = parsed_args["trialID"]
output_details = parsed_args["output_details"]

println("threadID: ", threadID)

@time main((L, d, rR_random, rH_random, rR_hosthost, rH_phagephage, rH_phagehost, rH_3body,
             N, J, Jprop, Nprop_phagedivhost, T, taskID, threadID, trialID,
             output_details, parsed_args["population_constraint"], parsed_args["soft_threshold"],
             parsed_args["save_dir"]))
