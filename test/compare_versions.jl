# =============================================================================
# test/compare_versions.jl
#
# Verifies that the refactored simulation() produces bit-for-bit identical
# time-series output to the original version, for every supported constraint
# mode, using a fixed random seed.
#
# The only logic difference between old and new:
#   OLD: get_margin_gene(G, Freq, L)       ← L is a free variable (global)
#   NEW: get_margin_gene(G, Freq, args["L"])← L read from args dict
# Both resolve to the same integer in normal usage; this script confirms it.
#
# Run with:   julia test/compare_versions.jl
# =============================================================================

include(joinpath(@__DIR__, "..", "src", "simulation_lightning.jl"))
using .simulation_light
using DataStructures, Random

# Global L so the old function's free-variable reference resolves correctly
L = 30   # must match args["L"] below

# ── Old simulation() — verbatim logic from git HEAD ───────────────────────────
function simulation_old(dict_genenotype2index, G, Freq_host, Freq_phage, args, T, relaxation_time;
                        stop=false, output_details=true, del_old=true)
    # NOTE: uses global `L` for the stop-condition margin call (original bug)
    if output_details
        TimeSeries_genotype_host  = zeros(Int32, T - relaxation_time, size(G, 1))
        TimeSeries_genotype_phage = zeros(Int32, T - relaxation_time, size(G, 1))
    end
    for t in 1:T
        NX = sum(Freq_host)
        if args["rR_random"] != 0
            G_host_new, FreqG_outflow_host, FreqG_inflow_host =
                recombination_random(dict_genenotype2index, G, Freq_host, args["L"], args["rR_random"], NX, del_old=del_old)
        end
        if args["rR_hosthost"] != 0
            G_host_new, FreqG_outflow_host, FreqG_inflow_host =
                recombination(dict_genenotype2index, G, Freq_host, args["L"], args["rR_hosthost"], NX, del_old=del_old)
        end
        NY = sum(Freq_phage)
        if args["rH_random"] != 0
            G_phage_new =
                recombination_random(dict_genenotype2index, G, Freq_phage, args["L"], args["rH_random"], NY, del_old=del_old)
        end
        if args["rH_phagephage"] != 0
            G_phage_new, FreqG_outflow_phage, FreqG_inflow_phage =
                recombination(dict_genenotype2index, G, Freq_phage, args["L"], args["rH_phagephage"], NY, del_old=del_old)
        end
        if args["rH_phagehost"] != 0
            G_phage_new, FreqG_outflow_phage, FreqG_inflow_phage =
                HorizontalGeneTransfer_2body(dict_genenotype2index, G, Freq_host, Freq_phage, args["L"], args["rH_phagehost"], NY, del_old=del_old)
        end
        if args["rH_3body"] != 0
            G_phage_new, G_phage_sampled =
                HorizontalGeneTransfer_3body_ImportantSampling(dict_genenotype2index, G, Freq_host, Freq_phage, args["L"], args["rH_3body"], NY)
        end
        if args["population_constraint"] == "hard"
            Fit_host, Fit_phage = cal_fitness_hard_constraint(Freq_host, Freq_phage, args["J"], args["Jprop"] * args["J"])
            Freq_host,  _ = select_hard_constraint(Freq_host,  Fit_host,  args["NX"])
            Freq_phage, _ = select_hard_constraint(Freq_phage, Fit_phage, args["NY"])
        elseif args["population_constraint"] == "soft"
            Fit_host, Fit_phage = cal_fitness_no_constraint(Freq_host, Freq_phage, args["J"], args["J"] * args["Jprop"], args["NX"], args["NY"])
            Freq_host,  _ = select_soft_constraint(Freq_host,  Fit_host,  args["soft_threshold"] * args["NX"])
            Freq_phage, _ = select_soft_constraint(Freq_phage, Fit_phage, args["soft_threshold"] * args["NY"])
        elseif args["population_constraint"] == "no"
            Fit_host, Fit_phage = cal_fitness_no_constraint(Freq_host, Freq_phage, args["J"], args["J"] * args["Jprop"], args["NX"], args["NY"])
            Freq_host,  _ = select_no_constraint(Freq_host,  Fit_host)
            Freq_phage, _ = select_no_constraint(Freq_phage, Fit_phage)
        end
        # ← original: free variable L (resolves to global)
        margin_host  = get_margin_gene(G, Freq_host,  L)
        margin_phage = get_margin_gene(G, Freq_phage, L)
        if output_details && t - relaxation_time > 0
            TimeSeries_genotype_host[t - relaxation_time, :]  = Freq_host
            TimeSeries_genotype_phage[t - relaxation_time, :] = Freq_phage
        end
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

# ── New simulation() — refactored version from updated main.jl ────────────────
function simulation_new(dict_genenotype2index, G, Freq_host, Freq_phage, args, T, relaxation_time;
                        stop=false, output_details=true, del_old=true)
    Lval = args["L"]   # ← fix: read from args, not free variable
    if output_details
        TimeSeries_genotype_host  = zeros(Int32, T - relaxation_time, size(G, 1))
        TimeSeries_genotype_phage = zeros(Int32, T - relaxation_time, size(G, 1))
    end
    for t in 1:T
        NX = sum(Freq_host)
        if args["rR_random"] != 0
            G_host_new, FreqG_outflow_host, FreqG_inflow_host =
                recombination_random(dict_genenotype2index, G, Freq_host, Lval, args["rR_random"], NX, del_old=del_old)
        end
        if args["rR_hosthost"] != 0
            G_host_new, FreqG_outflow_host, FreqG_inflow_host =
                recombination(dict_genenotype2index, G, Freq_host, Lval, args["rR_hosthost"], NX, del_old=del_old)
        end
        NY = sum(Freq_phage)
        if args["rH_random"] != 0
            G_phage_new =
                recombination_random(dict_genenotype2index, G, Freq_phage, Lval, args["rH_random"], NY, del_old=del_old)
        end
        if args["rH_phagephage"] != 0
            G_phage_new, FreqG_outflow_phage, FreqG_inflow_phage =
                recombination(dict_genenotype2index, G, Freq_phage, Lval, args["rH_phagephage"], NY, del_old=del_old)
        end
        if args["rH_phagehost"] != 0
            G_phage_new, FreqG_outflow_phage, FreqG_inflow_phage =
                HorizontalGeneTransfer_2body(dict_genenotype2index, G, Freq_host, Freq_phage, Lval, args["rH_phagehost"], NY, del_old=del_old)
        end
        if args["rH_3body"] != 0
            G_phage_new, G_phage_sampled =
                HorizontalGeneTransfer_3body_ImportantSampling(dict_genenotype2index, G, Freq_host, Freq_phage, Lval, args["rH_3body"], NY)
        end
        if args["population_constraint"] == "hard"
            Fit_host, Fit_phage = cal_fitness_hard_constraint(Freq_host, Freq_phage, args["J"], args["Jprop"] * args["J"])
            Freq_host,  _ = select_hard_constraint(Freq_host,  Fit_host,  args["NX"])
            Freq_phage, _ = select_hard_constraint(Freq_phage, Fit_phage, args["NY"])
        elseif args["population_constraint"] == "soft"
            Fit_host, Fit_phage = cal_fitness_no_constraint(Freq_host, Freq_phage, args["J"], args["J"] * args["Jprop"], args["NX"], args["NY"])
            Freq_host,  _ = select_soft_constraint(Freq_host,  Fit_host,  args["soft_threshold"] * args["NX"])
            Freq_phage, _ = select_soft_constraint(Freq_phage, Fit_phage, args["soft_threshold"] * args["NY"])
        elseif args["population_constraint"] == "no"
            Fit_host, Fit_phage = cal_fitness_no_constraint(Freq_host, Freq_phage, args["J"], args["J"] * args["Jprop"], args["NX"], args["NY"])
            Freq_host,  _ = select_no_constraint(Freq_host,  Fit_host)
            Freq_phage, _ = select_no_constraint(Freq_phage, Fit_phage)
        end
        margin_host  = get_margin_gene(G, Freq_host,  Lval)
        margin_phage = get_margin_gene(G, Freq_phage, Lval)
        if output_details && t - relaxation_time > 0
            TimeSeries_genotype_host[t - relaxation_time, :]  = Freq_host
            TimeSeries_genotype_phage[t - relaxation_time, :] = Freq_phage
        end
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

# ── Helper ─────────────────────────────────────────────────────────────────────
function run_both(constraint, seed=42; T=10000, relax=5000, N=1000000, stop=false)
    args = Dict{String,Any}(
        "L" => L, "d" => 2,
        "rR_hosthost" => 1e-4, "rR_random"    => 0.0,
        "rH_random"   => 0.0,  "rH_phagephage"=> 0.0,
        "rH_phagehost"=> 1e-4,  "rH_3body"     => 0.0,
        "NX" => N, "NY" => N,
        "J" => 0.01, "Jprop" => 1.0,
        "population_constraint" => constraint,
        "soft_threshold" => 5.0,
    )
    Random.seed!(seed)
    G1, d2i1, _, Fh1, Fp1 = initial_genotypes(L, N, N, 2)
    TSh_old, TSp_old = simulation_old(d2i1, G1, Fh1, Fp1, args, T, relax;
                                      stop=stop, output_details=true)
    Random.seed!(seed)
    G2, d2i2, _, Fh2, Fp2 = initial_genotypes(L, N, N, 2)
    TSh_new, TSp_new = simulation_new(d2i2, G2, Fh2, Fp2, args, T, relax;
                                      stop=stop, output_details=true)
    host_match  = TSh_old == TSh_new
    phage_match = TSp_old == TSp_new
    max_diff_h  = maximum(abs.(Int.(TSh_old) .- Int.(TSh_new)))
    max_diff_p  = maximum(abs.(Int.(TSp_old) .- Int.(TSp_new)))
    return host_match, phage_match, max_diff_h, max_diff_p
end

# ── Run comparisons across all constraint modes ────────────────────────────────
println("=" ^ 60)
println("  simulation() output comparison: old vs new")
println("  L=$L, N=1000000, T=100000, relaxation=5000, seed=42")
println("=" ^ 60)

results = Bool[]
for constraint in ("hard", "soft", "no")
    hm, pm, dh, dp = run_both(constraint)
    status = (hm && pm) ? "PASS ✓" : "FAIL ✗"
    push!(results, hm && pm)
    println("  constraint=$constraint  host=$(hm ? "match" : "MISMATCH (max_diff=$dh)")  phage=$(pm ? "match" : "MISMATCH (max_diff=$dp)")  → $status")
end

# Also test with stop=true (early-exit path)
let (hm, pm, dh, dp) = run_both("hard", 42; stop=true)
    status = (hm && pm) ? "PASS ✓" : "FAIL ✗"
    push!(results, hm && pm)
    println("  constraint=hard (stop=true)  host=$(hm ? "match" : "MISMATCH")  phage=$(pm ? "match" : "MISMATCH")  → $status")
end

println("=" ^ 60)
println(all(results) ? "All checks passed — outputs are identical." :
                       "Some checks FAILED — outputs differ!")
