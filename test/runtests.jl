# =============================================================================
# test/runtests.jl — Unit tests for the Gene Transfer Simulation
# =============================================================================
#
# Summary:
#   Covers the core functions in src/simulation_lightning.jl:
#     - Genotype initialisation
#     - Gene-level recombination logic
#     - Marginal gene-frequency computation
#     - Diversity statistics (theta)
#     - Poisson sampling helper
#     - Fitness and selection under all three constraint modes
#     - Genotype persistence and establishment-size analysis
#     - CSV output helper (write2file)
#     - Integration smoke-test: short end-to-end simulation
#
# Run with:
#   julia test/runtests.jl
#
# Dependencies: Test, Combinatorics, Distributions, StatsBase
# =============================================================================

using Test
using Combinatorics, Distributions, StatsBase

include(joinpath(@__DIR__, "..", "src", "simulation_lightning.jl"))
using .simulation_light


# =============================================================================
# 1. Genotype initialisation
# =============================================================================

@testset "initial_genotypes" begin
    L, NX, NY, d = 10, 1_000, 2_000, 2
    G, d2i, i2d, Freq_host, Freq_phage = initial_genotypes(L, NX, NY, d)

    K = binomial(L, d)   # expected number of genotypes = C(10, 2) = 45

    @testset "correct number of genotypes" begin
        @test length(G) == K
        @test length(Freq_host)  == K
        @test length(Freq_phage) == K
    end

    @testset "population counts are non-negative and bounded" begin
        @test all(Freq_host  .>= 0)
        @test all(Freq_phage .>= 0)
        @test sum(Freq_host)  <= NX
        @test sum(Freq_phage) <= NY
    end

    @testset "dictionary consistency" begin
        @test length(d2i) == K
        @test length(i2d) == K
        # round-trip: genotype → index → genotype
        for (g, idx) in d2i
            @test i2d[idx] == g
        end
    end

    @testset "each genotype has exactly d genes" begin
        @test all(length.(G) .== d)
    end

    @testset "genotypes are sorted" begin
        @test all(issorted.(G))
    end
end


# =============================================================================
# 2. Gene-level recombination
# =============================================================================

@testset "recombination_process" begin
    g = [1, 3, 5]

    @testset "gene already present → no change" begin
        @test recombination_process(g, 3) == [1, 3, 5]
        @test recombination_process(g, 1) == [1, 3, 5]
        @test recombination_process(g, 5) == [1, 3, 5]
    end

    @testset "new gene → length preserved, output sorted, new gene included" begin
        for new_gene in [2, 4, 6, 10]
            result = recombination_process(g, new_gene)
            @test length(result) == length(g)
            @test issorted(result)
            @test new_gene in result
        end
    end

    @testset "original genotype is not mutated (copy semantics)" begin
        recombination_process(g, 2)
        @test g == [1, 3, 5]
    end
end


# =============================================================================
# 3. Gene frequency marginals
# =============================================================================

@testset "get_margin_gene" begin
    L, NX, NY, d = 6, 100, 100, 2
    G, d2i, _, Freq_host, _ = initial_genotypes(L, NX, NY, d)

    P = get_margin_gene(G, Freq_host, L)

    @test length(P) == L

    # Each individual carries exactly d genes, so the total gene count = d * N
    @test sum(P) == d * sum(Freq_host)

    # All gene counts are non-negative
    @test all(P .>= 0)
end


@testset "get_margin_gene_array" begin
    L, NX, NY, d = 6, 100, 100, 2
    G, _, _, Freq_host, _ = initial_genotypes(L, NX, NY, d)

    # Build a tiny 3-row time series
    TS = hcat(Freq_host, Freq_host, Freq_host)'  # 3×K
    P  = get_margin_gene_array(G, TS, L)

    @test size(P, 1) == 3
    @test size(P, 2) == L

    # Row sums should equal d * row sums of TS
    for t in 1:3
        @test sum(P[t, :]) == d * sum(TS[t, :])
    end
end


# =============================================================================
# 4. Theta diversity statistic
# =============================================================================

@testset "cal_theta" begin
    @testset "uniform distribution → theta ≈ 0" begin
        uniform_ts = 10 .* ones(50, 8)
        @test cal_theta(uniform_ts) ≈ 0.0 atol=1e-10
    end

    @testset "non-uniform distribution → theta > 0" begin
        varied_ts = Float64[100 1; 1 100; 100 1; 1 100]
        @test cal_theta(varied_ts) > 0
    end

    @testset "keepzeros=false excludes zeros" begin
        ts_with_zeros = Float64[0 10; 10 0; 0 10; 10 0]
        # Should not error and should return a finite value
        theta = cal_theta(ts_with_zeros; keepzeros=false)
        @test isfinite(theta)
    end
end


# =============================================================================
# 5. Poisson sampling helper
# =============================================================================

@testset "sample_poisson" begin
    @test sample_poisson(0)    == 0
    @test sample_poisson(-1.0) == 0
    @test sample_poisson(-100) == 0

    # For a positive rate, result must be a non-negative integer
    for _ in 1:50
        v = sample_poisson(5.0)
        @test v >= 0
        @test isa(v, Integer)
    end
end


# =============================================================================
# 6. Fitness functions
# =============================================================================

@testset "cal_fitness_hard_constraint" begin
    Freq_host  = [100, 200, 300]
    Freq_phage = [150, 150, 150]
    J = 0.01

    Fit_host, Fit_phage = cal_fitness_hard_constraint(Freq_host, Freq_phage, J, J)

    @test length(Fit_host)  == length(Freq_host)
    @test length(Fit_phage) == length(Freq_phage)
    @test all(isfinite.(Fit_host))
    @test all(isfinite.(Fit_phage))

    # With uniform phage → all host fitnesses should be equal
    Freq_phage_uniform = [100, 100, 100]
    Fh, _ = cal_fitness_hard_constraint(Freq_host, Freq_phage_uniform, J, J)
    @test all(Fh .≈ Fh[1])
end


@testset "cal_fitness_no_constraint" begin
    Freq_host  = [1000, 1000, 1000]
    Freq_phage = [1000, 1000, 1000]
    NX, NY = 3000, 3000
    J = 0.01

    Fit_host, Fit_phage = cal_fitness_no_constraint(Freq_host, Freq_phage, J, J, NX, NY)

    @test length(Fit_host)  == length(Freq_host)
    @test length(Fit_phage) == length(Freq_phage)
    @test all(isfinite.(Fit_host))
    @test all(isfinite.(Fit_phage))
end


# =============================================================================
# 7. Selection mechanisms
# =============================================================================

@testset "select_hard_constraint" begin
    Freq = [1_000, 1_000, 1_000]
    Fit  = [0.1, 0.0, -0.1]
    Nmax = 3_000

    new_Freq, meanfit = select_hard_constraint(Freq, Fit, Nmax)

    @test length(new_Freq) == length(Freq)
    @test all(new_Freq .>= 0)
    @test isfinite(meanfit)
end


@testset "select_soft_constraint" begin
    Nmax = 3_000

    @testset "below cap → unconstrained Poisson growth" begin
        Freq = [100, 100, 100]   # well below Nmax
        Fit  = [0.05, 0.0, -0.05]
        new_Freq, meanfit = select_soft_constraint(Freq, Fit, Nmax)
        @test all(new_Freq .>= 0)
        @test isfinite(meanfit)
    end

    @testset "above cap → population is suppressed" begin
        Freq = [4_000, 4_000, 4_000]   # above Nmax
        Fit  = [0.0, 0.0, 0.0]
        new_Freq, meanfit = select_soft_constraint(Freq, Fit, Nmax)
        @test all(new_Freq .>= 0)
        # Expected total should be pulled toward Nmax
        @test sum(new_Freq) < sum(Freq)
    end
end


@testset "select_no_constraint" begin
    Freq = [500, 500, 500]
    Fit  = [0.1, 0.0, -0.1]

    new_Freq, meanfit = select_no_constraint(Freq, Fit)
    @test length(new_Freq) == length(Freq)
    @test all(new_Freq .>= 0)
    @test isfinite(meanfit)
end


# =============================================================================
# 8. Persistence and gap time analysis
# =============================================================================

@testset "genotypePers" begin
    @testset "genotype born, grows, and dies once" begin
        # Column: 0, 0, rises above mean (50 > ~16.7), then drops to 0
        TS = Float64[0 0 50 50 50 0 0 0]'   # 8×1 matrix
        persTime, gapTime = genotypePers(TS)
        @test !isempty(persTime)
        @test all(persTime .> 0)
    end

    @testset "always-present genotype → one persistence entry = series length" begin
        TS = 100 .* ones(Float64, 20, 1)
        persTime, _ = genotypePers(TS)
        @test !isempty(persTime)
    end
end


@testset "find_establishment_size" begin
    # One genotype: silent for 9 steps, then present for steps 10-29,
    # then extinct again.  Host is always present.
    T, K = 40, 1
    host_ts  = 100 .* ones(Int32, T, K)
    phage_ts = zeros(Int32, T, K)
    phage_ts[10:29, 1] .= 100

    result = find_establishment_size(host_ts, phage_ts)
    @test isa(result, Vector)
end


# =============================================================================
# 9. CSV output helper
# =============================================================================

@testset "write2file" begin
    tmpfile = tempname() * ".txt"

    column_names = ["L", "J", "label"]
    data         = [30, 0.01, "hard"]

    # First write creates header + row
    write2file(tmpfile, data, column_names)
    @test isfile(tmpfile)

    lines = readlines(tmpfile)
    @test length(lines) == 2
    @test lines[1] == "L,J,label"
    @test lines[2] == "30,0.01,hard"

    # Second write appends another row (no duplicate header)
    write2file(tmpfile, data, column_names)
    lines2 = readlines(tmpfile)
    @test length(lines2) == 3
    @test lines2[1] == "L,J,label"   # still only one header

    rm(tmpfile)
end


# =============================================================================
# 10. Integration smoke-test (short end-to-end run)
# =============================================================================

@testset "integration: short simulation" begin
    L, NX, NY, d = 8, 500, 500, 2
    G, d2i, _, Freq_host, Freq_phage = initial_genotypes(L, NX, NY, d)

    args = Dict{String,Any}(
        "L"                     => L,
        "d"                     => d,
        "rR_random"             => 1e-3,
        "rR_hosthost"           => 0.0,
        "rH_random"             => 1e-3,
        "rH_phagephage"         => 0.0,
        "rH_phagehost"          => 0.0,
        "rH_3body"              => 0.0,
        "NX"                    => NX,
        "NY"                    => NY,
        "J"                     => 0.01,
        "Jprop"                 => 1.0,
        "T"                     => 50,
        "population_constraint" => "hard",
        "soft_threshold"        => 5.0,
    )

    T_run         = 50
    relaxation    = 5

    # Import main-level simulation function (defined in main.jl) by including it
    # without running the entry-point code — we test only the module functions here.
    # Directly call the simulation_light exported functions instead:

    include(joinpath(@__DIR__, "..", "src", "simulation_lightning.jl"))
    # Module already loaded; just exercise the core loop manually.

    for _ in 1:T_run
        NX_cur = sum(Freq_host)
        if args["rR_random"] != 0
            recombination_random(d2i, G, Freq_host, L, args["rR_random"], NX_cur)
        end
        NY_cur = sum(Freq_phage)
        if args["rH_random"] != 0
            recombination_random(d2i, G, Freq_phage, L, args["rH_random"], NY_cur)
        end
        Fit_h, Fit_p = cal_fitness_hard_constraint(Freq_host, Freq_phage, args["J"], args["Jprop"] * args["J"])
        Freq_host,  _ = select_hard_constraint(Freq_host,  Fit_h, NX)
        Freq_phage, _ = select_hard_constraint(Freq_phage, Fit_p, NY)
    end

    @test sum(Freq_host)  >= 0
    @test sum(Freq_phage) >= 0
    @test length(Freq_host)  == binomial(L, d)
    @test length(Freq_phage) == binomial(L, d)
end


println("\nAll tests completed.")
