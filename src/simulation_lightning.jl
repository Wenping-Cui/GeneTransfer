# -*- coding: utf-8 -*-
# =============================================================================
# simulation_lightning.jl — Core simulation library for gene transfer dynamics
# =============================================================================
#
# Summary:
#   Implements the stochastic population-genetics engine for phage-host
#   co-evolution under horizontal gene transfer (HGT).  The module provides:
#     - Genotype initialisation and bookkeeping utilities
#     - Recombination mechanisms (random, pairwise host-host / phage-phage)
#     - Horizontal gene transfer (2-body phage-host, 3-body coinfection)
#     - Fitness evaluation and selection under hard / soft / no constraints
#     - Summary statistics: theta diversity, persistence times, establishment sizes
#
# Author: Wenping Cui
# Created: 2022-10-10
# =============================================================================

module simulation_light

    using Distributions, StatsBase, Combinatorics

    export initial_genotypes,
           recombination, recombination_random, recombination_process,
           HorizontalGeneTransfer_2body,
           HorizontalGeneTransfer_3body_Direct,
           HorizontalGeneTransfer_3body_ImportantSampling,
           get_margin_gene, get_margin_gene_array,
           add_newgenotypes,
           find_establishment_size, genotypePers,
           write2file,
           sample_poisson,
           cal_fitness_no_constraint, cal_fitness_hard_constraint,
           select_no_constraint, select_soft_constraint, select_hard_constraint,
           cal_theta, vecvec2matrix


    # -------------------------------------------------------------------------
    # I/O utilities
    # -------------------------------------------------------------------------

    """
        write2file(filename, save_data, column_name)

    Append a single data row to a CSV file.  If the file does not yet exist,
    a header line is written first.

    # Arguments
    - `filename`    : full path to the output file
    - `save_data`   : vector of scalar values to write (one per column)
    - `column_name` : vector of column header strings (same length as `save_data`)
    """
    function write2file(filename, save_data, column_name)
        if !isfile(filename)
            open(filename, "w") do io
                print(io, column_name[1])
                for d in column_name[2:end]
                    print(io, ",", d)
                end
                print(io, "\n")
            end
        end

        open(filename, "a") do io
            print(io, save_data[1])
            for d in save_data[2:end]
                print(io, ",", d)
            end
            print(io, "\n")
        end
    end


    # -------------------------------------------------------------------------
    # Array / matrix helpers
    # -------------------------------------------------------------------------

    """
        vecvec2matrix(M)

    Convert a vector-of-vectors `M` into a 2-D matrix by stacking rows.
    """
    function vecvec2matrix(M)
        return reduce(vcat, transpose.(M))
    end


    # -------------------------------------------------------------------------
    # Diversity statistics
    # -------------------------------------------------------------------------

    """
        cal_theta(TimeSeries; keepzeros=true)

    Compute the theta (θ) diversity statistic from a time-series matrix.

    θ is defined as the mean squared deviation from the mean, normalised by
    the mean: `mean((x - μ)² / μ)` over all elements.

    # Arguments
    - `TimeSeries` : T×K matrix of population counts (rows = time, cols = genotypes/genes)
    - `keepzeros`  : if `false`, zero entries are excluded before computing statistics
    """
    function cal_theta(TimeSeries; keepzeros=true)
        data = vec(TimeSeries[1:end, :])
        if !keepzeros
            data = data[data .> 0]
        end
        return mean((data .- mean(data)) .^ 2 ./ mean(data))
    end


    # -------------------------------------------------------------------------
    # Initialisation
    # -------------------------------------------------------------------------

    """
        initial_genotypes(L, NX, NY, d)

    Enumerate all C(L, d) genotypes and distribute the initial population
    uniformly across them.

    # Arguments
    - `L`  : total number of genes available
    - `NX` : initial host population size
    - `NY` : initial phage population size
    - `d`  : number of genes carried by each genotype

    # Returns
    - `G`                    : vector of genotypes (each is a sorted `d`-element vector)
    - `dict_genenotype2index`: Dict mapping genotype → frequency-vector index
    - `dict_index2genenotype`: Dict mapping frequency-vector index → genotype
    - `Freq_host`            : initial host population counts (one entry per genotype)
    - `Freq_phage`           : initial phage population counts

    # Example
        L, NX, NY, d = 10, 10^6, 10^6, 2
        G, d2i, i2d, Freq_host, Freq_phage = initial_genotypes(L, NX, NY, d)
    """
    function initial_genotypes(L, NX, NY, d)
        dict_genenotype2index = Dict{Array{Int64,1}, Int64}()
        dict_index2genenotype = Dict{Int64, Array{Int64,1}}()
        G = collect(combinations(1:L, d))
        K = size(G, 1)
        for i in 1:K
            dict_genenotype2index[G[i]] = i
            dict_index2genenotype[i]    = G[i]
        end
        Freq_host  = floor(Int64, NX / K) .* vec(ones(Int64, 1, K))
        Freq_phage = floor(Int64, NY / K) .* vec(ones(Int64, 1, K))
        return G, dict_genenotype2index, dict_index2genenotype, Freq_host, Freq_phage
    end


    # -------------------------------------------------------------------------
    # Analysis: persistence and establishment
    # -------------------------------------------------------------------------

    """
        find_establishment_size(TimeSeries_genotype_host, TimeSeries_genotype_phage)

    For each genotype, identify episodes where the phage population is born,
    exceeds the global mean, and then goes extinct.  Record the host population
    size at the birth time of each such episode.

    Returns a vector of establishment sizes (possibly empty).
    """
    function find_establishment_size(TimeSeries_genotype_host, TimeSeries_genotype_phage)
        est_size = []
        ave = 1.4 * mean(TimeSeries_genotype_host)
        for i in 1:size(TimeSeries_genotype_host, 2)
            dataX = TimeSeries_genotype_host[:, i]
            dataY = TimeSeries_genotype_phage[:, i]
            birthtime = 1
            avePass   = false
            birth     = false
            for t in 2:length(dataY)-1
                if dataY[t-1] <= 0.1 && dataY[t] > 0   # birth event
                    birthtime = t
                    birth     = true
                end
                if !avePass && birth && dataY[t] >= ave
                    avePass = true
                end
                if dataY[t] == 0.                        # death event
                    if birth && avePass && birthtime != 1
                        push!(est_size, dataX[birthtime])
                    end
                    avePass = false
                    birth   = false
                end
            end
        end
        return est_size
    end


    """
        genotypePers(TS)

    Compute persistence times (how long each genotype survives per episode)
    and gap times (time between consecutive extinctions of the same genotype).

    Only episodes that exceeded the global mean population are counted.

    Returns `(persTime, gapTime)` as two vectors of integers.
    """
    function genotypePers(TS)
        persTime = []
        gapTime  = []
        ave = mean(TS)
        for i in 1:size(TS, 2)
            data      = TS[:, i]
            birthtime = 1
            deathtime = 1
            avePass   = false
            birth     = false
            count_persTime = 0
            for t in 2:size(TS, 1)
                if data[t-1] <= 0.1 && data[t] > 0   # birth
                    birthtime = t
                    birth     = true
                end
                if !avePass && birth && data[t] >= ave
                    avePass = true
                end
                if data[t] == 0.                      # death
                    if birth && avePass
                        if deathtime != 1
                            push!(gapTime, birthtime - deathtime)
                        end
                        count_persTime += 1
                        push!(persTime, t - birthtime)
                        deathtime = t
                    end
                    avePass = false
                    birth   = false
                end
            end
            if count_persTime == 0
                push!(persTime, size(TS, 1))
            end
        end
        return persTime, gapTime
    end


    # -------------------------------------------------------------------------
    # Internal: genotype frequency bookkeeping
    # -------------------------------------------------------------------------

    """
        count_genotypes(dict_genenotype2index, Freq_like, G_sampled)

    Build a frequency vector from a list of sampled genotypes `G_sampled`,
    using `Freq_like` as a zero template.
    """
    function count_genotypes(dict_genenotype2index, Freq_like, G_sampled)
        Freq = zero(Freq_like)
        for (g, n) in countmap(G_sampled)
            Freq[dict_genenotype2index[g]] = n
        end
        return Freq
    end


    """
        add_newgenotypes(dict_genenotype2index, G_sampled, Freq)

    Increment `Freq` for each genotype in `G_sampled`.  Returns the subset
    of `G_sampled` that had a zero count before (i.e., newly appearing genotypes).
    """
    function add_newgenotypes(dict_genenotype2index, G_sampled, Freq)
        G_new = []
        for g in G_sampled
            i = dict_genenotype2index[g]
            if Freq[i] == 0
                push!(G_new, g)
            end
            Freq[i] += 1
        end
        return G_new
    end


    """
        delete_genotypes(dict_genenotype2index, G_sampled, Freq)

    Decrement `Freq` for each genotype in `G_sampled` (floor at zero).
    """
    function delete_genotypes(dict_genenotype2index, G_sampled, Freq)
        for g in G_sampled
            i = dict_genenotype2index[g]
            if Freq[i] > 0
                Freq[i] -= 1
            end
        end
        return []
    end


    # -------------------------------------------------------------------------
    # Gene frequency marginals
    # -------------------------------------------------------------------------

    """
        get_margin_gene(G, Freq, L)

    Compute the total population count for each individual gene (1:L),
    summing over all genotypes that carry it.

    Returns a length-L vector of counts.

    # Example
        get_margin_gene(G, Freq_host, L)
    """
    function get_margin_gene(G, Freq, L)
        K = length(Freq)
        P = zeros(L)
        for i in 1:K
            P[G[i]] .+= Freq[i]
        end
        return P
    end


    """
        get_margin_gene_array(G, FreqSeries, L)

    Time-series version of `get_margin_gene`.  Given a T×K frequency matrix
    `FreqSeries`, returns a T×L matrix of per-gene counts.
    """
    function get_margin_gene_array(G, FreqSeries, L)
        K = size(G, 1)
        T = size(FreqSeries, 1)
        P = zeros(eltype(FreqSeries), T, L)
        for i in 1:K
            P[:, G[i]] .+= FreqSeries[:, i]
        end
        return P
    end


    # -------------------------------------------------------------------------
    # Recombination: gene-level operation
    # -------------------------------------------------------------------------

    """
        recombination_process(g0, n)

    Replace one randomly chosen gene in genotype `g0` with gene `n`.
    If `n` is already present in `g0`, the genotype is returned unchanged.
    The result is always sorted.

    # Example
        recombination_process([1, 3, 5], 2)  # → [1, 2, 5] or [1, 3, 2] sorted
    """
    function recombination_process(g0, n)
        g = copy(g0)
        if n in g
            return g
        end
        i    = rand(1:length(g))
        g[i] = n
        return sort(g)
    end


    # -------------------------------------------------------------------------
    # Weighted population sampling
    # -------------------------------------------------------------------------

    # Private core shared by both public sampler variants.
    # with_deletion=true  → Freq is decremented in-place (sampling without replacement).
    # with_deletion=false → Freq is read-only; N is clamped to sum(Freq) to prevent errors.
    function _sample_binarysearch(G, Freq, N; with_deletion::Bool)
        total = sum(Freq)
        if total == 0 || N == 0
            return []
        end
        N       = min(N, total)
        cum_sum = cumsum(Freq)
        G_new   = []
        for g in sample(1:last(cum_sum), N, replace=false)
            i = searchsortedfirst(cum_sum, g)
            while Freq[i] == 0
                i += 1
            end
            with_deletion && (Freq[i] -= 1)
            push!(G_new, G[i])
        end
        return G_new
    end

    """
        sample_iterated_binarysearch_with_deletion(G, Freq, N)

    Sample `N` genotypes without replacement, weighted by `Freq`.
    Decrements `Freq` in-place as individuals are removed.
    Uses a binary-search on the cumulative-sum array for O(log K) per draw.
    """
    sample_iterated_binarysearch_with_deletion(G, Freq, N) =
        _sample_binarysearch(G, Freq, N; with_deletion=true)

    """
        sample_iterated_binarysearch_without_deletion(G, Freq, N)

    Sample `N` genotypes without replacement, weighted by `Freq`.
    Does NOT modify `Freq` (read-only sampling). `N` is clamped to `sum(Freq)`.
    """
    sample_iterated_binarysearch_without_deletion(G, Freq, N) =
        _sample_binarysearch(G, Freq, N; with_deletion=false)


    """
        sample_iterated(G, Freq, N)

    Sample `N` genotypes with replacement using a temporary copy of `Freq`
    for within-draw depletion.  Slower than the binary-search variants but
    kept for reference and testing.
    """
    function sample_iterated(G, Freq, N)
        G_new    = []
        Freq_tmp = copy(Freq)
        for _ in 1:N
            gi = sample(1:length(G), Weights(Freq_tmp))
            push!(G_new, G[gi])
            Freq_tmp[gi] -= 1
        end
        return G_new
    end


    # -------------------------------------------------------------------------
    # Recombination mechanisms
    # -------------------------------------------------------------------------

    """
        recombination_random(dict_genenotype2index, G, Freq, L, r, N; del_old=true)

    Random mutation: `floor(r*N)` individuals are replaced by randomly
    drawn genotypes (uniform over all C(L,d) genotypes).

    If `del_old=true`, the source individuals are first removed from `Freq`.

    Returns `(G_new, FreqG_outflow, FreqG_inflow)`.
    """
    function recombination_random(dict_genenotype2index, G, Freq, L, r, N; del_old=true)
        Freq0 = copy(Freq)
        if del_old
            # Modifies Freq in-place; return value intentionally discarded
            sample_iterated_binarysearch_with_deletion(G, Freq, min(sum(Freq), floor(Int64, r * N)))
        end

        FreqG_outflow = Freq0 .- Freq
        Freq0         = copy(Freq)

        sample_G_add = sample(G, floor(Int64, r * N))
        G_new        = add_newgenotypes(dict_genenotype2index, sample_G_add, Freq)

        FreqG_inflow = Freq .- Freq0
        return G_new, FreqG_outflow, FreqG_inflow
    end


    """
        recombination(dict_genenotype2index, G, Freq, L, r, N; del_old=true)

    Pairwise (sexual) recombination: `floor(r*N)` individuals gain a gene
    drawn from the population's marginal gene-frequency distribution.
    One gene is replaced in each recipient genotype.

    Returns `(G_new, FreqG_outflow, FreqG_inflow)`.
    """
    function recombination(dict_genenotype2index, G, Freq, L, r, N; del_old=true)
        Freq0 = copy(Freq)

        if del_old
            sample_G = sample_iterated_binarysearch_with_deletion(G, Freq, min(sum(Freq), floor(Int64, r * N)))
        else
            sample_G = sample_iterated_binarysearch_without_deletion(G, Freq, floor(Int64, r * N))
        end

        FreqG_outflow = Freq0 .- Freq

        # Marginal gene frequencies supply the donor alleles
        P      = get_margin_gene(G, Freq, L)
        gain_g = sample(1:L, Weights(P), floor(Int64, r * N))

        G_sampled = recombination_process.(sample_G, gain_g)

        Freq0        = copy(Freq)
        G_new        = add_newgenotypes(dict_genenotype2index, G_sampled, Freq)
        FreqG_inflow = Freq .- Freq0

        return G_new, FreqG_outflow, FreqG_inflow
    end


    # -------------------------------------------------------------------------
    # Horizontal gene transfer
    # -------------------------------------------------------------------------

    """
        HorizontalGeneTransfer_2body(dict_genenotype2index, G, Freq_host, Freq_phage,
                                     L, r, N; del_old=true)

    2-body HGT: `floor(r*N)` phages each acquire one gene drawn from the
    host marginal gene-frequency distribution (phage infects host).

    Only `Freq_phage` is modified; `Freq_host` is read-only.

    Returns `(G_new, FreqG_outflow, FreqG_inflow)`.
    """
    function HorizontalGeneTransfer_2body(dict_genenotype2index, G, Freq_host, Freq_phage,
                                          L, r, N; del_old=true)
        Freq0 = copy(Freq_phage)

        if del_old
            sample_G = sample_iterated_binarysearch_with_deletion(G, Freq_phage, min(sum(Freq_phage), floor(Int64, r * N)))
        else
            sample_G = sample_iterated_binarysearch_without_deletion(G, Freq_phage, floor(Int64, r * N))
        end

        FreqG_outflow = Freq0 .- Freq_phage

        # Donor genes come from the host population
        P      = get_margin_gene(G, Freq_host, L)
        gain_g = sample(1:L, Weights(P), floor(Int64, r * N))

        G_sampled = recombination_process.(sample_G, gain_g)

        Freq0        = copy(Freq_phage)
        G_new        = add_newgenotypes(dict_genenotype2index, G_sampled, Freq_phage)
        FreqG_inflow = Freq_phage .- Freq0

        return G_new, FreqG_outflow, FreqG_inflow
    end


    """
        HorizontalGeneTransfer_3body_Direct(dict_genenotype2index, G, Freq_host, Freq_phage,
                                            L, r, N)

    Direct 3-body coinfection: two phages and one host are sampled; gene
    transfer occurs only when both phages together carry all host genes
    (full coverage condition).

    Returns `(G_new, G_sampled)`.
    """
    function HorizontalGeneTransfer_3body_Direct(dict_genenotype2index, G, Freq_host, Freq_phage,
                                                  L, r, N)
        sample_phage1 = sample_iterated_binarysearch_without_deletion(G, Freq_phage, floor(Int64, r * N))
        sample_phage2 = sample_iterated_binarysearch_without_deletion(G, Freq_phage, floor(Int64, r * N))
        sample_host   = sample_iterated_binarysearch_without_deletion(G, Freq_host,  floor(Int64, r * N))

        G_sampled = []
        for i in 1:length(sample_phage1)
            append!(G_sampled, recombination_process_3body(sample_host[i], sample_phage1[i], sample_phage2[i]))
        end

        G_new = add_newgenotypes(dict_genenotype2index, G_sampled, Freq_phage)
        return G_new, G_sampled
    end


    """
        HorizontalGeneTransfer_3body_ImportantSampling(dict_genenotype2index, G,
                                                        Freq_host, Freq_phage, L, r, N)

    Importance-sampled 3-body coinfection (more efficient than direct sampling).
    Hosts are weighted by the probability that both coinfecting phages carry
    their genes, using marginal phage gene frequencies.

    Returns `(G_new, G_sampled)`.
    """
    function HorizontalGeneTransfer_3body_ImportantSampling(dict_genenotype2index, G,
                                                             Freq_host, Freq_phage, L, r, N)
        margin_phage      = get_margin_gene(G, Freq_phage, L) ./ sum(Freq_phage)
        Sampling_weights  = ones(length(G))
        for i in 1:length(G)
            for j in G[i]
                Sampling_weights[i] *= margin_phage[j]
            end
        end

        G_sampled = sample_iterated(G, Freq_host .* Sampling_weights, floor(Int64, r * N))
        G_new     = add_newgenotypes(dict_genenotype2index, G_sampled, Freq_phage)

        return G_new, G_sampled
    end


    """
        recombination_process_3body(g_host, g1, g2)

    Determine whether a 3-body coinfection event results in gene transfer.
    Transfer occurs when every host gene is present in at least one of the
    two phage genotypes (`g1` or `g2`).  If so, returns `[g_host]`; otherwise
    returns `[]`.
    """
    function recombination_process_3body(g_host, g1, g2)
        G_new  = []
        infect = true
        for i in g_host
            infect &= (i in g1 || i in g2)
        end
        if infect
            push!(G_new, g_host)
        end
        return G_new
    end


    # -------------------------------------------------------------------------
    # Poisson sampling helper
    # -------------------------------------------------------------------------

    """
        sample_poisson(rate)

    Draw one sample from Poisson(rate).  Returns 0 for non-positive rates.
    """
    function sample_poisson(rate)
        rate <= 0 ? 0 : rand(Poisson(rate))
    end


    # -------------------------------------------------------------------------
    # Fitness functions
    # -------------------------------------------------------------------------

    """
        cal_fitness_no_constraint(Freq_host, Freq_phage, JX, JY, NX, NY)

    Compute per-genotype fitness under unconstrained (exponential) growth.

    Fitness is frequency-dependent via a Lotka-Volterra-like interaction:
    - Host fitness decreases with phage abundance (predation).
    - Phage fitness increases with host abundance (resource).

    Returns `(Fit_host, Fit_phage)`.
    """
    function cal_fitness_no_constraint(Freq_host, Freq_phage, JX, JY, NX, NY)
        K = size(Freq_host, 1)
        return JX .- JX / (NY / K) * Freq_phage,
               JY / (NX / K) * Freq_host .- JY
    end


    """
        cal_fitness_hard_constraint(Freq_host, Freq_phage, JX, JY)

    Compute per-genotype fitness under a hard population-size constraint.
    Fitness is normalised by the mean genotype frequency, making it
    explicitly density-dependent.

    Returns `(Fit_host, Fit_phage)`.
    """
    function cal_fitness_hard_constraint(Freq_host, Freq_phage, JX, JY)
        return JX .- JX / mean(Freq_phage) * Freq_phage,
               JY / mean(Freq_host) * Freq_host .- JY
    end


    # -------------------------------------------------------------------------
    # Selection mechanisms
    # -------------------------------------------------------------------------

    """
        select_no_constraint(Freq, Fit)

    Poisson-sample the next generation without any population-size cap.
    Each genotype i produces `Poisson(Freq[i] * exp(Fit[i]))` offspring.

    Returns `(new_Freq, mean_fitness)`.
    """
    function select_no_constraint(Freq, Fit)
        return sample_poisson.(Freq .* exp.(Fit)), mean(Fit)
    end


    """
        select_soft_constraint(Freq, Fit, Nmax)

    Selection with a soft population-size cap `Nmax`.  When the current
    population `N > Nmax`, an additional logistic penalty is applied so
    that the expected total offspring equals approximately `Nmax`.

    Returns `(new_Freq, effective_mean_fitness)`.
    """
    function select_soft_constraint(Freq, Fit, Nmax)
        N = sum(Freq)
        if N > Nmax
            meanFit = sum(Fit .* Freq) / N          # N already cached above
            N       = sum(Freq .* (exp.(-meanFit .+ Fit)))
            alpha   = log(2)
            rates   = Freq .* (exp.((1 - N / Nmax) * alpha .+ (-meanFit .+ Fit)))
            Freq    = sample_poisson.(rates)
            fit_mean = -meanFit + (1 - N / Nmax) * alpha
        else
            Freq     = sample_poisson.(Freq .* exp.(Fit))
            fit_mean = mean(Fit)
        end
        return Freq, fit_mean
    end


    """
        select_hard_constraint(Freq, Fit, Nmax)

    Selection with a hard population-size constraint.  The expected total
    offspring is forced to equal `Nmax` via a logistic rescaling factor α.

    Returns `(new_Freq, effective_mean_fitness)`.
    """
    function select_hard_constraint(Freq, Fit, Nmax)
        meanFit = sum(Fit .* Freq) / sum(Freq)
        N       = sum(Freq .* (exp.(-meanFit .+ Fit)))
        alpha   = log(2)
        rates   = Freq .* (exp.((1 - N / Nmax) * alpha .+ (-meanFit .+ Fit)))
        Freq    = sample_poisson.(rates)
        return Freq, -meanFit + (1 - N / Nmax) * alpha
    end

end  # module simulation_light
