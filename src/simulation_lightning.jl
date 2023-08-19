# -*- coding: utf-8 -*-
"""
Created on Thu 10/10/2022
@author: Wenping Cui
"""

module simulation_light

	using Distributions, StatsBase, Combinatorics
	
	export initial_genotypes, recombination, HorizontalGeneTransfer_2body, HorizontalGeneTransfer_3body_Direct, HorizontalGeneTransfer_3body_ImportantSampling, get_margin_gene, get_margin_gene_array, sample_iterated_3d, add_newgenotypes, find_establishment_size, genotypePers, write2file, recombination_random, sample_poisson,cal_fitness_no_constraint, select_no_constraint, select_soft_constraint, cal_fitness_hard_constraint, select_hard_constraint, cal_theta, vecvec2matrix

	function write2file(filename, save_data, column_name)
		if !isfile(filename)
			open(filename,"w")  do io
			    print(io, column_name[1])
			    for d in column_name[2:end]
			        print(io, ",", d)
			    end
			    print(io, "\n")
			end
		end

		open(filename,"a")  do io
		    print(io, save_data[1])
		    for d in save_data[2:end]
		        print(io, ",", d)
		    end
		    print(io, "\n")
		end
	end

	function vecvec2matrix(M)
	    return reduce(vcat,transpose.(M))
	end


	function cal_theta(TimeSeries; keepzeros = true)
	    data = vec(TimeSeries[1:end, :])
	    if keepzeros
	    	data = data
	    else
	    	data = data[data .>0]
	    end
	    return mean((vec(data) .- mean(data)).^2 ./mean(data))
	end


	function initial_genotypes(L, NX, NY, d)
	    """
	    function: Initilize genotypes carrying d different genes 
	    Input:  
	        L Int, # of genes
	        d Int, # of genes one strain can carry
	        NX, NY, Int, total population size of host and phages
	    Output: 
	        G: vector of vector, representing different genotypes, its length is indentical to Freq

	        Freq_host, Freq_phage, population size for each genotypes

	        dict_genenotype2index: a hashmap: 
	        map genotypes to index in Freq_host and Freq_phage
	    
	        dict_index2genenotype: a hashmap: 
	        map freq index to genotypes
	        
	    test example: 
	        L, NX, NY, d = 10,  10^6, 10^6, 2
	        G, dict_genenotype2index, dict_index2genenotype, Freq_host, Freq_phage = initial_genotypes(L, NX, NY, d)
	    """
	    dict_genenotype2index = Dict{Array{Int64,1}, Int64}()
	    dict_index2genenotype = Dict{Int64, Array{Int64,1}}()
	    G = collect(combinations(1:L, d) )
	    K = size(G, 1)
	    for i in 1:K 
	        dict_genenotype2index[G[i]] = i
	        dict_index2genenotype[i] = G[i]
	    end
	    G = collect(combinations(1:L, d) )
	    Freq_host, Freq_phage =  floor(Int64, NX/K).*vec(ones(Int64, 1,K)), floor(Int64, NY/K).*vec(ones(Int64, 1,K))
	    return G, dict_genenotype2index, dict_index2genenotype, Freq_host, Freq_phage
	end

	function find_establishment_size(TimeSeries_genotype_host, TimeSeries_genotype_phage)
	    est_size= []
	    ave = 1.4*mean(TimeSeries_genotype_host)
	    for i in 1:size(TimeSeries_genotype_host, 2)
	        dataX = TimeSeries_genotype_host[:, i]
	        dataY = TimeSeries_genotype_phage[:, i]
	        birthtime = 1
	        avePass = false
	        birth = false
	        for t in 2:length(dataY)-1
	            if dataY[t-1] <= 0.1 && dataY[t] > 0 #birth
	                birthtime = t
	                birth = true
	            end

	            if !avePass && birth && dataY[t] >= ave 
	                avePass = true      
	            end

	            if dataY[t] == 0. #death
	                if birth && avePass
	                    if birthtime!=1
	                        push!(est_size, dataX[birthtime])
	                    end
	                end
	            avePass = false  
	            birth = false      
	            end
	        end
	    end
	    return est_size 
	end

	function genotypePers(TS)
	    persTime = [] 
	    gapTime = []
	    ave = mean(TS)
	    for i in 1:size(TS, 2)
	        data = TS[:, i]
	        birthtime = 1
	        deathtime = 1
	        avePass = false
	        birth = false
	        count_persTime = 0
	        for t in 2:size(TS, 1)
	            if data[t-1] <= 0.1 && data[t] > 0 #birth
	                birthtime = t
	                birth = true
	            end

	            if !avePass && birth && data[t] >= ave 
	                avePass = true      
	            end

	            if data[t] == 0. #death
	                if birth && avePass
	                    if deathtime!=1
	                        push!(gapTime, birthtime - deathtime)
	                    end
	                    count_persTime += 1
	                    push!(persTime, t-birthtime)
	                    deathtime = t
	                end
	            avePass = false  
	            birth = false      
	            end
	        end
	        if count_persTime == 0
	        	push!(persTime, size(TS, 1))
	        end
	    end
	    return persTime, gapTime
	end

	function count_genotypes(dict_genenotype2index, Freq_like, G_sampled)
	        Freq = zero(Freq_like)
	        for (g, n) in countmap(G_sampled)
	            Freq[dict_genenotype2index[g]] = n
	        end
	        return Freq
	end

	function recombination_random(dict_genenotype2index, G, Freq, L, r, N; del_old = true)
	    """
	    dict_genenotype2index: hashmap from genotype to index in Freq vector
	    
	    G: genotypes, updated in this function
	    
	    Freq: frequency of each genotype, update in this function
	    
	    r: recombination rate
	    
	    N: total number of species in the system
	    """
	    Freq0 = copy(Freq)
	    # random mutation from existed genotype to arbitary genotypes
	    if del_old
	    	sample_G_del = sample_iterated_binarysearch_with_deletion(G, Freq, min(sum(Freq), floor(Int64, r*N)) )
	    end
	 
	 	FreqG_outflow = Freq0 .- Freq

	 	Freq0 = copy(Freq)

	    sample_G_add = sample(G, floor(Int64, r*N))
	 
	    # add new genotypes
	    G_new = add_newgenotypes(dict_genenotype2index, sample_G_add, Freq);
	    
	    FreqG_inflow = Freq .- Freq0

	    return G_new, FreqG_outflow, FreqG_inflow
	    
	end

	function sample_iterated_binarysearch_with_deletion(G, Freq, N)
	    """
	    function: sample N genotypes from the pool weighted by its population size

	    Input: 
	        G, array of vector, represnnting different genotypes
	        Freq, array of Int, representing population size, its length is identical to G.
	        N, Int sampling size
	    Output:
	        G_new, array of vector with length N
	        Freq, the updated population frequency

	    test case: sample_G, Freq_host = sample_iterated(G, Freq_host, floor(Int64, r*N));
	    """
	    if sum(Freq)==0
	    	return []
	    end
	    cum_sum  = cumsum(Freq)
	    G_new = []
	    for g in sample(1:last(cum_sum), N, replace = false)
	        i = searchsortedfirst(cum_sum, g)
	        while Freq[i] == 0
	             i += 1
	        end
	        Freq[i] -= 1
	        push!(G_new, G[i])
	    end

	    return G_new
	end

	function sample_iterated_binarysearch_without_deletion(G, Freq, N)
	    """
	    function: sample N genotypes from the pool weighted by its population size

	    Input: 
	        G, array of vector, represnnting different genotypes
	        Freq, array of Int, representing population size, its length is identical to G.
	        N, Int sampling size
	    Output:
	        G_new, array of vector with length N
	        Freq, the updated population frequency

	    test case: sample_G, Freq_host = sample_iterated(G, Freq_host, floor(Int64, r*N));
	    """
	    cum_sum  = cumsum(Freq)
	    G_new = []
	    for g in sample(1:last(cum_sum), N, replace = false)
	        i = searchsortedfirst(cum_sum, g)
	        while Freq[i] == 0
	             i += 1
	        end
	        push!(G_new, G[i])
	    end

	    return G_new
	end

	function recombination(dict_genenotype2index, G, Freq, L, r, N; del_old = true)
	    """
	    dict_genenotype2index: hashmap from genotype to index in Freq vector
	    
	    G: genotypes, updated in this function
	    
	    Freq: frequency of each genotype, update in this function
	    
	    r: recombination rate
	    
	    N: total number of species in the system
	    """
	    
	    # sample hosts who are going to gain genes, 
	    # minus hosts elminated by recombination, and update the new population size

	    Freq0 = copy(Freq)

	    if del_old
	    	#sample_G = sample_iterated_deletion(G, Freq, floor(Int64, r*N))
	    	sample_G = sample_iterated_binarysearch_with_deletion(G, Freq, min(sum(Freq), floor(Int64, r*N)))
	    else 
	    	#sample_G = sample_iterated(G, Freq, floor(Int64, r*N))
	    	sample_G = sample_iterated_binarysearch_without_deletion(G, Freq, floor(Int64, r*N))
	    end
	    
	    FreqG_outflow = Freq0 .- Freq

	    # get marginal probability of pure genes
	    P =  get_margin_gene(G, Freq, L)
	    
	    # sample pure genes which hosts could gain
	    gain_g = sample(1:L, Weights(P), floor(Int64, r*N));
	    
	    # do recombination and gain new genotypes
	    G_sampled = recombination_process.(sample_G, gain_g);

	    #FreqG_inflow = count_genotypes(dict_genenotype2index, Freq, G_sampled)

	    Freq0 = copy(Freq)
	    
	    # add new genotypes
	    G_new = add_newgenotypes(dict_genenotype2index, G_sampled, Freq);
	    
	    FreqG_inflow = Freq .- Freq0

	    return G_new, FreqG_outflow, FreqG_inflow
	    
	end

	function sample_iterated_deletion(G, Freq, N)
	    """
	    function: sample N genotypes from the pool weighted by its population size
	    
	    Input: 
	        G, array of vector, represnnting different genotypes
	        Freq, array of Int, representing population size, its length is identical to G.
	        N, Int sampling size
	    Output:
	        G_new, array of vector with length N
	        Freq, the updated population frequency
	    
	    test case: sample_G, Freq_host = sample_iterated(G, Freq_host, floor(Int64, r*N));
	    """
	    G_new = []
	    
	    for i in 1:N
	        gi = sample(1:length(G), Weights(Freq))
	        push!(G_new, G[gi])
	        Freq[gi] -= 1 # avoid the number excesses the strain population in the current pool.
	    end
	    
	    return G_new
	end



	function get_margin_gene(G, Freq, L)
	    """
	    function: calculate the population size of each gene
	    
	    test case: get_margin_gene(G, Freq_host, L)
	    """
	    K = length(Freq)
	    P = zeros(L)
	    for i in 1:K
	        P[G[i]] .+= Freq[i]
	    end
	    return P
	end

	function get_margin_gene_array(G, FreqSeries, L)
	    """
	    function: calculate the population size of each gene
	    """
	    K = size(G, 1)
	    T = size(FreqSeries, 1)
	    P = zeros(eltype(FreqSeries), T, L)
	    for i in 1:K
	        P[:, G[i]] .+= FreqSeries[:, i]
	    end
	    return P
	end

	function recombination_process(g0, n)
	    """
	    g: genotype
	    
	    n: gained pure gene
	    
	    replace the gene in g with n randomly
	    
	    test case: recombination_process(G[1], 5)
	    """
	    g = copy(g0)
	    if n in g
	        return g
	    end
	    i = rand(1:length(g))
	    g[i] = n
	    return sort(g)
	end

	function add_newgenotypes(dict_genenotype2index, G_sampled, Freq)
	     """
	    if existed in the pool 
	        the corresponding frequency +1
	    
	    else frequency==0, means it does not existed in the pool
	        insert the genotype in the pool  
	        and return the new genotype not sampled before
	    """
	    G_new = []
	    
	    for g in G_sampled
	        i = dict_genenotype2index[g]
	        if Freq[i]==0
	            push!(G_new, g)
	        end 
	        Freq[i] += 1
	    end
	    
	    return G_new
	end

	function HorizontalGeneTransfer_2body(dict_genenotype2index, G, Freq_host, Freq_phage, L, r, N; del_old = true)
	    """
	    dict_genenotype2index: hashmap from genotype to index in Freq vector
	    
	    G: genotypes
	    
	    Freq_host, Freq_phage: frequency of each genotype, update in this function
	    
	    r: horizontal gene transfer rate
	    
	    N: total number of species in the system
	    """
	    Freq0 = copy(Freq_phage)

	    # sample phages who are going to gain genes
	    if del_old
	    	#sample_G = sample_iterated_deletion(G, Freq, floor(Int64, r*N))
	    	sample_G = sample_iterated_binarysearch_with_deletion(G, Freq_phage, min(sum(Freq_phage), floor(Int64, r*N)))
	    else 
	    	#sample_G = sample_iterated(G, Freq, floor(Int64, r*N))
	    	sample_G = sample_iterated_binarysearch_without_deletion(G, Freq_phage, floor(Int64, r*N))
	    end

	    #FreqG_outflow = count_genotypes(dict_genenotype2index, Freq_host, sample_G)
		FreqG_outflow = Freq0 .- Freq_phage  


	    # get marginal probability of pure genes from the hosts
	    P =  get_margin_gene(G, Freq_host, L)
	    
	    # sample pure genes which hosts could gain
	    gain_g = sample(1:L, Weights(P), floor(Int64, r*N))
	    
	    # do recombination and gain new genotypes
	    G_sampled = recombination_process.(sample_G, gain_g);
	    
	    #FreqG_inflow = count_genotypes(dict_genenotype2index, Freq_host, G_sampled)

	    # add new genotypes
	    Freq0 = copy(Freq_phage)

	    G_new = add_newgenotypes(dict_genenotype2index, G_sampled, Freq_phage);
	    
	    FreqG_inflow = Freq_phage .- Freq0

	    return G_new, FreqG_outflow, FreqG_inflow
	end

	function sample_iterated(G, Freq, N)
	    """
	    function: sample N genotypes from the pool weighted by its population size
	    
	    Input: 
	        G, array of vector, represnnting different genotypes
	        Freq, array of Int, representing population size, its length is identical to G.
	        N, Int sampling size
	    Output:
	        G_new, array of vector with length N
	        Freq, the updated population frequency
	    
	    test case: sample_G, Freq_host = sample_iterated(G, Freq_host, floor(Int64, r*N));
	    """
	    G_new = []
	    Freq_tmp = copy(Freq)
	    for i in 1:N
	        gi = sample(1:length(G), Weights(Freq_tmp))
	        push!(G_new, G[gi])
	        Freq_tmp[gi] -= 1 # avoid the number excesses the strain population in the current pool.
	    end
	    
	    return G_new
	end

	function delete_genotypes(dict_genenotype2index, G_sampled, Freq)
	     """
	    delete genotypes from the pool
	    
	    update the frequency of the population size
	    """
	    G_new = []
	    
	    for g in G_sampled
	        i = dict_genenotype2index[g]
			if Freq[i] > 0 
	        	Freq[i] -= 1
			end
	    end
	    
	    return G_new
	end

	function HorizontalGeneTransfer_3body_Direct(dict_genenotype2index, G, Freq_host, Freq_phage, L, r, N)
	    """
	    dict_genenotype2index: hashmap from genotype to index in Freq vector
	    
	    G: genotypes
	    
	    Freq_host, Freq_phage: frequency of each genotype, update in this function
	    
	    r: horizontal gene transfer rate
	    
	    N: total number of species in the system
	    """
	    # sample phages who are going to gain genes
	    sample_phage1 = sample_iterated_binarysearch_without_deletion(G, Freq_phage, floor(Int64, r*N))
	    
	    sample_phage2 = sample_iterated_binarysearch_without_deletion(G, Freq_phage, floor(Int64, r*N))
	    
	    sample_host = sample_iterated_binarysearch_without_deletion(G, Freq_host, floor(Int64, r*N))
	    
	    # do recombination and gain new genotypes
	    G_sampled = []
	    for i in 1:length(sample_phage1)
	        append!(G_sampled, recombination_process_3body(sample_host[i], sample_phage1[i], sample_phage2[i]))
	    end 
	    # add new genotypes
	    G_new = add_newgenotypes(dict_genenotype2index, G_sampled, Freq_phage);
	    
	    return G_new, G_sampled
	end

	function HorizontalGeneTransfer_3body_ImportantSampling(dict_genenotype2index, G, Freq_host, Freq_phage, L, r, N)
	    """
	    dict_genenotype2index: hashmap from genotype to index in Freq vector

	    G: genotypes

	    Freq_host, Freq_phage: frequency of each genotype, update in this function

	    r: horizontal gene transfer rate

	    N: total number of species in the system
	    """
	    # sample phages who are going to gain genes
	    margin_phage = get_margin_gene(G, Freq_phage, L)./sum(Freq_phage)
	    Sampling_weights = ones(length(G)) 
	    for i in 1:length(G)
	        p = 1
	        for j in G[i]
	            Sampling_weights[i] *= margin_phage[j]
	        end
	    end
	    # do recombination and gain new genotypes
	    G_sampled = sample_iterated(G, Freq_host.*Sampling_weights, floor(Int64, r*N))
	    # add new genotypes
	    G_new = add_newgenotypes(dict_genenotype2index, G_sampled, Freq_phage);

	    return G_new, G_sampled
	end


	function recombination_process_3body(g_host, g1, g2)
	    G_new = []
	    infect = true
	    for i in g_host
			infect &= (i in g1 || i in g2)
	    end
	    if infect
	    	push!(G_new, g_host)
	        # for i in g1
	        #     for j in g2
	        #         if i != j
	        #             push!(G_new, sort([i, j]))
	        #         end
	        #     end
	        # end
	    end  
	    return G_new
	end

	function sample_poisson(rate)
	    # sample a poisson distribution
	    if rate <= 0
	    	return 0
	    else 
	    	return rand(Poisson(rate))
	    end
	end



	function cal_fitness_no_constraint(Freq_host, Freq_phage, JX, JY, NX, NY)
	    """
	    calculate the fitness

	    Freq_host, Freq_phage : frequency of each genotype
	    """
	    return JX.-JX/(NY/size(Freq_host, 1))*Freq_phage, JY/(NX/size(Freq_host, 1))*Freq_host .- JY
	end


	function select_no_constraint(Freq, Fit)
	    """
	    selections without constraint

	    Freq_host, Freq_phage : frequency of each genotype
	    """
	    return sample_poisson.(Freq .*exp.(Fit)), mean(Fit)
	end


	function select_soft_constraint(Freq, Fit, Nmax)
	    """
	    selections with soft constraint

	    Freq_host, Freq_phage : frequency of each genotype

	    Nmax: the maxium 

	    """
	    N = sum(Freq)
	    if N>Nmax
	    	meanFit = sum(Fit.*Freq)/sum(Freq)
	    	N = sum(Freq.*(exp.(-meanFit.+Fit)))
	        alpha = log(2)
	        rates = Freq.*(exp.((1-N/Nmax)*alpha.+(-meanFit.+Fit)))
	        Freq = sample_poisson.(rates)
	        fit_mean = -meanFit+(1-N/Nmax)*alpha
	    else
	        Freq = sample_poisson.(Freq .*exp.(Fit))
	        fit_mean = mean(Fit)
	    end
	    return Freq, fit_mean
	end


	function cal_fitness_hard_constraint(Freq_host, Freq_phage, JX, JY)
	    """
	    calculate the fitness

	    Freq_host, Freq_phage : frequency of each genotype
	    """
	    return JX .-JX/mean(Freq_phage)*Freq_phage, JY/mean(Freq_host)*Freq_host .-JY
	end


	function select_hard_constraint(Freq, Fit, Nmax)
	    """
	    perform selection process
	    
	    Update genotype frequency according to 
	    poisson sampling results
	    """
	    meanFit = sum(Fit.*Freq)/sum(Freq)
	    N = sum(Freq.*(exp.(-meanFit.+Fit)))
	    alpha = log(2)
	    rates = Freq.*(exp.((1-N/Nmax)*alpha.+(-meanFit.+Fit)))
	    Freq = sample_poisson.(rates)
	    return Freq, -meanFit+(1-N/Nmax)*alpha
	end
end