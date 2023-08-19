# -*- coding: utf-8 -*-
"""
Created on Thu 10/10/2022
@author: Wenping Cui
example: julia main.jl  --T 200000 --L 40 --J 0.01 --N 10^6 --rR_random 3E-5 --rH_random 3E-5
"""

include("src/simulation_lightning.jl")
using .simulation_light 
using BenchmarkTools, ArgParse, Combinatorics, DataStructures, SpecialFunctions, Printf, StatsBase

"""
merge outputs: awk '!/^(L,)/' *.txt > output.txt

"""

function simulation(dict_genenotype2index, G, Freq_host, Freq_phage, args, T, relaxation_time; stop=false, output_details=true, del_old = true)

	if output_details
    	TimeSeries_genotype_host, TimeSeries_genotype_phage = zeros(Int32, T-relaxation_time, size(G,1)), zeros(Int32, T-relaxation_time, size(G,1))
    end
    
    for t in 1:T
    	NX = sum(Freq_host)
    	if args["rR_random"] != 0
    		G_host_new, FreqG_outflow_host, FreqG_inflow_host = recombination_random(dict_genenotype2index, G, Freq_host, args["L"], args["rR_random"], NX,del_old = del_old) 
        end

        if args["rR_hosthost"] != 0
        	G_host_new , FreqG_outflow_host, FreqG_inflow_host = recombination(dict_genenotype2index, G, Freq_host, args["L"], args["rR_hosthost"], NX, del_old = del_old) 
        end 


        NY = sum(Freq_phage)
        if args["rH_random"] != 0
            G_phage_new = recombination_random(dict_genenotype2index, G, Freq_phage, args["L"], args["rH_random"], NY,del_old = del_old) 
        end
        if args["rH_phagephage"] != 0
            G_phage_new, FreqG_outflow_phage, FreqG_inflow_phage = recombination(dict_genenotype2index, G, Freq_phage, args["L"], args["rH_phagephage"], NY,del_old = del_old) 
        end
        if args["rH_phagehost"] != 0
            G_phage_new, FreqG_outflow_phage, FreqG_inflow_phage = HorizontalGeneTransfer_2body(dict_genenotype2index, G, Freq_host, Freq_phage, args["L"], args["rH_phagehost"], NY,del_old = del_old)
        end
        if args["rH_3body"] != 0
            G_phage_new, G_phage_sampled = HorizontalGeneTransfer_3body_ImportantSampling(dict_genenotype2index, G, Freq_host, Freq_phage, args["L"], args["rH_3body"], NY)
        end

        
        if args["population_constraint"] =="hard"
			Fit_host, Fit_phage = cal_fitness_hard_constraint(Freq_host, Freq_phage, args["J"], args["Jprop"]*args["J"])
			Freq_host, meanfit_host = select_hard_constraint(Freq_host, Fit_host, args["NX"])
			Freq_phage, meanfit_phage = select_hard_constraint(Freq_phage, Fit_phage, args["NY"])

		elseif args["population_constraint"] =="soft"
			Fit_host, Fit_phage = cal_fitness_no_constraint(Freq_host, Freq_phage, args["J"], args["J"]*args["Jprop"], args["NX"], args["NY"])
        	Freq_host, meanfit_host = select_soft_constraint(Freq_host, Fit_host, args["soft_threshold"]*args["NX"])
        
        	Freq_phage, meanfit_phage = select_soft_constraint(Freq_phage, Fit_phage, args["soft_threshold"]*args["NY"])

		elseif args["population_constraint"] =="no"
			Fit_host, Fit_phage = cal_fitness_no_constraint(Freq_host, Freq_phage, args["J"], args["J"]*args["Jprop"], args["NX"], args["NY"])
        	Freq_host, meanfit_host = select_no_constraint(Freq_host, Fit_host)
        	Freq_phage, meanfit_phage = select_no_constraint(Freq_phage, Fit_phage)
        end

        margin_host, margin_phage = get_margin_gene(G, Freq_host, L), get_margin_gene(G, Freq_phage, L)

        if output_details && t-relaxation_time>0
        	TimeSeries_genotype_host[t-relaxation_time,:]=Freq_host
        	TimeSeries_genotype_phage[t-relaxation_time,:]=Freq_phage
	    end

        if stop && (sum(margin_host .==0.)>0 || sum(margin_phage .> 0)<3 )
        	if output_details
        		return TimeSeries_genotype_host[1:t-1,:], TimeSeries_genotype_phage[1:t-1,:]
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

function AnalyzeData(args, TimeSeries_genotype_host, TimeSeries_genotype_phage, G, L, save_dir, relaxation_time)

	T = size(TimeSeries_genotype_host, 1) + relaxation_time
	# established host population size
	mean_Host_N, mean_Phage_N = mean(sum(TimeSeries_genotype_host,dims=2)), mean(sum(TimeSeries_genotype_phage,dims=2))
	std_Host_N, std_Phage_N = std(sum(TimeSeries_genotype_host,dims=2)), std(sum(TimeSeries_genotype_phage,dims=2))

	estab_size = find_establishment_size(TimeSeries_genotype_host, TimeSeries_genotype_phage)
	estab_size_mean = isempty(estab_size) ? 0 : mean(estab_size) 

	# genotype persistence time and 
	persTimeHost, gapTimeHost = genotypePers(TimeSeries_genotype_host)
	persTimePhage, gapTimePhage = genotypePers(TimeSeries_genotype_phage)

	persTimeHost_mean = isempty(persTimeHost) ? 0 : mean(persTimeHost) 
	persTimePhage_mean = isempty(persTimePhage) ? 0 : mean(persTimePhage) 
	gapTimeHost_mean = isempty(gapTimeHost) ? 0 : mean(gapTimeHost) 
	gapTimePhage_mean = isempty(gapTimePhage) ? 0 : mean(gapTimePhage) 


	# calculatate theta absolute
	theta_GH = cal_theta(TimeSeries_genotype_host)
    
    TimeSeries_gene_host = get_margin_gene_array(G, TimeSeries_genotype_host, L)

    TimeSeries_genotype_host = 0

    theta_gH = cal_theta(TimeSeries_gene_host)



	theta_GP = cal_theta(TimeSeries_genotype_phage)

	TimeSeries_gene_phage = get_margin_gene_array(G, TimeSeries_genotype_phage, L)

	TimeSeries_genotype_phage = 0

	theta_gP = cal_theta(TimeSeries_gene_phage)



	save_data = [args["L"], args["d"], args["J"], args["Jprop"], args["NX"], args["NY"], args["rR_random"], args["rH_random"],args["rR_hosthost"], args["rH_phagephage"], args["rH_phagehost"], args["rH_3body"], args["T"],
	T, estab_size_mean, persTimeHost_mean, persTimePhage_mean, gapTimeHost_mean, gapTimePhage_mean, theta_GH, theta_GP, theta_gH, theta_gP, args["population_constraint"], args["soft_threshold"], mean_Host_N, std_Host_N, mean_Phage_N ,std_Phage_N]

	column_name = ["L", "d", "J", "Jprop", "NX", "NY", "rR_random", "rH_random", "rR_hosthost", "rH_phagephage", "rH_phagehost", "rH_3body", "simulation_time", "GenePersisentTime", "mean_establish_size", "persTimeHost", "persTimePhage", "gapTimeHost", "gapTimePhage", "theta_GH", "theta_GP", "theta_gH", "theta_gP", "constraint", "soft_threshold", "mean_host_N", "std_host_N", "mean_phage_N", "std_phage_N"]
	file_name = "L_"*string(args["L"])*"_N_"*@sprintf("%.2E",args["NX"])*"_J_"*string(args["J"])*"_outputTask"*string(args["taskID"])*"Thread"*string(args["threadID"])*".txt"
	write2file(joinpath(save_dir,file_name), save_data, column_name)
end

function main(pars; relaxation_time = 5*10^3)
	args = OrderedDict()
	L, d, rR_random, rH_random, rR_hosthost, rH_phagephage, rH_phagehost, rH_3body, N, J, Jprop, Nprop_phagedivhost, T, taskID, threadID, trialID, output_details, population_constraint, soft_threshold, save_dir= pars
	#L, r, N, J = 15, 1. *10^-5, 10^8, 5.

	L = Int64(L)
	N = Int64(N)
	T = Int64(T)
	trialID = Int64(trialID)
	args["L"] = L
	args["d"] = d
	args["rR_random"] = rR_random
	args["rR_hosthost"] = rR_hosthost
	args["rH_random"] = rH_random
	args["rH_phagephage"] = rH_phagephage
	args["rH_phagehost"] = rH_phagehost
	args["rH_3body"] = rH_3body
	args["NX"] = N
	args["NY"] = N*Nprop_phagedivhost
	args["J"] = J
	args["Jprop"] = Jprop
	args["T"] = T
	args["taskID"], args["threadID"], args["trialID"]  = taskID, threadID, trialID
	args["population_constraint"] = population_constraint
	args["soft_threshold"] = soft_threshold
	G, dict_genenotype2index, dict_index2genenotype, Freq_host, Freq_phage = initial_genotypes(args["L"], args["NX"], args["NY"], args["d"] )
	
	t = 0
	if output_details
		TimeSeries_genotype_host, TimeSeries_genotype_phage = simulation(dict_genenotype2index, G, Freq_host, Freq_phage, args, T, relaxation_time; stop=true, output_details=true, del_old = true)
		t = size(TimeSeries_genotype_host, 1)+relaxation_time
		AnalyzeData(args, TimeSeries_genotype_host, TimeSeries_genotype_phage, G, L, save_dir, relaxation_time)
	else
		t = simulation(dict_genenotype2index, G, Freq_host, Freq_phage, args, T, relaxation_time; stop=true, output_details=false, del_old = true)
		column_name = ["L", "d", "J", "Jprop", "NX", "NY", "rR_random", "rH_random", "rR_hosthost", "rH_phagephage", "rH_phagehost", "rH_3body","GenePersisentTime", "constraint", "soft_threshold"]
		save_data = [args["L"], args["d"], args["J"], args["Jprop"], args["NX"], args["NY"], args["rR_random"], args["rH_random"], args["rR_hosthost"], args["rH_phagephage"], args["rH_phagehost"], args["rH_3body"], t, args["population_constraint"], args["soft_threshold"]]
		file_name = "L_"*string(args["L"])*"_N_"*@sprintf("%.2E",args["NX"])*"_J_"*string(args["J"])*"_outputTask"*string(args["taskID"])*"Thread"*string(args["threadID"])*"_geneduration.txt"
		write2file(joinpath(save_dir,file_name), save_data, column_name)
	end

	println(args, "\n")
	println(t, "\n")
	
	return 
end 


function parse_commandline()
	s = ArgParseSettings()

	@add_arg_table s begin
		"--L"
			help = "total number of genes"
			arg_type = Int64
			default = 30
		"--d"
			help = "number of genes in each genotype"
			arg_type = Int64
			default = 2
		"--r"
			help = "host and phage recombination rate"
			arg_type = Float64
			default = 10^(-3)
		"--J"
			help = "interaction strength"
			arg_type = Float64
			default = 5.
		"--Jprop"
			help = "ratio between phage and host interaction strength"
			arg_type = Float64
			default = 1.
		"--N"
			help = "population size"
			arg_type = Float64
			default = 1. *10^6
		"--T"
			help = "simulation generations"
			arg_type = Float64
			default = 5. *10^5
		"--taskID"
			help = "task ID"
			arg_type = Int64
			default = 1
		"--threadID"
			help = "thread ID"
			arg_type = String
			default = "1"
		"--trialID"
			help = "trial ID"
			arg_type = Int64
			default = 1
		"--rR_hosthost"
			help = "host-host recombination rate"
			arg_type = Float64
			default = 0.
		"--rR_random"
			help = "host random recombination rate"
			arg_type = Float64
			default = 0.
		"--rH_random"
			help = "phage random recombination rate"
			arg_type = Float64
			default = 0.
		"--rH_phagephage"
			help = "phage-phage recombination rate"
			arg_type = Float64
			default = 0.
		"--rH_phagehost"
			help = "phage-host horizontal gene transfer rate"
			arg_type = Float64
			default = 0.
		"--rH_3body"
			help = "phage coinfection recombination rate"
			arg_type = Float64
			default = 0.
		"--Nprop_phagedivhost"
			help = "ratio between phage and host abundance"
			arg_type = Float64
			default = 1.
		"--r_prop"
			help = "ratio between phage and host recombination rate"
			arg_type = Float64
			default = 1.
		"--output_details"
			help = "save statistical properties"
			arg_type = Bool
			default = true
		"--population_constraint"
			help = "soft, hard, no"
			arg_type = String
			default = "soft"
		"--save_dir"
			help = "save data folder"
			arg_type = String
			default = "results_simulation_test"
		"--soft_threshold"
			help = "maximum number of population size"
			arg_type = Float64
			default = 5.
	end
	return parse_args(s)
end

#L, r, N, J, T = 120, 1*10^-4, 10^6, 5, 2*10^5

jobs = []

parsed_args = parse_commandline()


L, r, N, J ,T = parsed_args["L"], parsed_args["r"], parsed_args["N"], parsed_args["J"], floor(Int64, parsed_args["T"])
rR_random, rH_random, rR_hosthost, rH_phagephage, rH_phagehost, rH_3body =parsed_args["rR_random"],  parsed_args["rH_random"], parsed_args["rR_hosthost"], parsed_args["rH_phagephage"], parsed_args["rH_phagehost"], parsed_args["rH_3body"]
Jprop = parsed_args["Jprop"]
Nprop_phagedivhost  = parsed_args["Nprop_phagedivhost"]
d = parsed_args["d"]
r_prop = parsed_args["r_prop"]
rH_random, rH_phagephage, rH_phagehost, rH_3body = r_prop*rH_random, r_prop*rH_phagephage, r_prop*rH_phagehost, r_prop*rH_3body


isdir(parsed_args["save_dir"]) || mkdir(parsed_args["save_dir"])

taskID =  parsed_args["taskID"]
#taskID = parse(Int64, parsed_args["taskID"])
threadID = parse(Int64, parsed_args["threadID"])
trialID = parsed_args["trialID"]

println(threadID)
output_details = parsed_args["output_details"]
@time main((L, d, rR_random, rH_random, rR_hosthost, rH_phagephage, rH_phagehost, rH_3body, N, J, Jprop, Nprop_phagedivhost,T, taskID, threadID, trialID, output_details, parsed_args["population_constraint"],  parsed_args["soft_threshold"], parsed_args["save_dir"]))
