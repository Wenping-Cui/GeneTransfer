#!/bin/bash -l


SGE_TASK_ID=1

rR_array=()
for i in $(seq -6. 0.2 -1.5); do 
    output=$(bc -l <<<"e($i*l(10))")
    rR_array+=($output)
done

rR_array=( $(shuf -e ${rR_array[@]}) )

Ninput=(1000000 10000000)

Linput=$(seq 30 10 50)

Jinput=(0.005)



parallel --jobs 28 main.jl  --T 250000 --taskID $SGE_TASK_ID --threadID {%} --d 2 --J {1} --N {2} --rR_hosthost {3} --rH_phagehost {3} --L {4} --save_dir "results/"  --population_constraint "hard" --output_details true ::: ${Jinput[@]}   ::: ${Ninput[@]} ::: ${rR_array[@]} ::: ${Linput[@]} 