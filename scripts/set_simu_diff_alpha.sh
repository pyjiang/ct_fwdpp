##read from saved state of initial population
## set up simulations with different alpha values

args=$#
if [ $args != "7" ]; then
  echo "incorrect number of parameters! usage:  bash set_simu_diff_alpha.sh   ngen   stype index_begin index_end  init_ga recomb_rate alpha_list"
  exit
fi


ngen=$1
stype=$2
index_begin=$3
index_end=$4
ga=$5
recomb_rate=$6
alpha_list=$(echo $7 | tr "," "\n")  ##note: alpha_list is comma separated list


for index in  $(eval echo {$index_begin..$index_end})
do
    if [ "$stype" = "0" ]; then ## stabil
      general_name="stabil_nz_1000_Lv_10000_ga_${ga}_${recomb_rate}_${index}"
    else
      general_name="randinit_nz_1000_Lv_10000_ga_${ga}_${recomb_rate}_${index}"
    fi

    log_name="${general_name}_alphas.log"
    echo $log_name
    for i in $alpha_list
    do
      cmd="../ct_haploid_fixed_diff_rob_alpha 0 $ngen $RANDOM   ${general_name} $i &"
      echo $cmd >> $log_name ## save the command to log file

      eval $cmd ## execute the command
      sleep 1
    done
done
