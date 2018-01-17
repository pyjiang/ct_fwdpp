## modifed on 8.8.17, add ourput dir
## modified on 3.6.17
## added robust bits
## modified on 2.3.17
## do not pass scale_func_flag but pass scale_flag as parameter
## modified on 2.2.17
## add mu_rob into file_names
## currently only use 3 degrees of initial robustness
## see behavior with and without evolve robustness
## use new version of code to rerun basic version of the simulations (no evolution of robustness)
## for both stabil and randinit
## just to reproduce previous results
## take 4 parameters: stype(whether stabil or rand), robust_flag (robust_flag=2: control, robust_flag=1: additve robustness) low_scale high_scale

## number of arguments
args=$#
if [ $args != "10" ]; then
  echo "incorrect number of parameters! usage: bash set_init_diff_rob_evol.sh ngen stype ga robust_flag robust_bits shuffle_ratio mu_rob low_bound high_bound  scale_flag"
  exit
fi

ngen=$1
stype=$2
ga=$3
rob_flag=$4
robust_bits=$5
shuffle_ratio=$6
mu_rob=$7 ## mutation rate per robust bit
low_bound=$8
high_bound=$9
scale_flag=${10}



if [ "$rob_flag" = "1" ]; then
  rob="_1"
elif [ "$rob_flag" = "2" ]; then
  rob="_ctrl"
elif [ "$rob_flag" = "0" ]; then
  rob="_0"
fi


if [ "$stype" = "0" ]; then ## stabil
  general_name="stabil_${low_bound}_${high_bound}_sf${scale_flag}_${shuffle_ratio}_nz_1000_Lv_10000_ga_${ga}${rob}_${mu_rob}"
else
  general_name="randinit_${low_bound}_${high_bound}_sf${scale_flag}_${shuffle_ratio}_nz_1000_Lv_10000_ga_${ga}${rob}_${mu_rob}"
fi


log_name="${general_name}.log"
echo $log_name
for i in {1..20}
do
  cmd="../ct_haploid_diff_init_rob_evol 5000 200 $shuffle_ratio $ngen  $RANDOM $mu_rob 1000 10000 $robust_bits 0.5 0.5  $ga  $stype  ${general_name}_$i   $rob_flag   $low_bound $high_bound $scale_flag &"
  echo $cmd >> $log_name ## save the command to log file

  eval $cmd ## execute the command
  sleep 1
done
