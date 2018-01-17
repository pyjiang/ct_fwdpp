## set up initial population (scale set to 1)

args=$#
if [ $args != "5" ]; then
  echo "incorrect number of parameters! usage:  bash set_init_pop_for_scale.sh  recomb_rate     stype ga begin end"
  exit
fi


recomb_rate=$1
stype=$2
ga=$3
begin=$4
end=$5


if [ "$stype" = "0" ]; then ## stabil
  general_name="stabil_nz_1000_Lv_10000_ga_${ga}_${recomb_rate}"
else
  general_name="randinit_nz_1000_Lv_10000_ga_${ga}_${recomb_rate}"
fi


log_name="${general_name}_init.log"
echo $log_name
for i in $(eval echo {$begin..$end})
do
  cmd="../ct_haploid_fixed_rob_alpha_init 5000 200 $recomb_rate   $RANDOM 1000 10000  0.5 0.5 $ga $stype ${general_name}_$i 1   &"
  echo $cmd >> $log_name ## save the command to log file

  eval $cmd ## execute the command
  sleep 1
done
