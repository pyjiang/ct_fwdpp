## use new version of code to rerun basic version of the simulations (no evolution of robustness)
## for both stabil and randinit
## just to reproduce previous results
## take two parameters: ga, stype(whether stabil or rand)

## number of arguments
args=$#
if [ $args != "2" ]; then
  echo "incorrect number of parameters! usage: bash basic_new.sh gamma stype"
  exit
fi


ga=$1
stype=$2
if [ "$stype" = "0" ]; then ## stabil
  general_name="stabil_nz_1000_Lv_10000_ga_${ga}"
else
  general_name="randinit_nz_1000_Lv_10000_ga_${ga}"
fi

log_name="${general_name}.log"
for i in {1..20}
  do cmd="../ct_haploid_save_state_new 5000 200  16000000  $RANDOM 1000 10000 0.5 0.5 0 $ga $stype ${general_name}_$i &"
  # echo $cmd | tee -a $log_name
  echo $cmd >> $log_name
  eval $cmd
  sleep 1
  done
