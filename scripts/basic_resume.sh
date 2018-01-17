## haploid resume
## resume running simulations
## bash hap_resume.sh  common_name ga begin end prev_ngen ngens (begin, end are simulation indexes)

ga=$2
for i in $(eval echo {$3..$4})
do
  common_name="$1_ga_${ga}_$i"
  ../ct_haploid_resume_new ${common_name}  $5 $6 1  &
done
