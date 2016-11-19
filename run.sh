#!/bin/sh

mkdir "clade1234_marginal"
./pro_4comp.o -m kinship_impute_130.txt -p kinship_clade1234.txt -c Cov_miaoyan.txt -l country_member_130_1.txt -n ncount.txt -g sequence_130_1.txt -s SNPMAF_clade1234_unique_pattern_33610.txt -o "result_acc.txt" -b "trajectory_acc.txt" -r "chisq_bac.txt" -t "chisq_acc.txt"
#done


##mv "result_acc.txt" ./clade1234_marginal
##mv "trajectory_acc.txt" ./clade1234_marginal
##mv "error_acc.txt" ./clade1234_marginal
#mv "chisq_acc_allchromosome.txt" ./clade1234_marginal
