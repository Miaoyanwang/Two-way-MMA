#!/bin/sh

mkdir "clade1234_marginal"
./pro_4comp_interaction.o -m kinship_impute_130.txt -p kinship_clade1234.txt -c Cov_miaoyan.txt -n ncount.txt -g host_SNP.txt -s pathogen_SNP.txt
#done


##mv "result_acc.txt" ./clade1234_marginal
##mv "trajectory_acc.txt" ./clade1234_marginal
##mv "error_acc.txt" ./clade1234_marginal
#mv "chisq_acc_allchromosome.txt" ./clade1234_marginal
