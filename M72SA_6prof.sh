#!/bin/bash
#$ -cwd -V
#$ -q long.q
#$ -m beas
#$ -M c.weerasuriya@lshtm.ac.uk
#$ -N M72SA_xt2
#$ -l mem_free=2G,h_vmem=2G
#$ -R y
#$ -j y
#$ -o /home/lsh1604836/clusterOE
#$ -t 1-6

R CMD BATCH "$HOME"/M72CE/M72SA/modelfit_active_GSA_vHratio_M72.R "$HOME"/M72CE/M72SA/Vxoutput_M72_refit/"${SGE_TASK_ID}"_console.out
