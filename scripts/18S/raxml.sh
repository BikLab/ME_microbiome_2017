#!/bin/bash

#SBATCH --job-name="raxml"
#SBACTH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=taruna.aggarwal@ucr.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e raxml-18s.err-%N
#SBATCH -o raxml-18s.out-%N


module load RAxML

IN=/rhome/taruna/shared/taruna/memb/18S/analysis-results/uclust_openref99_MEMB1_18s_mergedONLY_26July2017/pynast_aligned_seqs
OUT=/rhome/taruna/shared/taruna/memb/18S/analysis-results/uclust_openref99_MEMB1_18s_mergedONLY_26July2017/raxml

raxmlHPC-PTHREADS-SSE3 -T 8 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 1000 -x 02938 -n 18S-nemasONLY -s $IN/rep_set_aligned_pfiltered_nemasONLY_w_info.fasta -w $OUT
