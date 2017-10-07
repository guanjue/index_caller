#PBS -l walltime=20:00:00
#PBS -l nodes=2:ppn=8
#PBS -l mem=80gb

cd /gpfs/home/gzx103/scratch/index_caller/try/test_data
time Rscript index_calling.R
