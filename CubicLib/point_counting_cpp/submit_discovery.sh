# loop over all fastq files in the directory, print the filename and submit the gzip jobs to SLURM
#
 
for FILE build/count_bigq_*; do
    echo ${FILE}
    sbatch -n 1 -t 1-00:00 --wrap="./${FILE}"
    sleep 1 # pause to be kind to the scheduler
done