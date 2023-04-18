
To run a job on discovery:
1. Compile it
2. `sbatch <your bash script here>`
3. Read the output/errors from the error files.

NOTE: `squeue` monitors your job status.

TODO:
We should discuss with the discovery people the following topics:
1. How to allow our jobs to run for a month without timing out.
2. Optimizing cluster usage (Node and jobs per node tasks)
3. Should we compile before/after batching?

Links to documentation:
https://services.dartmouth.edu/TDClient/1806/Portal/KB/?CategoryID=21663
https://services.dartmouth.edu/TDClient/1806/Portal/KB/ArticleDet?ID=132625

