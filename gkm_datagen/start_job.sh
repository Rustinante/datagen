#!/bin/bash

. /u/local/Modules/default/init/modules.sh
module load python/3.6.1

# python3 submit_massive_jobs.py dnase.np narrowfiles/dnase_narrowpeak/ $((${SGE_TASK_ID} - 1))

python3 submit_massive_jobs.py histone.np narrowfiles/histone_narrowpeak/ $((${SGE_TASK_ID} - 1))sh
