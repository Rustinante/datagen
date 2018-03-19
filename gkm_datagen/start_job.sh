#!/bin/bash

. /u/local/Modules/default/init/modules.sh
module load python/3.6.1

python3 submit_massive_jobs.py dnase narrowfiles/dnasenarrowpeak/ $((${SGE_TASK_ID} - 1))

