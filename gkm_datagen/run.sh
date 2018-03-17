#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 narrowpeak_filename output_dirname" >&2
  exit 1
fi

narrow_filename=$1
output_dirname=$2
echo "narrowpeak filename: ${narrow_filename}"
echo "output_dirname: ${output_dirname}"

. /u/local/Modules/default/init/modules.sh
module load python/3.6.1


python3 narrowpeak_to_fa.py $1 $2 || exit 1


cd ..
python3 -m gkm_datagen.species_letters_from_coord gkm_datagen/$2/$2.train.pos.coord gkm_datagen/$2 || exit 1
python3 -m gkm_datagen.species_letters_from_coord gkm_datagen/$2/$2.train.neg.coord gkm_datagen/$2 || exit 1
python3 -m gkm_datagen.species_letters_from_coord gkm_datagen/$2/$2.test.pos.coord  gkm_datagen/$2 || exit 1
python3 -m gkm_datagen.species_letters_from_coord gkm_datagen/$2/$2.test.neg.coord gkm_datagen/$2 || exit 1


echo "=> transporting files"
cd gkm_datagen
mkdir /u/home/a/aaronzho/project-ernst/lsgkm/tests/$2
python3 transport_files.py $2 /u/home/a/aaronzho/project-ernst/lsgkm/tests/$2 || exit 1


echo "=> submitting jobs"
cd /u/home/a/aaronzho/project-ernst/lsgkm/tests/${output_dirname}
echo "
#!/bin/bash
../task.sh $SGE_TASK_ID ${output_dirname}
" > jobarray.sh
qsub -cwd -V -N GKMJOB -pe shared 2 -l h_data=5G,h_rt=24:00:00,highp -M $HOME -m a -t 1-100:1 jobarray.sh