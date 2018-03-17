#!/bin/bash

species_index=$(($1 - 1))
dirname=$2

echo "=> species_index: ${species_index}"
echo "=> dirname ${dirname}"

cd ${dirname}

../../bin/gkmtrain -T25 -t0 -l10 -k6 -d3 -c1 ${dirname}.train.pos.coord.mult_species.at1.00/${species_index}_train_pos.fa ${dirname}.train.neg.coord.mult_species.at1.00/${species_index}_train_neg.fa species.${species_index} || exit 1

echo "=> Finished training species ${species_index}"

../../bin/gkmpredict -T25 ${dirname}.test.pos.coord.mult_species/${species_index}_test_pos.fa species_model/species.${species_index}.model.txt ${species_index}.pos_predict || exit 1

../../bin/gkmpredict -T25 ${dirname}.test.neg.coord.mult_species/${species_index}_test_neg.fa species_model/species.${species_index}.model.txt ${species_index}.neg_predict || exit 1

echo "=> Finished predicting species ${species_index}"

. /u/local/Modules/default/init/modules.sh
module load python/3.6.1

python3 ../roc.py ${species_index}.pos_predict ${species_index}.neg_predict ${species_index}.roc
echo "=> Finished calculating roc for ${species_index}"

