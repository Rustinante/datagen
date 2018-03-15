#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 narrowpeak_filename output_dirname" >&2
  exit 1
fi

narrow_filename=$1
output_dirname=$2
echo "narrowpeak filename: ${narrow_filename}"
echo "output_dirname: ${output_dirname}"

python3 narrowpeak_to_fa.py $1 $2

cd ..
python3 -m gkm_datagen.species_letters_from_coord gkm_datagen/$2/$2.train.pos.coord gkm_datagen/$2
python3 -m gkm_datagen.species_letters_from_coord gkm_datagen/$2/$2.train.neg.coord gkm_datagen/$2
python3 -m gkm_datagen.species_letters_from_coord gkm_datagen/$2/$2.test.pos.coord  gkm_datagen/$2
python3 -m gkm_datagen.species_letters_from_coord gkm_datagen/$2/$2.test.neg.coord gkm_datagen/$2

