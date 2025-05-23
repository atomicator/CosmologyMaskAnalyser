#!/bin/bash

file=bayesian_mark_two
mask="$1"
flat_prior=${2:-"False"}
target=80000

echo "${flat_prior}"

python3 ../code/binned_results/$file.py --path="${mask}_${target}_0.0.npy" --raise_dir=1 --processes=100 --flat_prior="${flat_prior}"  \
--threads=100 --realisations=100 --target="$target" --overdensity=0.0 --invert_bias=False --data_mask="${mask}"
python3 ../code/binned_results/$file.py --path="${mask}_${target}_0.1.npy" --raise_dir=1 --processes=100 --flat_prior="${flat_prior}"  \
--threads=100 --realisations=100 --target="$target" --overdensity=0.1 --invert_bias=False --data_mask="${mask}"
python3 ../code/binned_results/$file.py --path="${mask}_${target}_0.2.npy" --raise_dir=1 --processes=100 --flat_prior="${flat_prior}"  \
--threads=100 --realisations=100 --target="$target" --overdensity=0.2 --invert_bias=False --data_mask="${mask}"
python3 ../code/binned_results/$file.py --path="${mask}_${target}_-0.1.npy" --raise_dir=1 --processes=100 --flat_prior="${flat_prior}"  \
--threads=100 --realisations=100 --target="$target" --overdensity=0.1111111 --invert_bias=True --data_mask="${mask}"
python3 ../code/binned_results/$file.py --path="${mask}_${target}_-0.2.npy" --raise_dir=1 --processes=100 --flat_prior="${flat_prior}"  \
--threads=100 --realisations=100 --target="$target" --overdensity=0.25 --invert_bias=True --data_mask="${mask}"