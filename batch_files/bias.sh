#!/bin/bash

python3 ../code/bias_test/bias_test.py --save_path="random_act_400k.npy" --raise_dir=1  \
--threads=100 --iterations=100 --target=400000 --overdensity=0.0 --invert_bias=False --data_mask="sdss_act"
python3 ../code/bias_test/bias_test.py --save_path="random_act_100k.npy" --raise_dir=1  \
--threads=100 --iterations=100 --target=10000 --overdensity=0.0 --invert_bias=False --data_mask="sdss_act"
python3 ../code/bias_test/bias_test.py --save_path="random_planck_400k.npy" --raise_dir=1  \
--threads=100 --iterations=100 --target=400000 --overdensity=0.0 --invert_bias=False --data_mask="sdss_planck"
python3 ../code/bias_test/bias_test.py --save_path="random_planck_10k.npy" --raise_dir=1  \
--threads=100 --iterations=100 --target=10000 --overdensity=0.0 --invert_bias=False --data_mask="sdss_planck"

python3 ../code/bias_test/bias_test.py --save_path="normal_bias_act_400k.npy" --raise_dir=1  \
--threads=100 --iterations=100 --target=400000 --overdensity=0.1 --invert_bias=False --data_mask="sdss_act"
python3 ../code/bias_test/bias_test.py --save_path="normal_bias_act_100k.npy" --raise_dir=1  \
--threads=100 --iterations=100 --target=10000 --overdensity=0.1 --invert_bias=False --data_mask="sdss_act"
python3 ../code/bias_test/bias_test.py --save_path="normal_bias_planck_400k.npy" --raise_dir=1  \
--threads=100 --iterations=100 --target=400000 --overdensity=0.1 --invert_bias=False --data_mask="sdss_planck"
python3 ../code/bias_test/bias_test.py --save_path="normal_bias_planck_10k.npy" --raise_dir=1  \
--threads=100 --iterations=100 --target=10000 --overdensity=0.1 --invert_bias=False --data_mask="sdss_planck"

python3 ../code/bias_test/bias_test.py --save_path="invert_bias_act_400k.npy" --raise_dir=1  \
--threads=100 --iterations=100 --target=400000 --overdensity=0.11111 --invert_bias=True --data_mask="sdss_act"
python3 ../code/bias_test/bias_test.py --save_path="invert_bias_act_100k.npy" --raise_dir=1  \
--threads=100 --iterations=100 --target=10000 --overdensity=0.11111 --invert_bias=True --data_mask="sdss_act"
python3 ../code/bias_test/bias_test.py --save_path="invert_bias_planck_400k.npy" --raise_dir=1  \
--threads=100 --iterations=100 --target=400000 --overdensity=0.11111 --invert_bias=True --data_mask="sdss_planck"
python3 ../code/bias_test/bias_test.py --save_path="invert_bias_planck_10k.npy" --raise_dir=1  \
--threads=100 --iterations=100 --target=10000 --overdensity=0.11111 --invert_bias=True --data_mask="sdss_planck"
