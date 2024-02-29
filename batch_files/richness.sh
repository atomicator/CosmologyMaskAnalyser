#!/bin/bash

python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_r_bin_1.png" --min_r=0.000 --max_r=10.0 --mask_set="point" --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_r_bin_2.png" --min_r=10.0  --max_r=20.0 --mask_set="point" --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_r_bin_3.png" --min_r=20.0  --max_r=50.0 --mask_set="point"  --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_r_bin_4.png" --min_r=50.0 --max_r=100.0 --mask_set="point"  --data_mask="sdss_act" &

python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_r_min_1.png" --min_r=10.0   --mask_set="point" --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_r_min_2.png" --min_r=20.0  --mask_set="point"  --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_r_min_3.png" --min_r=50.0  --mask_set="point"  --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_r_min_4.png" --min_r=100.0 --mask_set="point"  --data_mask="sdss_act" &

python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_r_bin_1.png" --min_r=0.000 --max_r=10.0 --mask_set="point" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_r_bin_2.png" --min_r=10.0  --max_r=20.0 --mask_set="point" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_r_bin_3.png" --min_r=20.0  --max_r=50.0 --mask_set="point" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_r_bin_4.png" --min_r=50.0 --max_r=100.0 --mask_set="point" --data_mask="sdss_planck" &

python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_r_min_1.png" --min_r=10.0  --mask_set="point" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_r_min_2.png" --min_r=20.0 --mask_set="point" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_r_min_3.png" --min_r=50.0 --mask_set="point" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_r_min_4.png" --min_r=100.0 --mask_set="point" --data_mask="sdss_planck" &

wait
