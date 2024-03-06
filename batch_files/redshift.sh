#!/bin/bash

python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_z_bin_1.png" --min_z=0.000 --max_z=0.319 --mask_set="point" --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_z_bin_2.png" --min_z=0.319  --max_z=0.42 --mask_set="point" --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_z_bin_3.png" --min_z=0.42  --max_z=0.489 --mask_set="point"  --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_z_bin_4.png" --min_z=0.489 --max_z=0.5 --mask_set="point"  --data_mask="sdss_act" &

python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_z_min_1.png" --min_z=0.319   --mask_set="point" --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_z_min_2.png" --min_z=0.42  --mask_set="point"  --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_z_min_3.png" --min_z=0.489  --mask_set="point"  --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_z_min_4.png" --min_z=0.5 --mask_set="point"  --data_mask="sdss_act" &

python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_z_bin_1.png" --min_z=0.000 --max_z=0.319 --mask_set="point" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_z_bin_2.png" --min_z=0.319  --max_z=0.42 --mask_set="point" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_z_bin_3.png" --min_z=0.42  --max_z=0.489 --mask_set="point" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_z_bin_4.png" --min_z=0.489 --max_z=0.5 --mask_set="point" --data_mask="sdss_planck" &

python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_z_min_1.png" --min_z=0.319  --mask_set="point" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_z_min_2.png" --min_z=0.42 --mask_set="point" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_z_min_3.png" --min_z=0.489 --mask_set="point" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_z_min_4.png" --min_z=0.5 --mask_set="point" --data_mask="sdss_planck" &

wait
