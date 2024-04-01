#!/bin/bash

#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_r_bin_1.png" --min_r=0.000 --max_r=10.0 --mask_set="point" --data_mask="sdss_act" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_r_bin_2.png" --min_r=10.0  --max_r=20.0 --mask_set="point" --data_mask="sdss_act" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_r_bin_3.png" --min_r=20.0  --max_r=50.0 --mask_set="point"  --data_mask="sdss_act" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_r_bin_4.png" --min_r=50.0 --max_r=100.0 --mask_set="point"  --data_mask="sdss_act" &

#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_r_min_1.png" --min_r=10.0   --mask_set="point" --data_mask="sdss_act" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_r_min_2.png" --min_r=20.0  --mask_set="point"  --data_mask="sdss_act" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_r_min_3.png" --min_r=50.0  --mask_set="point"  --data_mask="sdss_act" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./act_r_min_4.png" --min_r=100.0 --mask_set="point"  --data_mask="sdss_act" &

#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_r_bin_1.png" --min_r=0.000 --max_r=10.0 --mask_set="point" --data_mask="sdss_planck" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_r_bin_2.png" --min_r=10.0  --max_r=20.0 --mask_set="point" --data_mask="sdss_planck" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_r_bin_3.png" --min_r=20.0  --max_r=50.0 --mask_set="point" --data_mask="sdss_planck" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_r_bin_4.png" --min_r=50.0 --max_r=100.0 --mask_set="point" --data_mask="sdss_planck" &

#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_r_min_1.png" --min_r=10.0  --mask_set="point" --data_mask="sdss_planck" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_r_min_2.png" --min_r=20.0 --mask_set="point" --data_mask="sdss_planck" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_r_min_3.png" --min_r=50.0 --mask_set="point" --data_mask="sdss_planck" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_r_min_4.png" --min_r=100.0 --mask_set="point" --data_mask="sdss_planck" &

python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_s_3.png" --min_r=4.794 --max_r=15.85 --mask_set="both" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  --weight_function="excess" --save_path="./act_s_3.png" --min_r=4.794 --max_r=15.85 --mask_set="both" --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_s_4.png" --min_r=15.85 --max_r=34.15 --mask_set="both" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  --weight_function="excess" --save_path="./act_s_4.png" --min_r=15.85 --max_r=34.15 --mask_set="both" --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_s_5.png" --min_r=34.15 --max_r=59.29 --mask_set="both" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  --weight_function="excess" --save_path="./act_s_5.png" --min_r=34.15 --max_r=59.29 --mask_set="both" --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_s_6.png" --min_r=59.29 --max_r=90.64 --mask_set="both" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  --weight_function="excess" --save_path="./act_s_6.png" --min_r=59.29 --max_r=90.64 --mask_set="both" --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_s_7.png" --min_r=90.64 --max_r=127.5 --mask_set="both" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  --weight_function="excess" --save_path="./act_s_7.png" --min_r=90.64 --max_r=127.5 --mask_set="both" --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./planck_s_8.png" --min_r=127.5 --mask_set="both" --data_mask="sdss_planck" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  --weight_function="excess" --save_path="./act_s_8.png" --min_r=127.5 --mask_set="both" --data_mask="sdss_act" &

wait
