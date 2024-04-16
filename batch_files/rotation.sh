#!/bin/bash

mask_name=$1

#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-1.png" --lon_shift=-1.0 --data_mask="$sdss_{mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-2.png" --lon_shift=-2.0 --data_mask="$sdss_{mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-3.png" --lon_shift=-3.0 --data_mask="$sdss_{mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-4.png" --lon_shift=-4.0 --data_mask="$sdss_{mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-5.png" --lon_shift=-5.0 --data_mask="$sdss_{mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &

python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-10.png" --lon_shift=-10.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-15.png" --lon_shift=-15.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-20.png" --lon_shift=-20.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
##python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-25.png" --lon_shift=-25.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-30.png" --lon_shift=-30.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#
wait
#
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-40.png" --lon_shift=-40.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-50.png" --lon_shift=-50.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-60.png" --lon_shift=-60.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-70.png" --lon_shift=-70.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-80.png" --lon_shift=-80.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#
#wait
#
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-90.png" --lon_shift=-90.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-100.png" --lon_shift=-100.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-110.png" --lon_shift=-110.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-120.png" --lon_shift=-120.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-130.png" --lon_shift=-130.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#
#wait
#
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-140.png" --lon_shift=-140.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-150.png" --lon_shift=-150.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-160.png" --lon_shift=-160.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-170.png" --lon_shift=-170.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
##python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_-180.png" --lon_shift=-180.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_1.png" --lon_shift=1.0 --data_mask="${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_2.png" --lon_shift=2.0 --data_mask="${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_3.png" --lon_shift=3.0 --data_mask="${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_4.png" --lon_shift=4.0 --data_mask="${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_5.png" --lon_shift=5.0 --data_mask="${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &

#wait
#
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_10.png" --lon_shift=10.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
##python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_15.png" --lon_shift=15.0 --data_mask="${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_20.png" --lon_shift=20.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
##python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_25.png" --lon_shift=25.0 --data_mask="${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_30.png" --lon_shift=30.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point" &
#
#wait
#
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_40.png" --lon_shift=40.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point"&
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_50.png" --lon_shift=50.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point"&
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_60.png" --lon_shift=60.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point"&
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_70.png" --lon_shift=70.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point"&
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_80.png" --lon_shift=80.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point"&
#
#wait
#
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_90.png" --lon_shift=90.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point"&
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_100.png" --lon_shift=100.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point"&
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_110.png" --lon_shift=110.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point"&
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_120.png" --lon_shift=120.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point"&
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_130.png" --lon_shift=130.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point"&
#
#wait
#
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_140.png" --lon_shift=140.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point"&
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_150.png" --lon_shift=150.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point"&
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_160.png" --lon_shift=160.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point"&
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_170.png" --lon_shift=170.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point"&
#python3 ../code/binned_results/debug_script.py --save_path="${mask_name}_rot_180.png" --lon_shift=180.0 --data_mask="sdss_${mask_name}" --raise_path=1 --catalogue="sdss_filtered" --min_r=20 --mask_set="point"&
#
#wait
#