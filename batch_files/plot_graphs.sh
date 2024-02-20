#!/bin/bash
# Generate random distributions
#python3 "./../code/binned_results/weighted_random_distribution.py" --raise_path=1 --nside=2 --target=3000000 --save_path="./../data/binned_random_catalogues/sdss_nside_2.fits" &
#python3 "./../code/binned_results/weighted_random_distribution.py" --raise_path=1 --nside=8 --target=3000000 --save_path="./../data/binned_random_catalogues/sdss_nside_8.fits" &
#python3 "./../code/binned_results/weighted_random_distribution.py" --raise_path=1 --nside=32 --target=3000000 --save_path="./../data/binned_random_catalogues/sdss_nside_32.fits" &
wait
# Write new graphs

#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss" --weight_function="excess_measurement" --save_path="./sdss_excess.pdf" --data_mask="sdss_planck" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="80k" --weight_function="excess_measurement" --save_path="./sdss_excess_random_80k.pdf" --data_mask="sdss_planck" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="400k" --weight_function="excess_measurement" --save_path="./sdss_excess_random_400k.pdf" --data_mask="sdss_planck" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="10m" --weight_function="excess_measurement" --save_path="./sdss_excess_random_10m.pdf" --data_mask="sdss_planck" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="planck_point_biased" --weight_function="excess_measurement" --save_path="./sdss_planck_point_biased.pdf" --data_mask="sdss_planck" &

python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss" --weight_function="excess_measurement" --save_path="./act_excess.pdf" --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="80k" --weight_function="excess_measurement" --save_path="./act_excess_random_80k.pdf" --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="400k" --weight_function="excess_measurement" --save_path="./act_excess_random_400k.pdf" --data_mask="sdss_act" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="10m" --weight_function="excess_measurement" --save_path="./act_excess_random_10m.pdf" --data_mask="sdss_act" &

python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./z_bin_1.png" --min_z=0.000 --max_z=0.319 --mask_set="point" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./z_bin_2.png" --min_z=0.319 --max_z=0.420 --mask_set="point" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./z_bin_3.png" --min_z=0.420 --max_z=0.489 --mask_set="point" &

python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./z_min_1.png" --min_z=0.319 --mask_set="point" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./z_min_2.png" --min_z=0.420 --mask_set="point" &
python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered" --weight_function="excess" --save_path="./z_min_3.png" --min_z=0.489 --mask_set="point" &

wait
