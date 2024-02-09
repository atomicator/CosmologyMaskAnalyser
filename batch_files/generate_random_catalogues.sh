#!/bin/bash
# Generate random distributions
python3 "./../code/binned_results/weighted_random_distribution.py" --raise_path=1 --target=80000 --save_path="./../code/binned_results/80k.fits" &
python3 "./../code/binned_results/weighted_random_distribution.py" --raise_path=1  --target=400000 --save_path="./../code/binned_results/400k.fits" &
python3 "./../code/binned_results/weighted_random_distribution.py" --raise_path=1  --target=10000000 --save_path="./../code/binned_results/10m.fits" &
wait

nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="80k"  \
--weight_function="excess_measurement" --save_path="./80k.pdf" &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="400k"  \
--weight_function="excess_measurement" --save_path="./400k.pdf" &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="10m"  \
--weight_function="excess_measurement" --save_path="./10m.pdf" &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="80k"  \
--weight_function="excess_measurement" --save_path="./act_80k.pdf" --data_mask=sdss_act &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="400k"  \
--weight_function="excess_measurement" --save_path="./act_400k.pdf" --data_mask=sdss_act &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="10m"  \
--weight_function="excess_measurement" --save_path="./act_10m.pdf" --data_mask=sdss_act &
wait
