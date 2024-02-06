#!/bin/bash
# Generate random distributions
python3 "./../code/binned_results/weighted_random_distribution.py" --raise_path=1 --target=80000 --save_path="./../data/binned_random_catalogues/80k.fits" &
python3 "./../code/binned_results/weighted_random_distribution.py" --raise_path=1  --target=400000 --save_path="./../data/binned_random_catalogues/400k.fits" &
python3 "./../code/binned_results/weighted_random_distribution.py" --raise_path=1  --target=10000000 --save_path="./../data/binned_random_catalogues/10m.fits" &
wait

nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="80k"  \
--weight_function="excess_measurement" --save_path="./80k.pdf" &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="400k"  \
--weight_function="excess_measurement" --save_path="./400k.pdf" &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="10m"  \
--weight_function="excess_measurement" --save_path="./10m.pdf"
wait
