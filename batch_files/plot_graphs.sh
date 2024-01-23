#!/bin/bash
# Generate random distributions
#python3 "./../code/binned_results/weighted_random_distribution.py" --raise_path=1 --nside=2 --target=3000000 --save_path="./../data/binned_random_catalogues/sdss_nside_2.fits"
#python3 "./../code/binned_results/weighted_random_distribution.py" --raise_path=1 --nside=8 --target=3000000 --save_path="./../data/binned_random_catalogues/sdss_nside_8.fits"
#python3 "./../code/binned_results/weighted_random_distribution.py" --raise_path=1 --nside=32 --target=3000000 --save_path="./../data/binned_random_catalogues/sdss_nside_32.fits"

# Write new graphs

python3 "./../code/binned_results/results_plot_new.py" --raise_path=1 --catalogue=sdss --weight_version=1 --save_path="./../graphs/binned_results/sdss_sample_1.pdf"
python3 "./../code/binned_results/results_plot_new.py" --raise_path=1 --catalogue=sdss --weight_version=2 --save_path="./../graphs/binned_results/sdss_sample_2.pdf"
python3 "./../code/binned_results/results_plot_new.py" --raise_path=1 --catalogue=sdss --weight_version=3 --save_path="./../graphs/binned_results/sdss_sample_3.pdf"

python3 "./../code/binned_results/results_plot_new.py" --raise_path=1 --catalogue="./../data/binned_random_catalogues/sdss_nside_8.fits" --weight_version=1 --save_path="./../graphs/binned_results/sdss_random_1.pdf"
python3 "./../code/binned_results/results_plot_new.py" --raise_path=1 --catalogue="./../data/binned_random_catalogues/sdss_nside_8.fits" --weight_version=2 --save_path="./../graphs/binned_results/sdss_random_2.pdf"
python3 "./../code/binned_results/results_plot_new.py" --raise_path=1 --catalogue="./../data/binned_random_catalogues/sdss_nside_8.fits" --weight_version=3 --save_path="./../graphs/binned_results/sdss_random_3.pdf"
