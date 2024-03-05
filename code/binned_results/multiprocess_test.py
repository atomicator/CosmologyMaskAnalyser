from toolkit import toolkit, weights
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import argparse
import multiprocessing.pool
import time

def test(x):
    while True:
        print(x)
        x += 1
        #time.sleep(10)
    return True


results = []
pool = multiprocessing.pool.ThreadPool(processes=5)
for i in range(10):
    results.append(pool.apply_async(test, (i,)))

for i in range(10):
    results[i].get()
