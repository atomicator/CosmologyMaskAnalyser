import numpy as np
import multiprocessing.pool
import os
import time

os.system("taskset -p 0xff %d" % os.getpid())

print(multiprocessing.cpu_count())

def test(_data):
    os.system("taskset -p 0xff %d" % os.getpid())
    data = np.random.randn(2, 1000)
    for nn in range(1000):
        for ii in data[0, :]:
            for jj in data[1, :]:
                ii * jj
    return True


results = []
pool = multiprocessing.pool.ThreadPool(processes=10)
for i in range(10):
    results.append(pool.apply_async(test, (i,)))
    print(i)

time.sleep(100)

pool.close()
pool.join()

for i in range(10):
    results[i].get()
