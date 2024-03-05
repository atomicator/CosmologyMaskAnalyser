import multiprocessing.pool
import os

os.system("taskset -p 0xff %d" % os.getpid())

def test(x):
    while True:
        #print(x)
        x += 1
        #time.sleep(10)
    return True


results = []
pool = multiprocessing.pool.ThreadPool(processes=10)
for i in range(10):
    results.append(pool.apply_async(test, (i,)))

pool.close()
pool.join()

for i in range(10):
    results[i].get()
