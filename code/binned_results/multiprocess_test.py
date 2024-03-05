import multiprocessing.pool

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

pool.close()
pool.join()

for i in range(10):
    results[i].get()
