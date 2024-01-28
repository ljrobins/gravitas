import time
# import os
# os.system('python setup.py clean --all && python setup.py build_ext --inplace')
import gravitas
import numpy as np

# r_ecef = np.tile(np.array([[7000, 0, 0.0]], dtype=np.float64), (1000,1))
r_ecef = np.tile(np.array([[7000, 0, 0.0]], dtype=np.float64), (10_000,1))

t1  = time.time()
g = gravitas.acceleration(r_ecef, 16, "EGM96")
t2 = time.time()
print(t2-t1)

print(g[0,:])