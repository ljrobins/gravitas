import os
os.system('python setup.py clean --all && python setup.py build_ext --inplace')

import gravitas as gr
import numpy as np
import matplotlib.pyplot as plt
import pyspaceaware as ps

# degrees = np.array([4, 100, 360])
degrees = np.arange(1, 360, 1)
rv = np.tile(np.array([[7000, 1, 1]]), (1,1))
acc = np.zeros((degrees.size, 3))
dt = np.zeros_like(degrees, dtype=np.float64)
for i,deg in enumerate(degrees):
    print(deg)
    ps.tic()
    for j in range(3):
        acci = gr.acceleration(rv, deg, use_model='EGM96')
        print(acci)
        magdiff = np.linalg.norm(acci - acc[[i],:])
        acc[i,:] = acci
    dt[i] = ps.toc()

# plt.plot(degrees, dt)
# plt.show()

print(acc)

plt.plot(degrees, np.abs(acc))
plt.yscale('log')
plt.show()