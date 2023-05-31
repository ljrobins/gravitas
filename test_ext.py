# import pyspaceaware as ps

import os
os.system('python setup.py clean --all && python setup.py build_ext --inplace')
import gravitas._grav
import numpy as np

# r_ecef = np.tile(np.array([[7000, 0, 0.0]], dtype=np.float64), (1000,1))
r_ecef = np.tile(np.array([[7000, 0, 0.0], [0, 7000, 0]], dtype=np.float64), (1,1))

# ps.tic()
g = gravitas.acceleration(r_ecef, 360, "EGM96")
# ps.toc()

print(g)