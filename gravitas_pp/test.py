import os
os.system('python setup.py clean --all && python setup.py build_ext --inplace')
import time
import gravitas_cpp
import numpy as np
gravitas_cpp.print_test()

import pytest

def test_over_nmax():
    with pytest.raises(SystemError):
        gravitas_cpp.acceleration(np.array([[0, 0, 0]]), 361, "EGM96")

def test_bad_name():
    with pytest.raises(SystemError):
        gravitas_cpp.acceleration(np.array([[0, 0, 0]]), 16, "EGM96s")

def test_negative_nmax():
    with pytest.raises(SystemError):
        gravitas_cpp.acceleration(np.array([[0, 0, 0]]), -1, "EGM96")

def test_value():
    r_ecef = np.tile(np.array([[7000, 0, 0.0]], dtype=np.float64), (1,1))
    y = gravitas_cpp.acceleration(r_ecef, 16, "EGM96")
    assert np.allclose(y, np.array([[-8.14575422e-03, -2.02473403e-08, 3.56963279e-08]]))

def test_rows():
    r_ecef = np.tile(np.array([[7000, 0, 0.0]], dtype=np.float64), (1000,1))
    y = gravitas_cpp.acceleration(r_ecef, 16, "EGM96")
    assert y.shape == (r_ecef.shape[0],3)
    for i in range(1,r_ecef.shape[0]):
        assert np.allclose(y[i,:], y[0,:])

if __name__ == "__main__":

    r_ecef = np.tile(np.array([[7000, 0, 0.0]], dtype=np.float64), (1000,1))

    t1 = time.time()
    y = gravitas_cpp.acceleration(r_ecef, 360, "EGM96")
    t2 = time.time()
    print(t2-t1)
    # print(y)