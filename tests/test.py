import time
import gravitas
import numpy as np
import pytest


def test_over_nmax():
    with pytest.raises(SystemError):
        gravitas.earth_acceleration(np.array([[0, 0, 0]]), 361)


def test_negative_nmax():
    with pytest.raises(SystemError):
        gravitas.earth_acceleration(np.array([[0, 0, 0]]), -1)


def test_value():
    r_ecef = np.tile(np.array([[7000, 0, 0.0]], dtype=np.float64), (1, 1))
    y = gravitas.earth_acceleration(r_ecef, 16)
    assert np.allclose(
        y, np.array([[-8.14575422e-03, -2.02473403e-08, 3.56963279e-08]])
    )


def test_rows():
    r_ecef = np.tile(np.array([[7000, 0, 0.0]], dtype=np.float64), (1000, 1))
    y = gravitas.earth_acceleration(r_ecef, 16)
    assert y.shape == (r_ecef.shape[0], 3)
    for i in range(1, r_ecef.shape[0]):
        assert np.allclose(y[i, :], y[0, :])


if __name__ == "__main__":
    r_ecef = np.tile(np.array([[7000, 0, 0.0]], dtype=np.float64), (1000, 1))

    t1 = time.time()
    y = gravitas.moon_acceleration(r_ecef, 120)
    t2 = time.time()
    print(t2 - t1)
    print(y)
