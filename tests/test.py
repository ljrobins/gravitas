import time
import gravitas
import numpy as np
import pytest


def test_over_nmax():
    """Test that nmax > 360 raises an error"""
    with pytest.raises(Exception):  # Rust panics, C++ raises SystemError
        gravitas.earth_acceleration(np.array([[0, 0, 0]], dtype=np.float64), 361)


def test_negative_nmax():
    """Test that negative nmax raises an error"""
    with pytest.raises(Exception):  # Rust raises OverflowError, C++ raises SystemError
        gravitas.earth_acceleration(np.array([[0, 0, 0]], dtype=np.float64), -1)


def test_value():
    """Test that earth_acceleration returns correct values for a known position"""
    r_ecef = np.tile(np.array([[7000, 0, 0.0]], dtype=np.float64), (1, 1))
    y = gravitas.earth_acceleration(r_ecef, 16)
    assert np.allclose(
        y, np.array([[-8.14575422e-03, -2.02473403e-08, 3.56963279e-08]])
    )


def test_rows():
    """Test that batch processing returns correct shape and consistent results"""
    r_ecef = np.tile(np.array([[7000, 0, 0.0]], dtype=np.float64), (1000, 1))
    y = gravitas.earth_acceleration(r_ecef, 16)
    assert y.shape == (r_ecef.shape[0], 3)
    for i in range(1, r_ecef.shape[0]):
        assert np.allclose(y[i, :], y[0, :])


def test_monopole():
    """Test monopole term (nmax=0) matches -mu/r^2"""
    mu = 398600.44  # km^3/s^2
    r = np.array([[6778.0, 0.0, 0.0]], dtype=np.float64)
    y = gravitas.earth_acceleration(r, 0)
    expected = -mu / (6778.0 * 6778.0)
    assert np.allclose(y[0, 0], expected, rtol=1e-10)


def test_symmetry():
    """Test that opposite positions give opposite accelerations (only for nmax=0)"""
    r_pos = np.array([[6778.0, 0.0, 0.0]], dtype=np.float64)
    r_neg = np.array([[-6778.0, 0.0, 0.0]], dtype=np.float64)
    
    # Only test symmetry for nmax=0 (monopole term)
    # For nmax > 0, the spherical harmonic terms break perfect symmetry
    y_pos = gravitas.earth_acceleration(r_pos, 0)
    y_neg = gravitas.earth_acceleration(r_neg, 0)
    assert np.allclose(y_pos[0, :] + y_neg[0, :], [0, 0, 0], rtol=1e-10)


def test_zero_position():
    """Test that zero position returns zero acceleration"""
    r = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
    y = gravitas.earth_acceleration(r, 10)
    assert np.allclose(y, [[0, 0, 0]])


if __name__ == "__main__":
    # Note: moon_acceleration is not implemented in the Rust version
    # r_ecef = np.tile(np.array([[7000, 0, 0.0]], dtype=np.float64), (1000, 1))
    # t1 = time.time()
    # y = gravitas.moon_acceleration(r_ecef, 120)
    # t2 = time.time()
    # print(t2 - t1)
    # print(y)
    pass
