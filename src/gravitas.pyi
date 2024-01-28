from typing import Any
from numpy import ndarray

def earth_acceleration(r_ecef_km: ndarray, nmax: int) -> float:
    """Acceleration due to the Earth's gravity field, as defined by the EGM96 model.

    :param r_ecef_km: Position in ECEF coordinates, in km.
    :type r_ecef_km: ndarray [nx3]
    :param nmax: Maximum degree of the spherical harmonic expansion, 0 <= nmax <= 360.
    :type nmax: int
    :return: Acceleration in ECEF coordinates, in km/s^2.
    :rtype: float
    """
    pass

def moon_acceleration(r_mcmf_km: ndarray, nmax: int) -> float:
    """Acceleration due to the Earth's gravity field, as defined by the GRGM360 model produced by the GRAIL mission.

    :param r_mcmf_km: Position in MCMF coordinates, in km.
    :type r_mcmf_km: ndarray [nx3]
    :param nmax: Maximum degree of the spherical harmonic expansion, 0 <= nmax <= 360.
    :type nmax: int
    :return: Acceleration in MCMF coordinates, in km/s^2.
    :rtype: float
    """
    pass