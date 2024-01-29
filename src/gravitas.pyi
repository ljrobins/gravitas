from __future__ import annotations
import numpy
__all__ = ['earth_acceleration', 'moon_acceleration']
def earth_acceleration(arg0: numpy.ndarray, arg1: int) -> numpy.ndarray:
    """
            Acceleration due to the Earth's gravity field, as defined by the EGM96 model.
    
            :param r_ecef_km: Position in ECEF coordinates, in km.
            :type r_ecef_km: ndarray [nx3]
            :param nmax: Maximum degree of the spherical harmonic expansion, 0 <= nmax <= 360.
            :type nmax: int
            :return: Acceleration in ECI coordinates, in km/s^2.
            :rtype: float
    """
def moon_acceleration(arg0: numpy.ndarray, arg1: int) -> numpy.ndarray:
    """
            Acceleration due to the Moon's gravity field, as defined by the GRGM360 model.
    
            :param r_mcmf_km: Position in MCMF coordinates, in km.
            :type r_mcmf_km: ndarray [nx3]
            :param nmax: Maximum degree of the spherical harmonic expansion, 0 <= nmax <= 360.
            :type nmax: int
            :return: Acceleration in MCI coordinates, in km/s^2.
            :rtype: float
    """
