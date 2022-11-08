import numpy as np
from scipy.spatial import distance


def calc_angle(coordinates, a, b, c):
    """
    Calculate the angle between three atoms in space
    :param coordinates: 2D numpy array of coordinates from cclib.io.ccread ORCA coordinates
    :param a: index atom 1
    :param b: index of the centre point of the angle (atom 2)
    :param c: index of atom 3
    :return: Angle in degrees
    """
    ba = coordinates[0][b] - coordinates[0][a]
    bc = coordinates[0][b] - coordinates[0][c]
    angle = np.arccos(np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc)))
    return np.degrees(angle)


def calc_distance(coordinates, a, b):
    """
    Calculate the distance between two points
    :param coordinates: 2D numpy array of coordinates from cclib.io.ccread ORCA coordinates
    :param a: index of atom 1
    :param b: index of atom 2
    :return:
    """
    return distance.euclidean(coordinates[0][a], coordinates[0][b])