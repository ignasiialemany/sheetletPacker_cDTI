import numpy as np
from scipy.spatial import distance_matrix


def midpoint(polygon):
    """Compute the polygon mid-point.
    This is the mid-point between bounds, NOT the centroid! When enclosing the
    object in a circle, this guarantees the smallest radius, unlike when using
    the centroid.
    """
    x0, y0, x1, y1 = polygon.bounds
    return (x0+x1)/2, (y0+y1)/2


def normalize(vector, axis=1):
    """Turn vector into unit vector."""
    norm = np.sqrt(np.sum(vector**2, axis=axis, keepdims=True))
    normalised = np.divide(vector, norm, out=np.zeros_like(vector), where=norm>0)
    return normalised


def mat2arr(mat):
    """Convert distance matrix to array format."""
    arr = mat[np.mask_indices(mat.shape[0], np.tril, -1)]  # use tril
    idx = np.transpose(np.tril_indices(mat.shape[0], -1))
    return arr, idx


def mindist(polygons):
    """Compute the minimum spacing between object pairs."""
    radii = np.array([maxradius(p) for p in polygons])
    rad_1 = radii.reshape((-1, 1))
    rad_2 = rad_1.reshape((1, -1))
    min_dist = rad_1 + rad_2  # sum of two objects' radii as matrix
    minimum_distance, indices = mat2arr(min_dist)
    return minimum_distance, indices


def pdist(x0y0):
    """Calculate the Euclidean pairwise distance."""
    center_distances = distance_matrix(x0y0, x0y0)  # matrix
    center_distances, _ = mat2arr(center_distances)  # reduce to array
    return center_distances


def maxradius(polygon):
    """Find the maximum distance of any vertex to the origin.
    This assumes the polygon is centred around the origin.
    """
    xy = np.transpose(polygon.exterior.coords.xy)
    radius = np.max(distance_matrix(xy, [[0,0]]))
    return radius