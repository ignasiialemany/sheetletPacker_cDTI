import numpy as np


def overlap(polygons, region=None, rel_abs='abs'):
    """Calculate the total overlap of all polygons with the region (also a polygon)."""
    # check if region was provided
    if region is None:
        return 0
    # calculate overlap
    overlap_area = 0
    for polygon in polygons:
        intersection = region.intersection(polygon)
        overlap_area += intersection.area
    if rel_abs == 'abs':
        return overlap_area
    elif rel_abs == 'rel':
        return overlap_area/region.area
    else:
        raise Exception("rel_abs must be 'rel' or 'abs'")


def update(packer, region):

    # polygons
    polygons = packer.get_polygons(FoR='global')

    # collisions
    intersection_matrix = packer.find_intersections()
    intersections = np.sum(np.triu(intersection_matrix, 1), 1)
    collisions = np.sum(intersections[:-1])  # last one is self, so omit it

    # density
    density = overlap(polys, region, rel_abs='rel')

    # return
    return polygons, density, collisions