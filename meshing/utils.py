from shapely.geometry import Polygon, MultiPolygon, GeometryCollection
from shapely.ops import unary_union
import numpy as np

def fix_self_intersections(polygon):
    """
    Fixes self-intersections in a polygon.
    :param polygon: Shapely Polygon object.
    :return: Valid Shapely Polygon object or None if cannot be fixed.
    """
    if not polygon.is_valid:
        fixed = make_valid(polygon)
        if isinstance(fixed, (MultiPolygon, GeometryCollection)):
            # Handle the case where fixed is a collection
            polygons = []
            for geom in fixed.geoms:
                if isinstance(geom, Polygon):
                    polygons.append(geom)
                elif isinstance(geom, MultiPolygon) or isinstance(geom, GeometryCollection):
                    polygons.extend([p for p in geom.geoms if isinstance(p, Polygon)])
            if not polygons:
                return None
            # Optionally, keep the largest polygon
            fixed = max(polygons, key=lambda p: p.area)
        elif isinstance(fixed, Polygon):
            return fixed
        else:
            return None
    return polygon if polygon.is_valid else None

def simplify_polygon(polygon, epsilon):
    """
    Simplifies a polygon using the Douglas-Peucker algorithm.
    :param polygon: Shapely Polygon object.
    :param epsilon: Tolerance for simplification.
    :return: Simplified Shapely Polygon object.
    """
    simplified = polygon.simplify(epsilon, preserve_topology=True)
    return simplified

def simplify_coordinates(vertices, epsilon):
    """
    Simplifies the coordinates of a polygon using the Ramer-Douglas-Peucker algorithm.
    :param vertices: Numpy array of vertices.
    :param epsilon: Tolerance for simplification.
    :return: Simplified numpy array of vertices.
    """
    from scipy.spatial import ConvexHull
    from shapely.geometry import LineString

    line = LineString(vertices)
    simplified_line = line.simplify(epsilon, preserve_topology=True)
    return np.array(simplified_line.coords)

def return_vertices(myocytes, centroids):
    polygons = []
    for myo in myocytes[0]:
        vertices = np.array(myo["vertices"] * 1000)
        vertices[:, 0] += centroids[0] * 1000
        vertices[:, 1] += centroids[1] * 1000
        
        # Remove points that are too close to each other
        vertices = simplify_coordinates(vertices, epsilon=0.1)
        
        # Ensure the polygon is closed by making the first and last points the same
        if not np.array_equal(vertices[0], vertices[-1]):
            vertices = np.vstack([vertices, vertices[0]])
        
        poly = Polygon(vertices)
        resolved_poly = fix_self_intersections(poly)
        if resolved_poly is not None:
            simplified_poly = simplify_polygon(resolved_poly, 0.1)
            #simplified_poly = simplify_polygon(simplified_poly.buffer(-0.01), 0.01)
            if simplified_poly.is_valid:
                polygons.append(simplified_poly)

    merged = unary_union(polygons)
    if isinstance(merged, (MultiPolygon, GeometryCollection)):
        final_polygons = [geom for geom in merged.geoms if isinstance(geom, Polygon)]
    else:
        final_polygons = [merged]

    return final_polygons