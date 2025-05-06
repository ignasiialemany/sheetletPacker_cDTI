import numpy as np
from shapely.affinity import translate
from shapely.geometry import Polygon

from .distances import midpoint, normalize, pdist, mindist


class PolyPacker:

    def __init__(self, N=0):
        """Initialise PolyPacker object.
        Use N to preallocate enough memory.
        """
        self.centers = np.empty((N, 2), dtype=float)  # [(x,y)_0, (x,y)_1, ...]
        self._polygons = np.empty(N, dtype=Polygon)
        # distance variables
        M = N*(N-1)//2  # number of elements in triangular matrix (without main diagonal)
        self._indices = np.empty((M, 2), dtype=int)
        self._minimum_distance = np.empty(M, dtype=float)

    def add_polygons(self, polygons):
        """Add polygons to the packer.
        Provide polygons in absolute coordinates.
        The polygons are a list of containing any combination of the following:
            - instances of shapely.geometry.Polygon objects,
            - objects with .vertices property OR method that returns xy data, or
            - xy data (a sequence of (x,y) pairs)
        """
        # process different data types
        polys = []
        for poly in polygons:
            if not isinstance(poly, Polygon):
                vert = getattr(poly, 'vertices', None)  # check if it has .vertices
                if vert:
                    if callable(vert):  # .vertices is a method
                        xy = vert()
                    else:  # .vertices is a property
                        xy = vert
                else:  # no attribute .vertices, so we assume it is plain xy data
                    xy = poly
                poly = Polygon(xy)  # construct the Polygon object
            polys.append(poly)
        # split representation of each polygon into local polygon (shape) and center (position)
        centers = [midpoint(p) for p in polys]
        _polygons = [translate(p, -x0, -y0) for p, (x0, y0) in zip(polys, centers)]
        self.centers = np.append(self.centers, centers, axis=0)
        self._polygons = np.append(self._polygons, _polygons)
        self.update_state()  # set internal state _minimum_distance, _indices

    # ============================== #
    #            Polygons            #
    # ============================== #

    def get_polygon(self, idx, FoR='local'):
        """Get polygon in the appropriate Frame-of-Reference.
        By default, returns in the local frame to avoid translation.
        """
        if FoR == 'global':
            return translate(self._polygons[idx], *self.centers[idx])
        elif FoR == 'local':
            return self._polygons[idx]
        else:
            raise Exception('invalid frame of reference (FoR)')

    def get_polygons(self, FoR='local'):
        """Return all polygons in the appropriate Frame-of-Reference.
        By default, returns in the local frame to avoid translation.
        """
        return [self.get_polygon(i, FoR=FoR) for i in range(self.N)]

    @property
    def N(self):
        return len(self._polygons)

    # We are keeping the actual data (_polygons) private, because of Frame-of-Reference
    polygons = property(get_polygons)  # get method using default arguments

    # ============================== #
    #            Distance            #
    # ============================== #

    def update_state(self):
        """Update the internal state.
        When updating _polygons, the other internal variables need to be updated.
        """
        # compute minimum clearance between objects
        #NOTE: uses _polygons (internal, translated into local frame)!
        self._minimum_distance, self._indices = mindist(self._polygons)

    def find_intersections(self):
        """Find all intersections between polygons.
        Returns the symmetric intersection matrix where
          [i,j] = True iff polygons i and j intersect.
        """

        # compute objects' pairwise distances
        distances = pdist(self.centers)

        # detect possible collisions based on spacing
        candidates = self._indices[distances <= self._minimum_distance]

        # check intersection and fill intersection matrix
        nPolys = len(self.polygons)
        intersection_matrix = np.zeros((nPolys, nPolys), dtype=bool)  #TODO: scipy.sparse?
        for i, j in candidates:
            # get polygons in their global position
            poly_i = self.get_polygon(i, FoR='global')
            poly_j = self.get_polygon(j, FoR='global')
            # check for intersection
            if poly_i.intersects(poly_j):  # did intersect
                intersection_matrix[i, j] = True  # \__ symmetric
                intersection_matrix[j, i] = True  # /

        return intersection_matrix

    # ============================== #
    #             Update             #
    # ============================== #

    def step(self, att=0, rep=0):
        """Update positions of polygons.
        Attract all non-overlapping polygons by att towards the origin,
        while repelling all overlapping polygons by rep.
        """

        # find intersections
        intersection_matrix = self.find_intersections()
        intersection_matrix = intersection_matrix[:, :, np.newaxis]  # reshape for broadcasting

        # repulsion
        xy1 = np.reshape(self.centers, (1, -1, 2))
        xy2 = np.reshape(xy1         , (-1, 1, 2))
        Unorm = normalize(xy1-xy2, axis=2)  # unit vectors to others
        Unorm = np.where(intersection_matrix, Unorm, 0)  # only intersections
        Usum = np.sum(Unorm, axis=0, keepdims=False)  # sum contributions
        unit_vector_rep = normalize(Usum, axis=1)

        # attraction
        unit_vector_att = normalize(-self.centers, axis=1)

        # compute change and apply
        d_xy = np.where(np.any(intersection_matrix, axis=0),
                        rep * unit_vector_rep,
                        att * unit_vector_att)
        self.centers += d_xy