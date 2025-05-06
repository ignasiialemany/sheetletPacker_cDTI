import os.path

import numpy as np
from shapely.geometry import Polygon


def write_polygons(filepath, polygons, polydata=None, comment='Polygons', binary=False):
    """Create a (legacy) VTK file with polygons."""
    extension = ".vtk"  # needs the period
    if os.path.splitext(filepath)[1] != extension:
        filepath += extension
    with open(filepath, 'w') as vtk_file:
        vtk_file.write('# vtk DataFile Version 2.0\n')
        vtk_file.write('%s\n' % comment)
        vtk_file.write('%s\n' % ('BINARY' if binary else 'ASCII'))
        add_polygons_to_vtk(vtk_file, polygons, celldata=polydata, binary=binary)


def add_polygons_to_vtk(vtk_file, polygons, celldata=None, binary=False):

    if binary:
        raise NotImplementedError('cannot write binary files yet')

    # size up the problem
    nPolygons = len(polygons)
    nPolyPts = np.zeros(nPolygons, dtype=int)
    points = np.empty((0, 2))
    for i in range(nPolygons):
        xy = np.array(polygons[i].exterior.xy).transpose()
        nPolyPts[i] = xy.shape[0]
        points = np.append(points, xy, axis=0)
    nPoints = points.shape[0]
    nPolyPts_ranges = np.insert(np.cumsum(nPolyPts), 0, 0)

    # write to file
    vtk_file.write('DATASET POLYDATA\n')

    # write all point coordinates (2D -> 3D with z=0)
    vtk_file.write('POINTS %d float\n' % nPoints)
    for x, y in points:
        vtk_file.write('%g %g %g\n' % (x, y, 0))

    # write all polygon connectivities
    size = np.sum(nPolyPts)+nPolygons  # dont forget the integer for poly size
    vtk_file.write('POLYGONS %d %d\n' % (nPolygons, size))
    for i in range(nPolygons):
        vtk_file.write('%d' % nPolyPts[i])
        for idx in range(nPolyPts_ranges[i], nPolyPts_ranges[i+1]):
            vtk_file.write(' %d' % idx)
        vtk_file.write('\n')

    # (optionally) write ID
    if celldata is not None:
        vtk_file.write('CELL_DATA %d\n' % nPolygons)
        if type(celldata) is str:
            name = celldata.pop('name', 'celldata')
            dat = np.array(celldata['data'])  # should be 1D or 2D data
        else:
            name = 'celldata'
            dat = celldata
        dat = dat[:, np.newaxis] if dat.ndim == 1 else dat
        dim = dat.shape[1]
        type_ = get_type(dat)
        vtk_file.write('SCALARS %s %s %d\n' % (name, type_, dim))
        vtk_file.write('LOOKUP_TABLE default\n')
        for dat_i in dat:
            line = ('{} '*dim)[:-1].format(*dat_i)
            vtk_file.write(line+'\n')


def get_type(array):
    name = array.dtype.name
    if name.startswith("int"):
        return "int"
    elif name.startswith("float"):
        return "float"


def read_polygons(filepath):
    """Read polygons from a (legacy) VTK file."""

    with open(filepath, 'r') as vtk_file:
        header = vtk_file.readline().split()
        if header[0] != '#' or header[1] != 'vtk' or header[2] != 'DataFile':  # ignore Version
            raise Exception('invalid VTK file provided')
        _ = vtk_file.readline()  # comment
        ascii_or_binary = vtk_file.readline().rstrip()
        if ascii_or_binary.upper() != 'ASCII':
            raise NotImplementedError('can only read ascii files at the moment')
        dataset = vtk_file.readline().split()  # DATASET POLYDATA
        if dataset[1].upper() != 'POLYDATA':
            raise Exception('expected POLYDATA, got %s' % dataset[1])

        # read points
        points = vtk_file.readline().split()  # POINTS <N> <type>
        if points[0].upper() != 'POINTS':
            raise Exception('expected POINTS, got %s' % points[0])
        num_pts = int(points[1])
        pts = []
        for i in range(num_pts):
            xyz = vtk_file.readline().split()
            pts.append((float(xyz[0]), float(xyz[1])))

        # read polygons
        polygons = vtk_file.readline().split()  # POLYGONS <N> <M>
        if polygons[0].upper() != 'POLYGONS':
            raise Exception('expected POLYGONS, got %s' % polygons[0])
        num_poly = int(polygons[1])
        poly = []
        for i in range(num_poly):
            vertices = vtk_file.readline().split()
            _ = int(vertices[0])  # num_vert
            xy = [pts[int(v)] for v in vertices[1:]]
            poly.append(Polygon(xy))

        # read cell data
        celldata = vtk_file.readline().split()  # CELL_DATA <N>
        if len(celldata)==0 or celldata[0].upper() != 'CELL_DATA':
            # either at EOF or else we ignore the rest of file
            polydata = None

        else:

            num_cells = int(celldata[1])
            if num_cells != num_poly:
                raise Exception('expected same number of cells as polygons (%d), but got %d' % (num_poly, num_cells))
            scalars = vtk_file.readline().split()  # SCALARS <name> <type> <N>
            name = scalars[1]
            type_ = scalars[2]
            num_scalars = int(scalars[3])
            if scalars[0].upper() != 'SCALARS':
                raise Exception('expected SCALARS, got %s' % scalars[0])
            table = vtk_file.readline().split()  # LOOKUP_TABLE default
            if table[0].upper() != 'LOOKUP_TABLE':
                raise Exception('expected LOOKUP_TABLE, got %s' % table[0])
            polydata = np.empty((num_cells, num_scalars), dtype=type_)
            for i in range(num_cells):
                vals = vtk_file.readline().rstrip()
                polydata[i,:] = vals.split()
            polydata = dict(name=name, data=polydata)

        # done
        return poly, polydata
