#read .mat file 
import scipy.io as sio
import numpy as np
import os
from shapely.geometry import Polygon, MultiPolygon,GeometryCollection
from shapely.ops import unary_union
from shapely.validation import make_valid
import matplotlib.pyplot as plt
from .mesh_generation import MeshGenerator
from .utils import fix_self_intersections, simplify_polygon, simplify_coordinates, return_vertices

#centroids = [
#     [0.022600500000000003, 0.2829195],
#     [0.24740850000000003, 0.0375915], 
#     [0.2508285, 0.0927675],
#     [0.1904085, 0.15341549999999998], 
#     [0.07093650000000001, 0.28497150000000004],
#     [0.1128885, 0.2409675],
#     [0.1445805, 0.3084555],
#     [0.20454450000000002, 0.3038955], 
#     [0.29414850000000003, 0.2918115], 
#     [0.37235250000000003, 0.1340355], 
#     [0.3728085, 0.28451550000000003], 
#     [0.44326050000000006, 0.2788155]]

def mesh_sheetlet(polys, folder_name, extrusion_layers=20):
    
    folder = f"meshing/{folder_name}/"
    #if folder does not exist, create it
    if not os.path.exists(folder):
        os.makedirs(folder)
    
    # Check for overlaps
    for i in range(len(polys)):
        for j in range(i+1, len(polys)):
            if polys[i].intersects(polys[j]):
                print(f"Overlap detected")    
                
    # Generate mesh if requested
    heights = np.random.uniform(114, 126, len(polys))
    lower_height = np.random.uniform(0.5, 1.5, len(polys))
    for i, poly in enumerate(polys):
        vertices = np.array(poly.exterior.coords.xy).T
        mesh_gen = MeshGenerator(vertices, 20, 126.6)
        h = np.round(heights[i], 2)
        l_h = np.round(lower_height[i], 2)
        mesh_gen.create_mesh(2, folder + f"myocyte_{i}.msh", h - l_h, z_offset=l_h, extrusion_layers=extrusion_layers)

    print("Mesh created")   
    return polys

