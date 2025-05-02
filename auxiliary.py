from shapely.ops import unary_union
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
from shapely.affinity import scale, translate
import numpy as np
from shapely.geometry import Polygon, MultiPolygon, LineString, GeometryCollection
from shapely.ops import unary_union, split

def generate_random_split_line(polygon):
    # Calculate the centroid of the polygon
    centroid = polygon.centroid
    
    # Generate a random angle in radians
    angle = np.random.uniform(0, 2 * np.pi)
    
    # Create a long line passing through the centroid at the random angle
    length = max(polygon.bounds[2] - polygon.bounds[0], polygon.bounds[3] - polygon.bounds[1])
    line = LineString([
        (centroid.x + length * np.cos(angle), centroid.y + length * np.sin(angle)),
        (centroid.x - length * np.cos(angle), centroid.y - length * np.sin(angle))
    ])
    
    return line

def split_polygon(polygon):
    # Generate a random split line
    split_line = generate_random_split_line(polygon)
    
    # Split the polygon using the line
    split_result = split(polygon, split_line)
    
    list_polygons = []
    if isinstance(split_result, GeometryCollection):
        for geom in split_result.geoms:
            if isinstance(geom, Polygon):
                if geom.area > 20 / (1000 * 1000):
                    geom = geom.simplify(1., preserve_topology=True)
                    list_polygons.append(geom)
            continue
            
    
    return list_polygons    
    
def solve_overlaps(cells, intersection_threshold, area_threshold):
    merged_cells = []
    
    while cells:
        cell = cells.pop(0)
        cell = cell.buffer(0)
        to_merge = [cell]
        
        for other_cell in cells[:]:
            if cell.intersects(other_cell):
                intersection_area = cell.intersection(other_cell).area
                min_area = min(cell.area, other_cell.area)
                if intersection_area > intersection_threshold:
                    to_merge.append(other_cell)
                    cells.remove(other_cell)
        
        if len(to_merge) > 1:
            merged_cell = unary_union(to_merge)
            if merged_cell.area > area_threshold:
                split_result = split_polygon(merged_cell)
                merged_cells.extend(split_result)
            else:
                merged_cells.append(merged_cell)
        else:
            if cell.area > area_threshold:
                split_result = split_polygon(cell)
                merged_cells.extend(split_result)
            else:
                merged_cells.append(cell)
    
    return merged_cells

# Example usage
# cells = list of your Polygon or MultiPolygon objects
# intersection_threshold = threshold value to decide if the intersection is large
# area_threshold = threshold value to decide if a cell should be split
# processed_cells = solve_overlaps(cells, intersection_threshold, area_threshold)

def get_cells_in_sheetlet(polys, sheetlet):
    inner_cells = []
    for poly in polys:
        centroid = poly.centroid
        random_scale = 1/1.35
        scaled_poly = scale(poly, xfact=random_scale, yfact=random_scale, origin=(0, 0))
        shifted_poly = translate(scaled_poly, xoff=centroid.x - scaled_poly.centroid.x, yoff=centroid.y - scaled_poly.centroid.y)
        poly = shifted_poly
        overlap = poly.intersection(sheetlet)
        if overlap.geom_type == 'Polygon':
            if overlap.area/poly.area > 0.5:
                inner_cells.append(overlap)
        elif overlap.geom_type == 'MultiPolygon':
            if overlap.area/poly.area > 0.5:
                inner_cells.append(overlap)
    
    # Now, handle overlaps between cells
    i = 0
    while i < len(inner_cells):
        cell = inner_cells[i]
        j = i + 1
        while j < len(inner_cells):
            other_cell = inner_cells[j]
            cell_overlap = cell.intersection(other_cell)
            if cell_overlap.area / min(cell.area, other_cell.area) > 0 and cell_overlap.area / min(cell.area, other_cell.area) > 0.7:
                # Merge the two cells if overlap is less than 20% of the smaller cell's area
                merged_cell = unary_union([cell, other_cell])
                inner_cells[i] = merged_cell  # Update the first cell with the merged one
                inner_cells.pop(j)  # Remove the second cell
            else:
                j += 1
        i += 1
    
    return inner_cells


def get_inner_cells(polys, sheetlet):
    inner_cells = []
    
    # First, handle intersection with the sheetlet
    for poly in polys:
        overlap = poly.intersection(sheetlet)
        if overlap.geom_type == 'Polygon':
            if overlap.area / poly.area > 0.5:
                inner_cells.append(overlap)
        elif overlap.geom_type == 'MultiPolygon':
            if overlap.area / poly.area > 0.5:
                inner_cells.extend(overlap.geoms)

    return inner_cells


def compute_num_cells(polys, sheetlet):
    area = sheetlet.area*1000*1000
    expected_cells = area/150
    cells = get_inner_cells(polys, sheetlet)
    return len(cells) >= expected_cells, expected_cells, len(cells)

def compute_ics_square(polys,square):
    cells  = get_inner_cells(polys, square)
    union_area = unary_union(cells).area
    return union_area/square.area 

def compute_overlap( polys):
    total_area = sum(poly.area for poly in polys)
    union_area = unary_union(polys).area
    overlap = total_area - union_area        
    return overlap/total_area

def compute_total_overlap(polygons):
    total_individual_area = sum(poly.area for poly in polygons)
    union_poly = unary_union(polygons)
    total_union_area = union_poly.area
    total_overlap_area = total_individual_area - total_union_area
    return total_overlap_area


def compute_ics(polys, sheetlet):
    inner_cells = get_inner_cells(polys, sheetlet)
    #print("THIS AMOUNT OF CELLS IN", len(inner_cells))
    total_cell_area = sum([cell.area for cell in inner_cells])
    sheetlet_area = sheetlet.area
    ics = total_cell_area / sheetlet_area
    return ics

def adaptive_attraction(current_attraction, ics_change, increment=0.0001, decrement=0.0003,threshold=0.0035):
    if abs(ics_change) > threshold:
        #print("Decreasing attraction")
        return max(0, current_attraction - decrement)
    else:
        #print("Increasing attraction")
        return current_attraction + increment
    
def adaptive_overlap_attraction(current_attraction, overlap_change, increment=0.0001, decrement=0.0003):
    if abs(overlap_change) > 0.03:
        #print("Decreasing attraction")
        return max(0, current_attraction - decrement)
    else:
        #print("Increasing attraction")
        return current_attraction + increment



def polygon_to_pixels(polygon, scale_x, scale_y, minx, miny, image_width, image_height):
    def scale_and_clip(coords):
        x, y = coords
        x_scaled = ((np.array(x) - minx) * scale_x).astype(int)
        y_scaled = ((np.array(y) - miny) * scale_y).astype(int)
        x_scaled = np.clip(x_scaled, 0, image_width - 1)
        y_scaled = np.clip(y_scaled, 0, image_height - 1)
        return list(zip(x_scaled, y_scaled))
    
    if isinstance(polygon, Polygon):
        # Process a single Polygon
        return scale_and_clip(polygon.exterior.xy)
    elif isinstance(polygon, MultiPolygon):
        #check if polygon is valid
        polygon = polygon.buffer(0)
        # Process each Polygon in the MultiPolygon
        all_coords = []
        for poly in polygon.geoms:
            all_coords.extend(scale_and_clip(poly.exterior.xy))
        return all_coords
    else:
        raise TypeError("Input must be a Polygon or MultiPolygon")





