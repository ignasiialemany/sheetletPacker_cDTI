from shapely.geometry import Polygon,Point,MultiPolygon
from shapely.affinity import translate, rotate, scale
from shapely.ops import unary_union
from polypacker import PolyPacker
from auxiliary import compute_ics, compute_overlap, adaptive_attraction, compute_total_overlap, get_inner_cells, polygon_to_pixels, compute_ics_square,get_cells_in_sheetlet,compute_num_cells
from auxiliary import solve_overlaps
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw
from scipy.io import savemat
import multiprocessing
import os

# Configuration
CONFIG = {
    'output_dir': 'output',
    'images_dir': 'images',
    'matlab_inputs_dir': 'matlab_inputs',
    'cell_data_file': 'polys.npy',
    'units_per_pixel': 0.06/1000,
    'plot_dpi': 500,
    'final_plot_dpi': 800,
    'plot_figsize': (7, 7),
    'font_size': 13,
    'colors': {
        'cells': '#8a0101',
        'fillers': 'black',
        'sheetlet': '#00d149'
    }
}

def get_output_filename(prefix, index, seed, file_type, stage=None):
    """Generate consistent output filenames.
    
    Args:
        prefix: Base prefix for the filename
        index: Index number for the file
        seed: Random seed used for generation
        file_type: Type of file ('npy', 'png', 'mat')
        stage: Optional stage identifier (e.g., 'init', 'final', 'repulsion')
    
    Returns:
        str: Generated filename
    """
    base = f"{prefix}_sheetlet_index_{index}_seed_{seed}"
    if stage:
        base = f"{base}_{stage}"
    
    if file_type == 'npy':
        return f"{base}.npy"
    elif file_type == 'png':
        return f"cellgrid_{base}.png"
    elif file_type == 'mat':
        return f"{base}.mat"
    else:
        raise ValueError(f"Unsupported file type: {file_type}")

# Create necessary directories
for dir_path in [CONFIG['output_dir'], CONFIG['images_dir'], CONFIG['matlab_inputs_dir']]:
    os.makedirs(dir_path, exist_ok=True)

class SheetletPacker:
    
    def __init__(self, sheetlet, sheetletIndex, celltype="real"):
        """Initialize the SheetletPacker with a sheetlet geometry.
        
        Args:
            sheetlet: Either a Shapely Polygon or path to a .png/.npy file
            sheetletIndex: Index of the sheetlet
            celltype: Type of cells to use ("real" or "ellipse")
        """
        self.celltype = celltype
        self.sheetlets = []
        self.index = sheetletIndex
        if isinstance(sheetlet, Polygon):
            if sheetlet.geom_type != 'Polygon':
                raise ValueError("Sheetlets must be Shapely polygons")
            self.centroid = sheetlet.centroid
            self.sheetlet = sheetlet
            self.sheetlet = translate(self.sheetlet, -self.centroid.x, -self.centroid.y)
        elif isinstance(sheetlet, str):
            self.loadSheetlet(sheetlet)
        else:
            raise ValueError("sheetlets must be a list of Shapely polygons, np.ndarray, or list of .png paths")
        
    def loadSheetlet(self, path):
        """Load a sheetlet from a file.
        
        Args:
            path: Path to the sheetlet file (.png or .npy)
        """
        if path.split('.')[-1] == 'png':
            self.sheetlets.append(self.loadPNG(path))
        elif path.split('.')[-1] == 'npy':
            cells, sheetlet, centroid = np.load(path, allow_pickle=True)
            self.centroid = centroid
            self.sheetlet = sheetlet
        else:
            raise ValueError("Sheetlet file must be .png or .npy")

    def loadPNG(self, path):
        """Load a sheetlet from a PNG file.
        
        Args:
            path: Path to the PNG file
        """
        pass   
    
    def init_cell(self, x, y, tol=0.1):
        """Initialize a cell at the specified coordinates.
        
        Args:
            x: x-coordinate
            y: y-coordinate
            tol: Tolerance for simplification
        
        Returns:
            Shapely Polygon: The initialized cell
        """
        cells = np.load("polys.npy", allow_pickle=True)
        isValid = False
        while not isValid:
            random_int = np.random.randint(0, len(cells))
            cell = cells[random_int]
            if self.celltype == "real":
                if cell.centroid.x != 0 or cell.centroid.y != 0:
                    cell = translate(cell, -cell.centroid.x, -cell.centroid.y)
                cell = cell.simplify(tol, preserve_topology=True)
            elif self.celltype == "ellipse":
                minx, miny, maxx, maxy = cell.bounds
                a = (maxx - minx) / 2
                b = (maxy - miny) / 2
                cell = Point(0, 0).buffer(1)
                cell = scale(cell, xfact=a, yfact=b)
            if cell.is_valid and cell.geom_type == 'Polygon':
                isValid = True
        random_scale = 1.38
        cell = scale(cell, xfact=random_scale, yfact=random_scale, origin=(0, 0))
        cell = scale(cell, xfact=1/1000, yfact=1/1000, origin=(0, 0))
        cell = rotate(cell, np.random.uniform(0, 180))
        cell = translate(cell, x, y)
        return cell
    
    def init_packer(self, sheetlet, isPlot=False):
        """Initialize the packer with cells and fillers.
        
        Args:
            sheetlet: The sheetlet geometry
            isPlot: Whether to generate plots
        
        Returns:
            PolyPacker: The initialized packer
        """
        centroid_sheetlet = sheetlet.centroid
        minx, miny, maxx, maxy = sheetlet.bounds
        cells, fillers = [], []
        big_length = 1.8*max(maxx - minx, maxy - miny) / 2
        small_length = 1.1*max(maxx - minx, maxy - miny) / 2
        
        self.large_square = Polygon([(-big_length, -big_length), (-big_length, big_length), (big_length, big_length), (big_length, -big_length)])
        self.small_square = Polygon([(-small_length, -small_length), (-small_length, small_length), (small_length, small_length), (small_length, -small_length)])
        
        if isPlot:
            plt.plot(*self.large_square.exterior.xy)
            plt.plot(*self.small_square.exterior.xy)
            plt.plot(*sheetlet.exterior.xy)
            plt.show()
        
        grid_size = int((2*big_length)/(15/1000))
        points = []
        for i in range(grid_size):
            for j in range(grid_size):
                x = (i + 0.5) * 15/1000 - big_length
                y = (j + 0.5) * 15/1000 - big_length
                points.append([x,y])
        
        counter = 0
        area_cells = 0
        
        while points:
            point = points.pop()
            cell = self.init_cell(point[0], point[1])            
            area_cells += cell.area
            counter += 1
            cells.append(cell)
            self.large_square = self.large_square.difference(cell)
          
        ecv = 0.00045
        total_filler_area = area_cells * (ecv/(1-ecv))
        area_cells = 0
        area_fillers = 0
        
        while area_fillers < total_filler_area:
            x, y = np.random.uniform(-small_length, small_length), np.random.uniform(-small_length,small_length)
            radii = np.random.uniform(1,3.5)
            filler = Point(x, y).buffer(radii / 1000)
            if self.large_square.contains(filler.centroid):
                area_fillers += filler.area
                fillers.append(filler)
        
        self.n_cells = len(cells)
        self.n_fillers = len(fillers)
        
        packer = PolyPacker()
        packer.add_polygons(cells)
        packer.add_polygons(fillers)
        
        compute_ics(cells,sheetlet)
        
        import matplotlib
        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{siunitx} \sisetup{detect-all}'
        matplotlib.rcParams['font.size'] = 13
        matplotlib.rcParams['axes.titlesize'] = 13
        matplotlib.rcParams['axes.labelsize'] = 13
        matplotlib.rcParams['xtick.labelsize'] = 13
        matplotlib.rcParams['ytick.labelsize'] = 13
        matplotlib.rcParams['legend.fontsize'] = 12
        matplotlib.rcParams['figure.titlesize'] = 13

        if isPlot:
            plt.figure(figsize=CONFIG['plot_figsize'])
            first_cell = True
            first_filler = True
            for poly in packer.get_polygons(FoR='global')[:self.n_cells]:
                if first_cell:
                    plt.fill(poly.exterior.xy[0], poly.exterior.xy[1], 
                            color=CONFIG['colors']['cells'], label="Cells")
                    first_cell = False
                else:
                    plt.fill(poly.exterior.xy[0], poly.exterior.xy[1], 
                            color=CONFIG['colors']['cells'])
            for poly in packer.get_polygons(FoR='global')[self.n_cells:]:
                if first_filler:
                    plt.fill(*poly.exterior.xy, 
                            color=CONFIG['colors']['fillers'], label="Fillers")
                    first_filler = False
                else:
                    plt.fill(*poly.exterior.xy, 
                            color=CONFIG['colors']['fillers'])
            
            plt.plot(*self.sheetlet.exterior.xy, 
                    color=CONFIG['colors']['sheetlet'], label="Sheetlet")
            plt.xlabel('X axis ($\si{\milli\meter}$)')
            plt.ylabel('Y axis ($\si{\milli\meter}$)')
            plt.axis('equal')
            plt.tight_layout()
            plt.legend(loc="upper right")
            plt.savefig(os.path.join(CONFIG['images_dir'], "cellgrid_init.png"), 
                       dpi=CONFIG['plot_dpi'])
            plt.close()
        return packer

    def pack(self, filename, seed, isPlot=False):
        """Pack cells into the sheetlet.
        
        Args:
            filename: Base filename for outputs
            seed: Random seed for reproducibility
            isPlot: Whether to generate plots
        """
        packer = self.init_packer(self.sheetlet, isPlot=isPlot)
        
        iterations = []
        ics_values = []
        steps = 10
        rep = 0.001 
        
        for i in range(steps):
            packer.step(0.0000001, rep)
            if isPlot and i % 5 == 0:
                print("Initial Phase: Step {0}".format(i))
            ics = compute_ics_square(packer.get_polygons(FoR="global")[:self.n_cells], self.sheetlet)
            ics_values.append(ics)
            iterations.append(i)

        if isPlot:
            plt.figure(figsize=CONFIG['plot_figsize'])
            for poly in packer.get_polygons(FoR='global')[:self.n_cells]:
                plt.fill(poly.exterior.xy[0], poly.exterior.xy[1], 
                        color=CONFIG['colors']['cells'])
            for poly in packer.get_polygons(FoR='global')[self.n_cells:]:
                plt.fill(*poly.exterior.xy, 
                        color=CONFIG['colors']['fillers'])
        
            plt.plot(*self.sheetlet.exterior.xy, 
                    color=CONFIG['colors']['sheetlet'])
            plt.xlabel('X axis ($\si{\milli\meter}$)')
            plt.ylabel('Y axis ($\si{\milli\meter}$)')
            plt.axis('equal')
            plt.tight_layout()
            plt.savefig(os.path.join(CONFIG['images_dir'], 
                                   get_output_filename(filename, self.index, seed, 'png', 'repulsion')),
                       dpi=CONFIG['plot_dpi'])
            plt.close()
                
        for i in range (1000000):
            rep = 1/1000
            attraction = rep*0.83
            packer.step(attraction,rep)
            ics = compute_ics_square(packer.get_polygons(FoR="global")[:self.n_cells], self.sheetlet)
            iterations.append(iterations[-1]+1)
            ics_values.append(ics)
            print("ICS {0}, step {1}".format(ics,i))
            if ics > 0.75:
                break
            if isPlot and i % 100 == 0:
                plt.figure(figsize=CONFIG['plot_figsize'])
                for poly in packer.get_polygons(FoR='global')[:self.n_cells]:
                    plt.fill(poly.exterior.xy[0], poly.exterior.xy[1], 
                            color=CONFIG['colors']['cells'])
                for poly in packer.get_polygons(FoR='global')[self.n_cells:]:
                    plt.fill(*poly.exterior.xy, 
                            color=CONFIG['colors']['fillers'])
                plt.plot(*self.sheetlet.exterior.xy, 
                        color=CONFIG['colors']['sheetlet'])
                plt.xlabel('X axis ($\si{\milli\meter}$)')
                plt.ylabel('Y axis ($\si{\milli\meter}$)')
                plt.axis('equal')
                plt.tight_layout()
                plt.savefig(os.path.join(CONFIG['images_dir'], 
                                       get_output_filename(filename, self.index, seed, 'png', f'step_{i}')),
                           dpi=CONFIG['plot_dpi'])
                plt.close()
                    
        all_polys = packer.get_polygons(FoR="global")
        inner_cells = get_cells_in_sheetlet(all_polys[:self.n_cells], self.sheetlet)
        
        if isPlot:
            plt.figure()
            plt.plot(iterations, ics_values, "-o", markersize=3.5, linewidth=0.5)
            plt.xlabel("Iterations")
            plt.ylabel("Intracellular Area Ratio")
            plt.savefig(os.path.join(CONFIG['images_dir'], 
                                   get_output_filename(filename, self.index, seed, 'png', 'ics')),
                       dpi=CONFIG['plot_dpi'])
            plt.close()
            
            plt.figure(figsize=CONFIG['plot_figsize'])
            for poly in inner_cells:
                if isinstance(poly,Polygon):
                    plt.plot(poly.exterior.xy[0], poly.exterior.xy[1],
                            color=CONFIG['colors']['cells'], linewidth=0.5)
                elif isinstance(poly,MultiPolygon):
                    for p in poly.geoms:
                        plt.plot(p.exterior.xy[0], p.exterior.xy[1],
                                color=CONFIG['colors']['cells'], linewidth=0.5)
            plt.plot(*self.sheetlet.exterior.xy)
            plt.axis("equal")
            plt.savefig(os.path.join(CONFIG['images_dir'], 
                                   get_output_filename(filename, self.index, seed, 'png', 'final')),
                       dpi=CONFIG['final_plot_dpi'])
            plt.close()
                
        np.save(os.path.join(CONFIG['output_dir'], 
                           get_output_filename(filename, self.index, seed, 'npy')),
                np.array([inner_cells, self.sheetlet, self.centroid], dtype=object))

        self.generate_matlab_images(inner_cells, self.sheetlet, filename, seed)
        
    def generate_matlab_images(self, inner_cells, sheetlet, filename, seed, units_per_pixel=None):
        """Generate images for MATLAB processing.
        
        Args:
            inner_cells: List of cells inside the sheetlet
            sheetlet: The sheetlet geometry
            filename: Base filename for outputs
            seed: Random seed
            units_per_pixel: Resolution for image generation
        """
        if units_per_pixel is None:
            units_per_pixel = CONFIG['units_per_pixel']

        minx, miny, maxx, maxy = sheetlet.bounds
        width = maxx - minx
        height = maxy - miny
        resolution_x = int(width / units_per_pixel)
        resolution_y = int(height / units_per_pixel)
        resolution = (resolution_x, resolution_y)
        scale_x = resolution_x / width
        scale_y = resolution_y / height

        sheetlet_image = Image.new('1', resolution, 0)
        cells_image = Image.new('1', resolution, 0)

        def invert_y_coordinates(pixels, height):
            return [(x, height - y) for x, y in pixels]

        sheetlet_pixels = polygon_to_pixels(sheetlet, scale_x, scale_y, minx, miny, resolution_x, resolution_y)
        sheetlet_pixels = invert_y_coordinates(sheetlet_pixels, resolution_y)
        draw = ImageDraw.Draw(sheetlet_image)
        draw.polygon(sheetlet_pixels, outline=0, fill=1)

        draw = ImageDraw.Draw(cells_image)
        for cell in inner_cells:
            cell_pixels = polygon_to_pixels(cell, scale_x, scale_y, minx, miny, resolution_x, resolution_y)
            cell_pixels = invert_y_coordinates(cell_pixels, resolution_y)
            draw.polygon(cell_pixels, outline=0, fill=1)

        sheetlet_image.save(os.path.join(CONFIG['images_dir'], 
                                       get_output_filename(filename, self.index, seed, 'png', 'sheetlet')))
        cells_image.save(os.path.join(CONFIG['images_dir'], 
                                    get_output_filename(filename, self.index, seed, 'png', 'cells')))
        
        sheetlet_array = np.array(sheetlet_image)
        cells_array = np.array(cells_image)

        savemat(os.path.join(CONFIG['matlab_inputs_dir'], 
                           get_output_filename(filename, self.index, seed, 'mat', 'sheetlet')), 
                {'sheetlet_image': sheetlet_array})
        savemat(os.path.join(CONFIG['matlab_inputs_dir'], 
                           get_output_filename(filename, self.index, seed, 'mat', 'cells')), 
                {'cells_image': cells_array})
    
    def load_npy(self, path):
        """Load cell data from an .npy file.
        
        Args:
            path: Path to the .npy file
        """
        cells, sheetlet, centroid = np.load(path, allow_pickle=True)
        new_cells = solve_overlaps(cells, 10/(1000*1000), (250)/(1000*1000))
        
        for i in range(len(new_cells)):
            new_cells[i] = new_cells[i].simplify(0.5/1000, preserve_topology=True)

        self.generate_matlab_images(new_cells, sheetlet, units_per_pixel=0.06/1000)