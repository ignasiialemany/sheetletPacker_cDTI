import numpy as np
import gmsh
import meshio
import trimesh

def find_border_edges(triangles):
    edge_set = set()
    for tri in triangles:
        edges = [(tri[i], tri[(i + 1) % 3]) for i in range(3)]
        for edge in edges:
            edge_sorted = tuple(sorted(edge))
            if edge_sorted in edge_set:
                edge_set.remove(edge_sorted)
            else:
                edge_set.add(edge_sorted)
    return edge_set

def remove_duplicate_vertices(mesh):
    vertices, inverse_indices = np.unique(mesh.points, axis=0, return_inverse=True)
    mapping = inverse_indices
    new_cells = []
    for cell_block in mesh.cells:
        new_cell_data = mapping[cell_block.data]
        new_cells.append(meshio.CellBlock(cell_block.type, new_cell_data))
    mesh.points = vertices
    mesh.cells = new_cells
    return mesh

class MeshGenerator:
    
    def __init__(self, vertices, mesh_length, height=0.0):
        self.vertices = [vertex for vertex in vertices if not np.isnan(vertex).any()]
        self.close_vertices()
        self.length = mesh_length
        self.height = height
    
    def close_vertices(self):
        if any(self.vertices[0] != self.vertices[-1]):
            self.vertices.append(self.vertices[0])

    def create_mesh(self, dim, filename, extrusion_h, z_offset=0.0, extrusion_layers=10, verbose=False):
        gmsh.initialize()
        gmsh.model.add("test")
        
        # Suppress GMSH output
        gmsh.option.setNumber("General.Terminal", 0)
        
        options = {
            'Geometry.Tolerance': 0.0000001,
            'Mesh.CharacteristicLengthMin': 0.0001,
            'Mesh.Algorithm': 1,
            'Mesh.RecombinationAlgorithm': 1,  # Simple algorithm
            'Mesh.ElementOrder': 1,  # First order elements
            'Mesh.Optimize': 1,  # Optimize the mesh
            'Mesh.OptimizeNetgen': 1  # Use Netgen optimizer
        }
        
        for mo, val in options.items():
            gmsh.option.setNumber(mo, val)
        
        char_length = self.length
        
        if all(self.vertices[0] == self.vertices[-1]):
            self.vertices = self.vertices[:-1]
        
        # Create points
        points = [gmsh.model.geo.addPoint(v[0], v[1], z_offset, char_length) for v in self.vertices]
        
        # Create lines with consistent orientation
        lines = []
        for j in range(len(points)):
            # Ensure counter-clockwise orientation for outward normals
            lines.append(gmsh.model.geo.addLine(points[j], points[(j + 1) % len(points)]))
        
        curve_loop = gmsh.model.geo.addCurveLoop(lines)
        plane_surface = gmsh.model.geo.addPlaneSurface([curve_loop])
        
        gmsh.model.geo.synchronize()
        
        # Extrude with consistent orientation
        extrusion = gmsh.model.geo.extrude([(2, plane_surface)], 0, 0, extrusion_h, [extrusion_layers])
        
        # Add physical groups with consistent orientation
        gmsh.model.addPhysicalGroup(2, [extrusion[0][1]], 1)  # top surface
        side_surfaces = [extrusion[i][1] for i in range(2, len(extrusion))]
        gmsh.model.addPhysicalGroup(2, side_surfaces, 2)  # side surfaces
        gmsh.model.addPhysicalGroup(2, [plane_surface], 3)  # bottom surface
        gmsh.model.addPhysicalGroup(3, [extrusion[1][1]], 4)  # volume
        
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(2)
        gmsh.write(filename)
        
        if verbose:
            print(f"Saved initial {filename}")
        
        gmsh.model.remove()
        gmsh.finalize()
        
        mesh = meshio.read(filename)
        mesh = remove_duplicate_vertices(mesh)
        
        filename_base = filename.split(".")[0]
        if "triangle" in mesh.cells_dict:
            triangles = mesh.cells_dict["triangle"]
            border_edges = find_border_edges(triangles)
            
            if border_edges:
                if verbose:
                    print(f"Mesh has border edges: {border_edges}")
                return RuntimeError
            else:
                mesh.write(filename_base + ".ply", file_format="ply")
                mesh_check = trimesh.load(filename_base + ".ply")
                
                # Check if the mesh is watertight
                if not mesh_check.is_watertight:
                    if verbose:
                        print("Mesh is not watertight.")
                    return RuntimeError("Mesh is not watertight.")
                
                # Ensure normals are consistent and pointing outward
                if not mesh_check.is_winding_consistent:
                    if verbose:
                        print("Fixing winding consistency...")
                    mesh_check.fix_normals()
                
                # Make sure volume is positive (indicates outward normals)
                if mesh_check.volume < 0:
                    if verbose:
                        print("Inverting mesh to ensure outward normals...")
                    mesh_check.invert()
                
                # Additional checks and fixing for normals
                mesh_check.remove_degenerate_faces()
                mesh_check.remove_duplicate_faces()
                
                # Verify normal orientation
                if mesh_check.volume < 0:
                    if verbose:
                        print("Warning: Mesh volume is still negative after fixing")
                    mesh_check.invert()  # Try one more time
                
                # Check the watertightness and consistency again after fixing
                if not mesh_check.is_watertight:
                    if verbose:
                        print("Mesh is not watertight after fixing.")
                    return RuntimeError("Mesh is not watertight after fixing.")
                
                if not mesh_check.is_winding_consistent:
                    if verbose:
                        print("Mesh winding is not consistent after fixing.")
                    return RuntimeError("Mesh winding is not consistent after fixing.")
                
                # Final verification of normal orientation
                if mesh_check.volume < 0:
                    if verbose:
                        print("Error: Could not fix normal orientation")
                    return RuntimeError("Could not fix normal orientation")
    
                # Export to other formats
                mesh_check.export(filename_base + ".off")
                if verbose:
                    print(f"Mesh successfully exported to {filename_base}.off with outward-pointing normals")
                return True
        else:
            if verbose:
                print("Mesh does not contain triangular faces.")
            return None

