# Sheetlet Packer

This project provides tools for packing cells within sheetlets, processing them through MATLAB, and generating meshes for simulation.

## Project Structure

- `sheetlet_packer.py`: Main class for cell packing
- `run_simulation.py`: Script to run the packing process
- `calling_matlab.py`: Interface for MATLAB processing
- `create_mesh.py`: Mesh generation from processed geometries
- `matlab_scripts/`: MATLAB scripts for cell processing
- `output/`: Directory for output files
- `images/`: Directory for visualization files
- `matlab_inputs/`: Directory for MATLAB input files
- `matlab_outputs/`: Directory for MATLAB output files
- `meshing/`: Directory for generated meshes

## Setup

### Dependencies
- Python 3.x
- Required Python packages:
  - numpy
  - shapely
  - scipy
  - matplotlib
  - PIL
  - polypacker
- MATLAB R2023a or later

### Configuration
The main configuration is defined in `sheetlet_packer.py`:

```python
CONFIG = {
    'output_dir': 'output/',
    'images_dir': 'images/',
    'matlab_inputs_dir': 'matlab_inputs/',
    'matlab_scripts_dir': 'matlab_scripts/',
    'plot_settings': {
        'dpi': 300,
        'format': 'png',
        'transparent': True
    },
    'colors': {
        'cells': '#1f77b4',
        'fillers': '#ff7f0e',
        'sheetlet': '#2ca02c'
    }
}
```

## File Outputs

### From sheetlet_packer.py
- `output/` directory:
  - `output_sheetlet_index_{index}_seed_{seed}.npy`: Initial cell arrangement
  - `lessfiller_checked_output_{seed}__{index}.npy`: Processed cell arrangement
- `images/` directory:
  - Various visualization files showing the packing progress
- `matlab_inputs/` directory:
  - Input files for MATLAB processing

### From MATLAB Processing
- `matlab_outputs/` directory:
  - `checked_geo_{seed}_{index}.mat`: Final processed cell geometries
  - `beforethicken_output_sheetlet_index_{index}_seed_{seed}.mat`: Pre-thickening state
  - `output_sheetlet_index_{index}_seed_{seed}.mat`: Post-thickening state

### From Mesh Generation
- `meshing/packed_geo/` directory:
  - `myocyte_{index}.msh`: Generated mesh files for each cell

## Usage

1. Run the packing process:
   ```bash
   python run_simulation.py
   ```

2. Process with MATLAB:
   ```bash
   python calling_matlab.py
   ```

3. Generate meshes:
   ```bash
   python create_mesh.py
   ```

## Visualization Options

The `isPlot` parameter in `sheetlet_packer.py` controls the generation of visualization files:
- When `isPlot=True` (default), additional files are generated showing the progress of the packing process and iterations
- When `isPlot=False`, only essential binary images for MATLAB processing are created

## Notes

- The process uses a seed-based system for reproducibility
- Each sheetlet is processed independently
- The MATLAB processing includes cell thickening and overlap removal
- Mesh generation creates 3D geometries from the 2D cell arrangements

## Visualization

The script generates several visualization files:
- Initial cell arrangement
- Repulsion phase results
- Final cell arrangement
- ICS (Intracellular Space) values over iterations
- Cell area distributions
- Area-to-perimeter squared ratios

All visualizations are saved in the `images/` directory with appropriate labels and scales. 