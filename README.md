# Sheetlet Packer

This project generates and processes cell geometries for sheetlet packing simulations. It creates a grid of cells and fillers within a sheetlet boundary, processes them through MATLAB for thickening, and generates meshes for further analysis.

## Project Structure

```
.
├── sheetlet_packer.py      # Main cell packing and processing script
├── calling_matlab.py       # Interface for MATLAB processing
├── create_mesh.py          # Mesh generation from processed cells
├── auxiliary.py            # Helper functions
├── matlab_scripts/         # MATLAB processing scripts
├── output/                 # Output .npy files
├── images/                 # Visualization outputs
├── matlab_inputs/          # Input files for MATLAB
├── matlab_outputs/         # Output files from MATLAB
└── meshing/                # Generated meshes
```

## Setup

1. **Dependencies**:
   - Python 3.x
   - MATLAB R2023a or later
   - Required Python packages:
     - numpy
     - shapely
     - scipy
     - matplotlib
     - PIL
     - polypacker

2. **MATLAB Setup**:
   - Ensure MATLAB is installed and accessible from the command line
   - MATLAB scripts should be in the `matlab_scripts/` directory

3. **Directory Structure**:
   - The script will automatically create necessary directories:
     - `output/`
     - `images/`
     - `matlab_inputs/`
     - `matlab_outputs/`

## Configuration

The main configuration is in `sheetlet_packer.py`:

```python
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
```

## Visualization Options

The `isPlot` parameter controls the generation of visualization files. When `isPlot=True`, the following images are generated:

### During Initial Packing (`sheetlet_packer.py`):
- **Initial Phase**:
  - `cellgrid_output_sheetlet_index_{index}_seed_{seed}_init.png`: Shows the initial cell arrangement with:
    - Cells in red (#8a0101)
    - Fillers in black
    - Sheetlet boundary in green (#00d149)
    - Includes legend and axis labels in millimeters

- **Repulsion Phase**:
  - `cellgrid_output_sheetlet_index_{index}_seed_{seed}_repulsion.png`: Shows the cell arrangement after the initial repulsion phase
  - Generated every 5 steps during the repulsion phase

- **Main Packing Phase**:
  - `cellgrid_output_sheetlet_index_{index}_seed_{seed}_step_{i}.png`: Shows the cell arrangement during the main packing phase
  - Generated every 100 steps
  - Includes cells, fillers, and sheetlet boundary

- **Final Results**:
  - `cellgrid_output_sheetlet_index_{index}_seed_{seed}_final.png`: Shows the final cell arrangement
  - `cellgrid_output_sheetlet_index_{index}_seed_{seed}_ics.png`: Plots the Intracellular Space (ICS) values over iterations
  - `cellgrid_output_sheetlet_index_{index}_seed_{seed}_sheetlet.png`: Binary image of the sheetlet boundary
  - `cellgrid_output_sheetlet_index_{index}_seed_{seed}_cells.png`: Binary image of the cell arrangement

## File Outputs for sheetlet_packer.py

### 1. Initial Packing (`sheetlet_packer.py`)
- **Output Directory**: `output/`
  - `output_sheetlet_index_{index}_seed_{seed}.npy`: Contains packed cells, sheetlet, and centroid
  - `lessfiller_checked_output_{seed}__{index}.npy`: Initial cell data before thickening

- **Images Directory**: `images/`
  - `cellgrid_output_sheetlet_index_{index}_seed_{seed}_init.png`: Initial cell arrangement
  - `cellgrid_output_sheetlet_index_{index}_seed_{seed}_repulsion.png`: After repulsion phase
  - `cellgrid_output_sheetlet_index_{index}_seed_{seed}_final.png`: Final cell arrangement
  - `cellgrid_output_sheetlet_index_{index}_seed_{seed}_ics.png`: ICS values plot

- **MATLAB Inputs**: `matlab_inputs/`
  - `output_sheetlet_index_{index}_seed_{seed}_cells.mat`: Cell data for MATLAB
  - `output_sheetlet_index_{index}_seed_{seed}_sheetlet.mat`: Sheetlet boundary for MATLAB

### 2. MATLAB Processing (`calling_matlab.py`)
- **MATLAB Outputs**: `matlab_outputs/`
  - `checked_geo_{seed}_{index}.mat`: Processed cell geometries
  - `beforethicken_output_sheetlet_index_{index}_seed_{seed}.mat`: Cells before thickening
  - `output_sheetlet_index_{index}_seed_{seed}.mat`: Final thickened cells

### 3. Mesh Generation (`create_mesh.py`)
- **Mesh Outputs**: `meshing/packed_geo/`
  - `myocyte_{counter}.msh`: Individual cell meshes
- **Analysis Files**:
  - `areas_nothicken_{seed}.npy`: Cell areas before thickening
  - `check_thicken_{seed}.npy`: Cell areas after thickening

## Usage

1. **Run Packing**:
```python
python sheetlet_packer.py
```

2. **Process with MATLAB**:
```python
python calling_matlab.py
```

3. **Generate Meshes**:
```python
python create_mesh.py
```

## Notes

- The script uses a seed-based system for reproducibility
- Each sheetlet (index) is processed independently
- The process includes:
  1. Initial cell packing
  2. Repulsion phase
  3. MATLAB processing for thickening
  4. Mesh generation
  5. Analysis of cell properties

## Visualization

The script generates several visualization files:
- Initial cell arrangement
- Repulsion phase results
- Final cell arrangement
- ICS (Intracellular Space) values over iterations
- Cell area distributions
- Area-to-perimeter squared ratios

All visualizations are saved in the `images/` directory with appropriate labels and scales. 