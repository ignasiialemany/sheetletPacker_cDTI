import multiprocessing
import numpy as np
from sheetlet_packer import SheetletPacker

# Configuration for the main execution
PROCESS_CONFIG = {
    'indices': [0],  # List of sheetlet indices to process
    'seeds': [10],   # List of random seeds to use
    'input_file_pattern': 'sheetlet_{}.npy',  # Pattern for input files
    'output_prefix': 'output',  # Prefix for output files
    'use_multiprocessing': False,  # Whether to use parallel processing
    'plot_results': False  # Whether to generate plots
}

def process_sheetlet(i, seed, input_pattern, output_prefix, plot_results):
    """Process a single sheetlet with the given seed.
    
    Args:
        i: Sheetlet index
        seed: Random seed
        input_pattern: Pattern for input files
        output_prefix: Prefix for output files
        plot_results: Whether to generate plots
    """
    np.random.seed(i + seed)
    sheetlet_packer = SheetletPacker(
        input_pattern.format(i),
        sheetletIndex=i,
        celltype="real"
    )
    sheetlet_packer.pack(output_prefix, seed=seed, isPlot=plot_results)

def main():
    """Main function to run the simulation."""
    if PROCESS_CONFIG['use_multiprocessing']:
        processes = []
        for seed in PROCESS_CONFIG['seeds']:
            for i in PROCESS_CONFIG['indices']:
                p = multiprocessing.Process(
                    target=process_sheetlet,
                    args=(
                        i,
                        seed,
                        PROCESS_CONFIG['input_file_pattern'],
                        PROCESS_CONFIG['output_prefix'],
                        PROCESS_CONFIG['plot_results']
                    )
                )
                processes.append(p)
                p.start()

            for p in processes:
                p.join()
    else:
        for i in PROCESS_CONFIG['indices']:
            for seed in PROCESS_CONFIG['seeds']:
                process_sheetlet(
                    i,
                    seed,
                    PROCESS_CONFIG['input_file_pattern'],
                    PROCESS_CONFIG['output_prefix'],
                    PROCESS_CONFIG['plot_results']
                )

if __name__ == '__main__':
    main() 