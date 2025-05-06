import subprocess
import os

def call_matlab_script(script_path, *args):
    try:
        # Get the directory and function name from the script path
        folder_path = os.path.dirname(script_path)
        function_name = os.path.splitext(os.path.basename(script_path))[0]

        # Define the path to MATLAB's bin directory
        matlab_bin_path = '/Applications/MATLAB_R2023a.app/bin/'

        # Convert each argument to a string and format them correctly for MATLAB
        formatted_args = ', '.join(f"'{arg}'" if isinstance(arg, str) else str(arg) for arg in args)
        
        # Construct the MATLAB command to change directory, call function, and quit
        matlab_command = f"cd('{folder_path}'); {function_name}({formatted_args}); quit;"
        
        # Print the actual command to be executed for debugging
        print("Executing MATLAB command:", matlab_command)
        
        # Full command to run MATLAB with the specified script
        command_run = f"./matlab -nodisplay -nosplash -r \"{matlab_command}\""
        
        # Execute the MATLAB command in the specified directory
        output_run = subprocess.check_output(command_run, shell=True, cwd=matlab_bin_path, text=True)
        print("MATLAB script executed successfully. Output:")
        print(output_run)
               
    except subprocess.CalledProcessError as e:
        print("Failed to execute MATLAB script:")
        print(e.output)
    except Exception as e:
        print("An error occurred:", str(e))


if __name__ == "__main__":
    # Specify the path to your MATLAB script
    matlab_script_path = os.path.abspath('matlab_scripts/grow_cells_uniform.m')

    #target ECVs for each sheetlet (this is to match the real cell data)
    targetECVs = [0.125349613088517, 0.0952281420263555, 0.145679632336482, 0.126471095960046, 0.125071471442900, 0.102200599822306, 0.0898216394190672, 0.125595408290213, 0.131664551644232, 0.129589715653970, 0.128288007746796, 0.110550719827034]

    """input1 = f"../matlab_inputs/cells_newapproach4_0.mat"
    input2 = f"../matlab_inputs/sheetlet_newapproach4_0.mat"
    input3 = 0.06/1000 
    input4 = 0.06/1000
    input5 = 0.125"""

    #call_matlab_script(matlab_script_path, input1, input2, input3, input4, input5,"heh")
    import concurrent.futures

    # Define a function to call the MATLAB script with the input
    def call_matlab_parallel(matlab_script_path, input1, input2, input3, input4, input5, index):
        call_matlab_script(matlab_script_path, input1, input2, input3, input4, input5, index)

    # Create a ThreadPoolExecutor with the desired number of workers
    seeds= [10]

    for seed in seeds:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            # Iterate over the range of 12
            for index in [0]:
                input1 = f"../matlab_inputs/output_sheetlet_index_{index}_seed_{seed}_cells.mat"
                input2 = f"../matlab_inputs/output_sheetlet_index_{index}_seed_{seed}_sheetlet.mat"
                input3 = 0.06/1000
                input4 = 0.06/1000
                input5 = targetECVs[index]
                output_filename = f"output_sheetlet_index_{index}_seed_{seed}"
                # Submit the call_matlab_parallel function to the executor
                executor.submit(call_matlab_parallel, matlab_script_path, input1, input2, input3, input4, input5, output_filename)
