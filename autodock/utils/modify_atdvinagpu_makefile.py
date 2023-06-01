"""
This file is used when installing Vina-GPU.
It modifies the makefile to compile Vina-GPU with the necessary values.
"""

import sys

def modifyATDVinaGPUMakefile(filePath, boostLibPath, openCLLibPath, openCLVersion, gpuPlatform):
    """ This function receives a makefile and the 4 values neccessary for compiling AutoDock-VinaGPU. """
    # Reading makefile
    with open(filePath, 'r') as f:
        lines = f.readlines()
    
    # Defining the new values for the variables
    newValues = {
        'BOOST_LIB_PATH': boostLibPath,
        'OPENCL_LIB_PATH': openCLLibPath,
        'OPENCL_VERSION': openCLVersion,
        'GPU_PLATFORM': gpuPlatform
    }

    # Modifying the lines
    for i, line in enumerate(lines):
        if line.startswith(tuple(newValues.keys())):
            key = line.split('=')[0]
            lines[i] = f'{key}={newValues[key]}\n'

    # Writing the modified lines back to the file
    with open(filePath, 'w') as f:
        f.writelines(lines)

if __name__ == "__main__":
    """ This is the main execution funcion. Captures the arguments and calls the makefile modification function. """
    # Checking param length
    if len(sys.argv) != 6:
        print(f"Usage: {sys.argv[0]} filePath boostLibPath openCLLibPath openCLVersion gpuPlatform")
        sys.exit(1)

    # Calling modification function
    modifyATDVinaGPUMakefile(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])