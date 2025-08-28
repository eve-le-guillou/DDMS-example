import subprocess
from paraview.simple import *

def read_file(input_file):
    extension = input_file.split(".")[-1]
    if extension == "raw":
        extent, dtype = input_file.split(".")[0].split("_")[-2:]
        extent = [int(dim) for dim in extent.split("x")]

        dtype_pv = {
            "uint8": "unsigned char",
            "int16": "short",
            "uint16": "unsigned short",
            "float32": "float",
            "float64": "double",
        }

        raw = ImageReader(FileNames=[input_file])
        raw.DataScalarType = dtype_pv[dtype]
        raw.DataByteOrder = "LittleEndian"
        raw.DataExtent = [0, extent[0] - 1, 0, extent[1] - 1, 0, extent[2] - 1]
        return raw

    return None

filename = "backpack_512x512x373_uint16.raw"
bashCommand = f"wget http://klacansky.com/open-scivis-datasets/backpack/{filename} -P ."

result = subprocess.Popen(bashCommand, shell=True, universal_newlines=True)
return_code = result.wait()
if return_code != 0:
    print("Error for command: "+bashCommand)

data = read_file(filename)

if "uint16" in filename:
	calculator1 = Calculator(registrationName='Calculator1', Input=data)
	calculator1.Function = ''

	# Properties modified on calculator1
	calculator1.ResultArrayName = 'ImageFile'
	calculator1.Function = 'ImageFile'
	calculator1.ResultArrayType = 'Int'
	    
	SaveData("dataset.pvd", calculator1)
else:
	SaveData("dataset.pvd", data)
