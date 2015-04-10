using FiltEST_VTI
using MixedVoxels
import ProgressBar
import DataFrames

include("src/Permeability.jl")
include("src/utils.jl")


# read config into global variables
include("pleat_config.jl")

# load Gadfly only if needed
if plot
	println("Loading ... (Gadfly)!")
	import Gadfly
end

for Phi_0_ in 𝚽_0_, h in dVoxel, method in methods[h]

	# update progress bar
	message = "$Phi_0_ $phi $theta $h $method"
	ProgressBar.next!(progressbar, message)
	message = ""

	# take time of execution
	tic()

	# set method for treatment of mixed voxels
	mixedvoxel = mixedvoxelmethod[method]

	# dimensions in voxels (without in-/outlet, wall)
	Nx = round(Int64, L_x/h)
	Ny = round(Int64, L_y/h)
	Nz = round(Int64, L_z/h)

	#

	# create FiltEST-VTI-File
	housing = FiltEST_VTIFile()

	# set voxel length
	housing.spacing = h

	# set origin

	# add inlet

	# add outlet

	# symmetry boundary condition

	# 
	add_geometry(housing, inner_arc)
	add_geometry(housing, outer_arc)
	add_geometry(housing, inner_plane_left)
	add_geometry(housing, outer_plane_left)
	add_geometry(housing, inner_plane_right)
	add_geometry(housing, outer_plane_right)

	generate(housing, method, "Material", "Permeability")

end

# Finished
ProgressBar.next!(progressbar, "Finished successfully!")