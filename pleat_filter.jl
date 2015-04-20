using FiltEST_VTI
using MixedVoxels
import ProgressBar
import DataFrames

include("src/Permeability.jl")
include("src/utils.jl")


# read config into global variables
if !isdefined(:config)
	include("pleat_config.jl")
end

# load Gadfly only if needed
if plot
	println("Loading ... (Gadfly)!")
	import Gadfly
end

# DataFrames table for results
table = DataFrames.DataFrame(Phi_0_ = typeof(𝚽_0_[1])[],
				h = typeof(dVoxel[1])[],
				method = String[],
				pressuredrop = Float64[],
				runtime = Float16[],
				FiltEST_runtime = Float16[])

# calculate number of iterations
const n = length(𝚽_0_) * sum(a -> length(a), [methods[b] for b in dVoxel])
message = ("Ready to calculate $n settings!")
println(message)

# initialize progress bar
# minumum update interval: 1 second
progressbar = ProgressBar.Progress(n+2, 1, "Simulate rotated filter... ", 30)

# do calculations for all combinations of parameters
for Phi_0_ in 𝚽_0_, h in dVoxel, method in methods[h]

	# update progress bar
	message = "$Phi_0_ $phi $theta $h $method"
	ProgressBar.next!(progressbar, message)
	message = ""

	# take time of execution
	tic()

	# set method for treatment of mixed voxels
	mixedvoxel = mixedvoxelmethod[method]

	# generate housing.vti

	# dimensions in voxels (without in-/outlet, wall)
	Nx = round(Int64, L_x/h)
	Ny = round(Int64, L_y/h)
	Nz = round(Int64, L_z/h)

	#

	# create FiltEST-VTI-File
	housing = FiltEST_VTIFile()

	# set voxel length
	housing.spacing = h

	# create data arrays
	# FIXME: is Wall symmetrie?
	# symmetry boundary condition
	material = DataArray(Material, 1, Base.fill(Mt["Wall"], Nxt, Nyt, Nzt))

	# TODO: Float64 vector, 6 Komponenten, aber in *.vti nur skalar? Nein
	permeability = DataArray(Float64, 1, Base.fill(1.0, Nxt, Nyt, Nzt))

	# add fluid domain

	# set origin

	# add inlet

	# add outlet


	# add data arrays to FiltEST-VTI-File
	add_data(housing, "Material", material)
	add_data(housing, "Permeability", permeability)

	# filter is a pleat
	add_geometry(housing, inner_arc)
	add_geometry(housing, outer_arc)
	add_geometry(housing, inner_plane_left)
	add_geometry(housing, outer_plane_left)
	add_geometry(housing, inner_plane_right)
	add_geometry(housing, outer_plane_right)

	# do voxelisation
	generate(housing, method)

	# write FiltEST-VTI-File
	if vtifilename != ""
		write_file(housing, vtifilename, zipVTI)
	end

	# write SOL-File

	# runtime of rotated_filter for that model
	runtime = float64(toq())


	# do simulation in FiltEST
	if runFiltEST == false
		continue
	end

	# redirect STDOUT and STDERR to capture FiltEST output
	stdout = STDOUT
	stderr = STDERR
	(outread, outwrite) = redirect_stdout()
	(errorread, errorwrite) = redirect_stderr()

	# call FiltEST
	# TODO: call FiltEST directly without bash script, problem: Umgebungsvariablen für Pfade
	try
		run(`$filtest flowSimInput.xml`)
	catch
		# something went wrong when running FiltEST

		# get error output
		filtest_error = readavailable(errorread)

		# log error
		message = "Error in $Phi_0_ $h $method"

		# try next parameter set
		continue
	end

	# capture STDOUT and STDERR
	filtest_output = readavailable(outread)
	filtest_error = readavailable(errorread)

	# stop capturing of STDOUT and STDERR
	close(outwrite)
	close(errorwrite)
	close(outread)
	close(errorread)

	# restore original STDOUT and STDERR
	redirect_stdout(stdout)
	redirect_stderr(stderr)

	# extract pressure drop
	m = match(r"PressureDrop: ([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)", filtest_output)
	pressuredrop = float64(m.captures[1])

	# TODO: extract runtime of FiltEST
	# cat flowSummary.fest  | grep "Total Solver runtime" >> ../results.txt
	m = match(r"Total Solver runtime: ([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?) \(user clock\)", filtest_output)
	FiltEST_runtime = float64(m.captures[1])
	
	# save results into table (as row)
	DataFrames.push!(table, [Phi_0_ h method pressuredrop runtime FiltEST_runtime])

end

# plot

ProgressBar.next!(progressbar, "Finished calculations. Write table.")

# write table to file
if writetable == true && tablefilename != ""

	# open file
	tablefile = open(tablefilename, "w")

	# write header
	write(tablefile, "Pleated Filter (dFilter = $dFilter, K_0 = $K_0) for fluid with viscosity = $𝜇, density = $𝜌\n" * "="^20)

	# strip first line "10x2 DataFrame"
	#write(tablefile, replace(string(table),r"^[^\n]+\n","",1))
	# write full table in Markdown style, not only chunks, without header
	DataFrames.showall(tablefile, table, false, symbol("Row"), false)

	# close file
	close(tablefile)
end

# write table to CSV-file. Can easily be reloaded
if writeCSV == true && csvfilename != ""
	DataFrames.writetable(csvfilename, table)
end

# Finished
ProgressBar.next!(progressbar, "Finished successfully!")