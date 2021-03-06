using FiltEST_VTI
using MixedVoxels
using SaveK
import ProgressBar
import DataFrames

include("src/Permeability.jl")
include("src/utils.jl")

# read config into global variables
if !isdefined(:config)
	include("config_normal_filter.jl")
end

# DataFrames table for results
table = DataFrames.DataFrame(Phi_0_ = typeof(𝚽_0_[1])[],
				h = typeof(dVoxel[1])[],
				method = String[],
				xi = typeof(𝜉[1])[],
				pressuredrop = Float64[],
				pressuredrop_ideal = Float64[],
				relerror = Float64[],
				runtime = Float16[],
				FiltEST_runtime = Float16[])

# calculate number of iterations
const n = length(𝚽_0_) * length(𝜉) * sum(a -> length(a), [methods[b] for b in dVoxel])
message = ("Ready to calculate $n settings!")
println(message)

# initialize progress bar
# minumum update interval: 1 second
progressbar = ProgressBar.Progress(n+2, 1, "Simulate normal filter... ", 30)

# do calculations for all combinations of parameters
for Phi_0_ in 𝚽_0_, h in dVoxel, method in methods[h], xi in 𝜉

	# update progress bar
	message = "$Phi_0_ $h $method $xi"
	ProgressBar.next!(progressbar, message)
	message = ""

	# take time of execution
	tic()

	# set method for treatment of mixed voxels
	function mixedvoxel(xi)
		if xi > 1-𝜀
			return Mt["Fluid"], 1.0
		elseif xi > 𝜀
			return mixedvoxelmethod[method](xi, Phi_0_)
		else
			return Mt["Porous"], K_0
		end
	end

	# calculate ideal exact pressure drop
	# velocity in mm/s
	u = 0.0384*1e6/60/(L_y * L_z)
	# pressure drop in kPa given by theoretical considerations
	pressuredrop_ideal = 𝜇 * 1e-3 * u/K_0 * (dFilter + (1 - xi) * h)

	# generate housing.vti
	
	# dimensions in voxels (without in-/outlet, wall)
	Nx = round(Int64, L_x/h)
	Ny = round(Int64, L_y/h)
	Nz = round(Int64, L_z/h)

	# extensions of housing (with without in-/outlet, wall)
	Nyt = Ny + 2 * NWall
	Nzt = Nz + 2 * NWall

	# add layers for inlet and outlet in flow direction
	Nxt = Nx + 2 * NWall + NInlet + NOutlet


	# create FiltEST-VTI-File
	housing = FiltEST_VTIFile()

	# set geometry origin in lower back left corner of fluid area (after Inlet)
	NOffOrigin = [-(NWall+NInlet), -NWall, -NWall]
	housing.origin = h * NOffOrigin

	# set voxel length
	housing.spacing = h

	# create data arrays
	material = DataArray(Material, 1, Base.fill(Mt["Symmetry"], Nxt, Nyt, Nzt))

	# TODO: Float64 vector, 6 Komponenten, aber in *.vti nur skalar? Nein
	# permeability of wall =1, because not porous
	permeability = DataArray(Float64, 1, Base.fill(1.0, Nxt, Nyt, Nzt))

	# add fluid domain
	material.data[NWall+1:Nxt-NWall, NWall+1:Nyt-NWall, NWall+1:Nzt-NWall] = Mt["Fluid"]

	# add inlet
	material.data[NWall+1:NWall+NInlet, NWall+1:Nyt-NWall, NWall+1:Nzt-NWall] = Mt["Inflow"]

	# add outlet
	material.data[Nxt-NOutlet-NWall+1:Nxt-NWall, NWall+1:Nyt-NWall, NWall+1:Nzt-NWall] = Mt["Outflow"]

	# add data to FiltEST-VTI-File
	add_data(housing, "Material", material)
	add_data(housing, "Permeability", permeability)

	# filter is limited by two planes
	geom = Geometry([Nx,Ny,Nz], h)

	# center/pivot filter medium
	pivot = 0.5 * [L_x, L_y, L_z]

	# normal vector of planes (||.|| = 1)
	n_p = [1, 0, 0]

	# position vectors of planes
	# distances of planes from origin
	d1 = pivot[1] - dFilter * 0.5
	d2 = pivot[1] + dFilter * 0.5 + h * (1-xi)

	# add planes to geometry
	plane1 = Plane(-n_p, -d1)
	plane2 = Plane(n_p, d2)

	add!(geom, plane1)
	add!(geom, plane2)

	# generate voxels
	discretise(geom, NOffOrigin, mixedvoxel,
		housing.voxeldata["Material"],
		housing.voxeldata["Permeability"])


	# write FiltEST-VTI-File
	if vtifilename != ""
		write_file(housing, vtifilename, zipVTI)
	end

	# write SOL-File
	if solfilename != ""
		saveK(permeability.data, solfilename)
	end

	# runtime of rotated_filter for that model
	runtime = float64(toq())


	# do simulation in FiltEST
	if runFiltEST == false
		# continue with generation of VTI for next parameter set
		continue
	end

	# redirect STDOUT and STDERR to capture FiltEST output
	stdout = STDOUT
	stderr = STDERR
	(outread, outwrite) = redirect_stdout()
	(errorread, errorwrite) = redirect_stderr()

	# call FiltEST
	try
		run(`$filtest flowSimInput.xml`)
	catch
		# something went wrong when running FiltEST

		# get error output
		close(outwrite)
		close(errorwrite)
		filtest_error = readavailable(errorread)

		# tidy up
		close(outread)
		close(errorread)
		redirect_stdout(stdout)
		redirect_stderr(stderr)

		# log error
		message = "Error in $Phi_0_ $h $method $xi"
		ProgressBar.logerror(progressbar, message)

		# extract first error
		m = match(r".*ERROR: (.+)\n", filtest_output)
		if m != nothing
			ProgressBar.logerror(progressbar, m.captures[1])
		end
		m = match(r".*ERROR: (.+)\n", filtest_error)
		if m != nothing
			ProgressBar.logerror(progressbar, m.captures[1])
		end

		# try next parameter set
		continue
	end

	# stop capturing of STDOUT and STDERR
	# ensures that any buffered writes are flushed and available for reading
	close(outwrite)
	close(errorwrite)

	# capture STDOUT and STDERR
	# TODO: readavailable/readall?
	filtest_output = readavailable(outread)
	filtest_error = readavailable(errorread)

	close(outread)
	close(errorread)

	# restore original STDOUT and STDERR
	redirect_stdout(stdout)
	redirect_stderr(stderr)

	# extract pressure drop
	m = match(r"PressureDrop: ([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)", filtest_output)
	pressuredrop = float64(m.captures[1])

	# extract runtime of FiltEST
	m = match(r"Total Solver runtime: ([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?) \(user clock\)", filtest_output)
	FiltEST_runtime = float64(m.captures[1])
	
	relerror = abs((pressuredrop - pressuredrop_ideal)/pressuredrop_ideal)

	# save results into table (as row)
	DataFrames.push!(table, [Phi_0_ h method xi pressuredrop pressuredrop_ideal relerror runtime FiltEST_runtime])

end # end for

ProgressBar.next!(progressbar, "Finished calculations. Write table.")

# write table to file
if writetable == true && tablefilename != ""

	# open file
	tablefile = open(tablefilename, "w")

	# write header
	write(tablefile, "Normal Filter (dFilter = $dFilter, K_0 = $K_0) for fluid with viscosity = $𝜇, density = $𝜌\n" * "="^20 * "\n")

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
