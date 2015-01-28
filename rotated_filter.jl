using Permeability
using FiltEST_VTI
using MixedVoxels
import ProgressBar
import DataFrames

include("src/utils.jl")

# config-file as command line parameter
if length(ARGS) == 0
	config = "config.jl"
else
	config = ARGS[1]
end

# read config into global variables
include("config.jl")

if plot
	println("Load Gadfly!")
	import Gadfly
end

# DataFrames table for results
table = DataFrames.DataFrame(c_Phi_0_ = typeof(Phi_0__[1])[],
				c_phi = typeof(phi_[1])[],
				c_theta = typeof(theta_[1])[],
				c_dVoxel = typeof(dVoxel_[1])[],
				c_method = ASCIIString[],
				c_pressuredrop = Float64[])
# TODO: runtime
# TODO: know table size before: more efficient?

# number of iterations
n = length(Phi_0__) * length(phi_) * length(theta_) * length(dVoxel_) * sum(a -> length(a), [method_[b] for b in dVoxel_])
println("Ready to calculate $n settings!")

# initialize progress bar
# minumum update interval: 1 second
progressbar = ProgressBar.Progress(n, 1, "Simulate rotated filter... ", 30)

# do calculations for all combinations of parameters
for Phi_0_ in Phi_0__, phi in phi_
	for theta in theta_, dVoxel in dVoxel_, method in method_[dVoxel]

		# TODO: display current parameters, needs change of ProgressMeter.jl

		# convert angles to rad
		theta = deg2rad(theta)
		phi = deg2rad(phi)

		# set method for treatment of mixed voxels
		mixedvoxel = mixedvoxelmethod[method]

		# generate housing.vti
		
		# dimensions in voxels (without in-/outlet, wall)
		Nx = round(Int64, dX/dVoxel)
		Ny = round(Int64, dY/dVoxel)
		Nz = round(Int64, dZ/dVoxel)

		# extensions of housing (with without in-/outlet, wall)
		Nyt = Ny + 2 * NWall
		Nzt = Nz + 2 * NWall

		# add layers for inlet and outlet in flow direction
		Nxt = Nx + 2 * NWall + NInlet + NOutlet


		# create FiltEST-VTI-File
		housing = FiltEST_VTIFile()

		# set origin in lower back left corner of fluid area (after Inlet)
		housing.origin = [-dVoxel*(NWall+NInlet), -dVoxel*NWall, -dVoxel*NWall]

		# set voxel length
		housing.spacing = dVoxel

		# creta data arrays
		material = DataArray(Material, 1, Base.fill(Mt["Wall"], Nxt, Nyt, Nzt))

		# TODO: Float64 vector, 6 Komponenten, aber in *.vti nur skalar? Nein
		permeability = DataArray(Float64, 1, Base.fill(1.0, Nxt, Nyt, Nzt))

		# add fluid area
		material.data[NWall+1:Nxt-NWall, NWall+1:Nyt-NWall, NWall:Nzt-NWall] = Mt["Fluid"]

		# add inlet
		material.data[:, :, 1:NInlet] = Mt["Inflow"]

		# add outlet
		material.data[:, :, end-NOutlet:end] = Mt["Outflow"]

		# filter is limited by two planes

		# center/pivot filter medium
		pivot = 0.5 * [dX, dY, dZ]

		# normal vector of planes (||.||=1)
		n_p = [-cos(phi)*sin(theta), cos(theta), sin(phi)*sin(theta)]

		# position vectors of planes
		m1 = pivot + dFilter * 0.5 * n_p
		m2 = pivot - dFilter * 0.5 * n_p

		# distances of planes from origin
		d1 = dot(m1, n_p)
		d2 = dot(m2, n_p)

		# determine material and permeability for each interior voxel
		# TODO: Reihenfolge?
		# cf. origin
		for z in 1:Nz, y in 1:Ny, x in 1:Nx
			# origin of voxel
			x0 = [x-1, y-1, z-1] * dVoxel

			# center of voxel
			center = ([x-1, y-1, z-1] .+ 0.5) * dVoxel

			# transform plane into voxel coordinate system (einheitsvoxel)
			# distance between voxel origin and plane, divided by voxel length
			dist1 = (dot(x0, n_p) - d1) / dVoxel
			dist2 = (dot(x0, n_p) - d2) / dVoxel

			# intersect Voxel with first and second plane
			for d in (dist1, dist2)
				# "-" because plane is moved against normal direction
				intersectionpoints = intersection_points(n_p, -d)

				# treat mixed voxels
				if length(intersectionpoints) > 0
					# calculate volume fraction
					xi = 1 - volumefraction(intersectionpoints, n_p, -d)

					if xi > 1-eps
						# fluid voxel
						material.data[x,y,z] = Mt["Fluid"]
						continue
					elseif xi > eps
						# mixed porous voxel
						m, p = mixedvoxel(xi)
						material.data[x,y,z] = m
						permeability.data[x,y,z] = p
						continue
					else
						# full porous voxel
						material.data[x,y,z] = Mt["Porous"]
						permeability.data[x,y,z] = K_0
						continue
					end
				end
			end

			# voxel is entirely contained in filter medium: full porous voxel
			if dot(center, n_p) - d1 < 0 && dot(center, n_p) - d2 > 0
				material.data[x,y,z] = Mt["Porous"]
				permeability.data[x,y,z] = K_0
				continue
			end

			# else voxel is fluid
			material.data[x,y,z] = Mt["Fluid"]
		end

		# add data to FiltEST-VTI-File
		add_data(housing, "Material", material)
		add_data(housing, "Permeability", permeability)


		# write FiltEST-VTI-File
		if vtifilename != ""
			write_file(housing, vtifilename, zip)
		end
		
		# write SOL-File
		if solfilename != ""
			# create file
			solfile = open(solfilename, "w")

			# write only permeability of porous voxels
			write(solfile, filter(x -> x < 1.0, permeability.data))

			# close file
			close(solfile)
		end

		# do simulation in FiltEST
		if runFiltEST == false
			ProgressBar.next!(progressbar)
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
			Progressbar.next!(progressbar, "Error in $Phi_0_ $phi $theta $dVoxel $method")

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

		# TODO: extract runtime
		# TODO: runtime of rotated_filter for that model
		# cat flowSummary.fest  | grep "Total Solver runtime" >> ../results.txt
		
		# # tidy up
		# rm "flowSummary.fest"

		# save results into table (as row)
		# TODO: runtime
		DataFrames.push!(table, [Phi_0_ phi theta dVoxel method pressuredrop])


		# update progress bar
		ProgressBar.next!(progressbar)

	end # inner for

	# TODO: create plots, call maybe in background, so that we can already continue
	if plot == true
		# TODO: subset in one step
		subset = table[table[:Phi_0_] .== Phi_0_, :]
		subset = subset[subset[:phi] .== phi, :]
		Gadfly.plot(subset, {:x => "theta", :y => "pressuredrop"})
		# TODO: save
	end

end # outer for

# write table
if writetable == true && tablefilename != ""
	tablefile = open(tablefilename, "w")
	# strip first line "10x2 DataFrame"
	write(tablefile, replace(string(table),r"^[^\n]+\n","",1))
	close(tablefile)
end
#tableFileHandle.write("""\n| theta   | phi     | Mode    | Pressure Drop (num.) | Solver Runtime (User) |
#| ------- | ------- | ------- | -------------------- | --------------------- |
#""")
	#ableFileHandle.write("| {0:<7} | {1:<7} | {2:<7} | {3:<20} | {4:<21} |\n".format(theta, phi, matchObj.group(1), numerical, matchObj.group(5)+" s"))



# Finished
println("\nFinished successfully!")
# show log