#using Permeability
using FiltEST_VTI
using MixedVoxels
import ProgressBar
import DataFrames

include("src/Permeability.jl")
include("src/utils.jl")

# config-file as command line parameter
if length(ARGS) == 0
	config = "config.jl"
else
	config = ARGS[1]
end

# read config into global variables
include(config)

if plot
	println("Loading ... (Gadfly)!")
	import Gadfly
end

# DataFrames table for results
table = DataFrames.DataFrame(Phi_0_ = typeof(Phi_0__[1])[],
				phi = typeof(phi_[1])[],
				theta = typeof(theta_[1])[],
				dVoxel = typeof(dVoxel_[1])[],
				method = ASCIIString[],
				pressuredrop = Float64[],
				runtime = Float16[],
				FiltEST_runtime = Float16[])

# TODO: know table size before: more efficient?

# number of iterations
n = length(Phi_0__) * length(phi_) * length(theta_) * sum(a -> length(a), [method_[b] for b in dVoxel_])
message = ("Ready to calculate $n settings!")
println(message)

# initialize progress bar
# minumum update interval: 1 second
progressbar = ProgressBar.Progress(n+2, 1, "Simulate rotated filter... ", 30)

# do calculations for all combinations of parameters
for Phi_0_ in Phi_0__, phi in phi_
	for theta in theta_, dVoxel in dVoxel_, method in method_[dVoxel]

		# TODO: display current parameters, needs change of ProgressMeter.jl
		# update progress bar
		message = "$Phi_0_ $phi $theta $dVoxel $method"
		ProgressBar.next!(progressbar, message)
		message = ""

		# take time of execution
		tic()

		# convert angles to rad
		theta_rad = deg2rad(theta)
		phi_rad = deg2rad(phi)

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
		NOrigin = [-(NWall+NInlet), -NWall, -NWall]
		housing.origin = dVoxel * NOrigin

		# set voxel length
		housing.spacing = dVoxel

		# creta data arrays
		material = DataArray(Material, 1, Base.fill(Mt["Solid"], Nxt, Nyt, Nzt))

		# TODO: Float64 vector, 6 Komponenten, aber in *.vti nur skalar? Nein
		permeability = DataArray(Float64, 1, Base.fill(1.0, Nxt, Nyt, Nzt))

		# add fluid area
		material.data[NWall+1:Nxt-NWall, NWall+1:Nyt-NWall, NWall+1:Nzt-NWall] = Mt["Fluid"]

		# add inlet
		material.data[NWall+1:NWall+NInlet, NWall+1:Nyt-NWall, NWall+1:Nzt-NWall] = Mt["Inflow"]

		# add outlet
		material.data[Nxt-NOutlet-NWall+1:Nxt-NWall, NWall+1:Nyt-NWall, NWall+1:Nzt-NWall] = Mt["Outflow"]

		# filter is limited by two planes

		# center/pivot filter medium
		pivot = 0.5 * [dX, dY, dZ]

		# normal vector of planes (||.||=1)
		n_p = [-cos(phi_rad)*sin(theta_rad), cos(theta_rad), sin(phi_rad)*sin(theta_rad)]

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
			xyz = [x,y,z]

			# origin of voxel
			x0 = (xyz .- 1) * dVoxel

			# center of voxel
			center = (xyz .- 0.5) * dVoxel

			pos = xyz - NOrigin

			# transform plane into voxel coordinate system (einheitsvoxel)
			# distance between voxel origin and plane, divided by voxel length
			dist1 = (dot(x0, n_p) - d1) / dVoxel
			dist2 = (dot(x0, n_p) - d2) / dVoxel

			# intersect Voxel with first and second plane
			for (d, inv) in ((dist1, false), (dist2, true))

				# "-" because plane is moved against normal direction
				intersectionpoints = intersection_points(n_p, -d)

				# treat mixed voxels
				if length(intersectionpoints) > 0

					# calculate volume fraction
					if inv
						# volume on opposite normal side
						xi = volumefraction(intersectionpoints, n_p, -d)
					else
						xi = 1 - volumefraction(intersectionpoints, n_p, -d)
					end

					if xi > 1 - eps
						# fluid voxel
						material.data[pos...] = Mt["Fluid"]
						permeability.data[pos...] = 1.0
						break
					elseif xi > eps
						# mixed porous voxel
						m, p = mixedvoxel(xi, Phi_0_, center)
						material.data[pos...] = m
						permeability.data[pos...] = p
						break
					else
						# full porous voxel
						material.data[pos...] = Mt["Porous"]
						permeability.data[pos...] = K_0
						break
					end
				elseif dot(center, n_p) - d1 < 0 && dot(center, n_p) - d2 > 0
					# voxel is entirely contained in filter medium: full porous voxel
					material.data[pos...] = Mt["Porous"]
					permeability.data[pos...] = K_0
				else
					# else voxel is fluid
					material.data[pos...] = Mt["Fluid"]
					permeability.data[pos...] = 1.0
				end
			end			
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
			#write(solfile, filter(x -> x < 1.0, permeability.data))
			for x in filter(x -> x < 1.0, permeability.data)
				for i in 1:6
					write(solfile, x)
				end
			end

			# close file
			close(solfile)
		end

		# runtime of rotated_filter for that model
		runtime = float64(toq())


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
			message = "Error in $Phi_0_ $phi $theta $dVoxel $method"

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
		DataFrames.push!(table, [Phi_0_ phi theta dVoxel method pressuredrop runtime FiltEST_runtime])


	end # inner for

	# tidy up
	# FIXME: causes halt when does not exist
	#rm("flowSummary.fest")

	# TODO: create plots, call maybe in background, so that we can already continue
	if plot == true
		# TODO: subset in one step
		# subset(iris, :(Species .== "setosa"))
		#subset = table[table[:Phi_0_] .== Phi_0_, :]
		#subset = subset[subset[:phi] .== phi, :]
		# need to convert boolean DataFrames.DataArray into normal array (replacing NA with 0) in order to apply elementwise AND
		subset = table[DataFrames.array(table[:Phi_0_] .== Phi_0_, 0) & DataFrames.array(table[:phi] .== phi, 0), :]
		# TODO: plot title, legends, Syntax
		# TODO: "Splatting", in Julia.md
		#plot = Gadfly.plot([layer(y = subset[array(subset[:dVoxel] .== dVoxel, 0) & array(subset[:method] .== method, ""), :], x = theta_,
		#	Gadfly.Geom.line) for dVoxel in dVoxel_, method in method_[dVoxel]]...)
			#Theme(default_color = color(["red" "blue" "green" "cyan" "magenta" "yellow"][i%6+1]))
#p = plot([[ ]...)

		layers = Gadfly.Layer[]
		#for dVoxel in dVoxel_, method in method_[dVoxel]
		#	plotdata = subset[DataFrames.array(subset[:dVoxel] .== dVoxel, 0) & DataFrames.array(subset[:method] .== method, 0), :]
		#	push!(layers, Gadfly.layer(plotdata, x ="theta", y ="pressuredrop", Gadfly.Geom.point, Gadfly.Geom.line)[1])
		#end

		#red: #e41a1c
		#blue: #377eb8
		#green: #4daf4a
		#purple: #984ea3
		#orange: #ff7f00
		#yellow: #ffff33
		#brown: #a65628
		#pink: #f781bf

		colors = [Gadfly.color(c) for c in ["orange", "blue", "red", "purple", "green", "yellow"]]
#=plot( dataframe,
    layer( x=:X_Symbol, y=:Y_Symbol, Geom.line, color=["Series Label"],
    layer( x=:X_Symbol2, y=:Y_Symbol2, Geom.line, color=["Series Label 2"],
    Scale.discrete_color_manual(color_arr...)
)=#

		# split table into rows with same dVoxel, method
		for plotdata in DataFrames.groupby(subset, [:dVoxel, :method])
			# add plot layer
			push!(layers, Gadfly.layer(plotdata,
									x ="theta",
									y ="pressuredrop",
									Gadfly.Geom.point, Gadfly.Geom.line,
									color=["$(plotdata[:method][1]) $(plotdata[:dVoxel][1])"])...)
		end
		# plot
		p = Gadfly.plot(layers,
			Gadfly.Scale.color_discrete_manual(colors...),
			Gadfly.Guide.xlabel("𝜃"),
			Gadfly.Guide.ylabel("pressuredrop 𝚫p"),
			Gadfly.Guide.title("Rotated Filter: Permeability Scaling for mixed Voxels\n 𝜑 = $phi"))

		
		#plot = Gadfly.plot(subset, {:x => "theta", :y => "pressuredrop"},
		#	Gadfly.Guide.xlabel("𝜃"),
		#	Gadfly.Guide.ylabel("pressuredrop 𝚫p"),
		#	Gadfly.Guide.title("Rotated Filter: Permeability Scaling for mixed Voxels\n phi𝜑=$phi"))
		# save plot
		# TODO: PGF
		image = Gadfly.PDF("Phi_0_$(Phi_0_)_phi_$(phi).pdf", 12Gadfly.cm, 7.5Gadfly.cm)
		Gadfly.draw(image, p)
		#Gadfly.finish(image)
	end

end # outer for

ProgressBar.next!(progressbar, "Finished calculations. Write table.")

# write table to file
if writetable == true && tablefilename != ""

	# open file
	tablefile = open(tablefilename, "w")

	# write header
	write(tablefile, "Rotated Filter (dFilter = $dFilter, K_0 = $K_0) for fluid with viscosity = $mu, density = $rho\n" * "="^20)

	# strip first line "10x2 DataFrame"
	write(tablefile, replace(string(table),r"^[^\n]+\n","",1))

	# close file
	close(tablefile)
end

# write table to CSV-file. Can easily be reloaded
if writecsv == true && csvfilename != ""
	DataFrames.writetable(csvfilename, table)
end

# FIXME: Colors in Shell on Linux
# Finished
ProgressBar.next!(progressbar, "Finished successfully!")
