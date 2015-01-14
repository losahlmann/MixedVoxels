include("utils.jl")
include("Permeability.jl")
using FiltEST_VTI
using MixedVoxels

# global variables

# permeability of porous media in mm^-2
const K_0 = 7.0e-6

# solid volume fraction in porous media
const Phi_0_ = 0.1

# dimension of a voxel
const dVoxel = 0.4

# dimensions of housing
const dX = 26.4
const dY = 8.0
const dZ = 8.0
const dFilter = 1.6

# sizes in number of voxels
const NWall = 2
const NInlet = 4
const NOutlet = 4

# dimensions in voxels (without in-/outlet, wall)
# TODO: Rundung
const Nx = int64(dX/dVoxel)
const Ny = int64(dY/dVoxel)
const Nz = int64(dZ/dVoxel)

# extensions of housing (with without in-/outlet, wall)
const Nyt = Ny + 2 * NWall
const Nzt = Nz + 2 * NWall

# add layers for inlet and outlet in flow direction
const Nxt = Nx + 2 * NWall + NInlet + NOutlet

# read input parameters
theta = deg2rad(90)
phi = deg2rad(0)

# K_JJ, K_xi, fill, remove
# mode = K_xi

# convert angles to rad

# create FiltEST_VTIFile
housing = FiltEST_VTIFile()

# set origin in lower back left corner of fluid area (after Inlet)
housing.origin = [-dVoxel*(NWall+NInlet), -dVoxel*NWall, -dVoxel*NWall]

housing.spacing = dVoxel

# creta data arrays
material = DataArray(Material, 1, fill(Mt["Wall"], Nxt, Nyt, Nzt))

# TODO: Float64 vector, 6 Komponenten, aber in *.vti nur skalar?
permeability = DataArray(Float64, 6, fill(0.5e-6, Nxt, Nyt, Nzt))

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
for z in 1:Nz
	for y in 1:Ny
		for x in 1:Nx
			# origin voxel
			x0 = [x-1, y-1, z-1] * dVoxel

			# center voxel
			center = ([x-1, y-1, z-1] .+ 0.5) * dVoxel

			# transform plane into voxel coordinate system (einheitsvoxel)
			# distance between voxel origin and plane, divided by voxel length
			dist1 = (dot(x0, n_p) - d1) / dVoxel
			dist2 = (dot(x0, n_p) - d2) / dVoxel

			# intersect Voxel with first plane
			# "-" because plane is moved against normal direction
			intersectionpoints = intersection_points(n_p, -dist1)

			# calculate volume fraction
			if length(intersectionpoints) > 0
				# calculate volume fraction
				xi = 1 - volumefraction(intersectionpoints, n_p, -dist1)

				if (xi) # TODO
					# fluid voxel
				elseif (xi)
					# mixed porous voxel
				else
					# full porous voxel
				end
			end

			# intersection with second plane
			intersectionpoints = intersection_points(n_p, -dist2)

			# calculate volume fraction

			# voxel is entirely contained in filter medium: full porous voxel
			if dot(center, n_p)-d1 < 0 && dot(center, n_p)-d2 > 0
			end

			# else voxel is fluid
		end
	end
end

#
add_data(housing, "Material", material)
add_data(housing, "Permeability", permeability)


# write FiltEST_VTIFile
write_file(housing, "housing.vti", false)

# write SOL-File
solfile = open("savedK.sol", "w")

# MAYBE:
#for k in filter(x -> x < 1.0, permeability.data)
#	write(solfile, k)
#end

# write only permeability of porous voxels
write(solfile, filter(x -> x < 1.0, permeability.data))

# close file
close(solfile)
