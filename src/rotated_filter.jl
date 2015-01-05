using FiltEST_VTI

# global variables

#
K_0 = 7.0e-6

#
Phi_0_ = 0.1

#
dVoxel = 0.4

dX = 26.4
dY = 8.0
dZ = 8.0
dFilter = 1.6

# sizes in number of voxels
NWall = 2
NInlet = 4
NOutlet = 4

# dimensions in voxels (without in-/outlet, wall)
# TODO: Rundung
Nx::Int = dX/dVoxel
Ny::Int = dY/dVoxel
Nz::Int = dZ/dVoxel

# extensions of housing (with without in-/outlet, wall)
Nxt = Nx + 2 * NWall
Nyt = Ny + 2 * NWall
Nzt = Nz + 2 * NWall

# add layers for inlet and outlet in flow direction
Nxt += NInlet + NOutlet

# read input parameters
theta
phi
mode

# convert angles to rad

# create FiltEST_VTIFile
housing = FiltEST_VTIFile()

# set origin in lower back left corner of fluid area (after Inlet)
housing.origin = [-dVoxel*(NWall+NInlet), -dVoxel*NWall, -dVoxel*NWall]

housing.spacing = dVoxel

# creta data arrays
material = DataArray(Material, 1, fill(Mt["Wall"], Nxt, Nyt, Nzt))
# TODO: Float64 vector
permeability = DataArray(Float64, 6, fill(1.0, Nxt, Nyt, Nzt))

# add fluid area

# add inlet

# add outlet

# filter is limited by two planes

# center/pivot filter medium

# normal vector of planes (||.||=1)

# position vectors of planes

# determine material and permeability for each interior voxel

#
add_data(housing, "Material", material)
add_data(housing, "Permeability", permeability)


# write FiltEST_VTIFile
write_file(housing, "housing.vti", false)

# if a voxel is porous, write its permeability into the SOLFile
# TODO: use filter

# write SOL-File