config = true

# FiltEST location
filtest="/u/s/sahlmann/Documents/Fachpraktikum/FiltEST/filtest.sh"

# fluid properties
𝜇 = 1.81e-5 # viscosity
𝜌 = 1.2 # density

# permeability of porous media in mm^-2
K_0 = 7.0e-6

# solid volume fraction in porous media
𝚽_0_ = [0.1]

# minimum volume fraction for mixed voxel
𝜀 = 1e-3

# dimension of a voxel
dVoxel = [0.4]

# dimensions of housing
L_x = 26.4
L_y = 8.0
L_z = 8.0
dFilter = 1.6

# sizes in number of voxels
NWall = 2
NInlet = 4
NOutlet = 4

# read input parameters
𝜃 = [70]
𝜑 = [20]

#
methods = {
	0.4 => ["K_JJ"]
}

# compress VTI-file?
zipVTI = true

# do simulation or only write VTI-file?
runFiltEST = false

# plot simulation results?
plot = false

# write simulation results into files?
writetable = false
writeCSV = false

# filenames, need to be consistent with flowSimInput.xml
# leave empty if no file should be written
vtifilename = "test/housing.vti"
solfilename = "test/savedK.sol"
tablefilename = "results.txt"
csvfilename = "results.csv"

println("TEST creation of Rotated Filter VTI-File: 0.1 70 20 0.4 K_JJ")
include("../rotated_filter.jl")


housing_ref = read_file("test/housing_rotated_filter_ref.vti")
housing_c_ref = read_file("test/housing_rotated_filter_c_ref.vti")
housing_test = read_file(vtifilename)
@test housing_ref.voxeldata["Material"].data == housing_test.voxeldata["Material"].data
@test housing_c_ref.voxeldata["Material"].data == housing_test.voxeldata["Material"].data
@test_approx_eq_eps housing_ref.voxeldata["Permeability"].data housing_test.voxeldata["Permeability"].data 1e-10
@test_approx_eq_eps housing_c_ref.voxeldata["Permeability"].data housing_test.voxeldata["Permeability"].data 1e-10


savedK_ref = readK("test/savedK_rotated_filter_ref.sol", 2434)
savedK_cref = readK("test/savedK_rotated_filter_c_ref.sol", 2434)
savedK_test = readK(solfilename, 2434)

@test_approx_eq_eps savedK_ref savedK_cref 1e-12
@test_approx_eq_eps savedK_cref savedK_test 1e-12

rm(vtifilename)
rm(solfilename)