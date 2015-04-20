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
eps = 1e-3

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

println("TEST creation of Rotated Filter VTI-File: 0.1 20 70 0.4 K_JJ")
include("../rotated_filter.jl")

vtireffile = open("test/housing_rotated_filter_ref.vti")
vticreffile = open("test/housing_rotated_filter_c_ref.vti")
vtitestfile = open(vtifilename)

housing_ref = readall(vtireffile)
housing_cref = readall(vticreffile)
housing_test = readall(vtitestfile)

housing_ref = replace(housing_ref, r"<Date>(.+)<\/Date>", "<Date></Date>")
housing_cref = replace(housing_cref, r"<Date>(.+)<\/Date>", "<Date></Date>")
housing_test = replace(housing_test, r"<Date>(.+)<\/Date>", "<Date></Date>")

@test housing_ref == housing_cref
@test housing_ref == housing_test

close(vtireffile)
close(vticreffile)
close(vtitestfile)

housing_ref = read_file("test/housing_rotated_filter_ref.vti")
housing_test = read_file(vtifilename)
@test housing_ref.voxeldata["Material"].data == housing_test.voxeldata["Material"].data
@test housing_ref.voxeldata["Permeability"].data == housing_test.voxeldata["Permeability"].data

savedK_ref = readK("test/savedK_rotated_filter_ref.sol", 2434)
savedK_cref = readK("test/savedK_rotated_filter_c_ref.sol", 2434)
savedK_test = readK(solfilename, 2434)

@test savedK_ref == savedK_test
@test savedK_ref == savedK_cref

rm(vtifilename)
rm(solfilename)