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
vtifilename = "housing.vti"
solfilename = "savedK.sol"
tablefilename = "results.txt"
csvfilename = "results.csv"

println("TEST creation of Rotated Filter VTI-File: 0.1 20 70 0.4 K_JJ")
include("../rotated_filter.jl")

vtireffile = open("test/housing_rotated_filter_ref.vti")
vticreffile = open("test/housing_rotated_filter_c_ref.vti")
vtitestfile = open("housing.vti")
solreffile = open("test/savedK_rotated_filter_ref.sol")
solcreffile = open("test/savedK_rotated_filter_c_ref.sol")
soltestfile = open("savedK.sol")

housing_ref = readall(vtireffile)
housing_cref = readall(vticreffile)
housing_test = readall(vtitestfile)
savedK_ref = readall(solreffile)
savedK_cref = readall(solcreffile)
savedK_test = readall(soltestfile)
housing_ref = replace(housing_ref, r"<Date>(.+)<\/Date>", "<Date></Date>")
housing_cref = replace(housing_cref, r"<Date>(.+)<\/Date>", "<Date></Date>")
housing_test = replace(housing_test, r"<Date>(.+)<\/Date>", "<Date></Date>")

@test housing_ref == housing_cref
@test savedK_ref == savedK_cref
@test housing_ref == housing_test
@test savedK_ref == savedK_test

close(vtireffile)
close(vticreffile)
close(vtitestfile)
close(solreffile)
close(solcreffile)
close(soltestfile)
#rm("housing.vti")
#rm("savedK.sol")