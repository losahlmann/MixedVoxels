# TODO: Parameters in flowSimInput.xml

# FiltEST location
const filtest="/u/s/sahlmann/Documents/Fachpraktikum/FiltEST/filtest.sh"

# fluid properties
# viscosity in Pa.s
const 𝜇 = 1.81e-5
# density in kg/m^3
const 𝜌 = 1.2

# permeability of porous media in mm^2
const K_0 = 7.0e-6

# solid volume fraction in porous media
const 𝚽_0_ = [0.1, 0.2]

# minimum volume fraction for mixed voxel
const eps = 1e-3

# dimension of a voxel
const dVoxel = [0.4, 0.1]

# dimensions of housing
const L_x = 26.4
const L_y = 8.0
const L_z = 8.0
const dFilter = 1.6

# sizes in number of voxels
const NWall = 2
const NInlet = 4
const NOutlet = 4

# read input parameters
const 𝜉 = 0.0:0.1:1.0

#
const methods = {
	0.4 => ["K_JJ", "K_KC", "K_xi", "midpoint"],
	0.1 => ["K_JJ", "K_KC", "K_xi", "midpoint"]
}

# compress VTI-file?
const zipVTI = true

# do simulation or only write VTI-file?
const runFiltEST = true

# plot simulation results?
const plot = false

# write simulation results into files?
const writetable = true
const writeCSV = true

# filenames, need to be consistent with flowSimInput.xml
# leave empty if no file should be written
const vtifilename = "results/housing.vti"
const solfilename = "results/savedK.sol"
const tablefilename = "results/results.txt"
const csvfilename = "results/results.csv"