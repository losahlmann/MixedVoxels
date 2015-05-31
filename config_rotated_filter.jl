# TODO: Parameters in flowSimInput.xml

# FiltEST location
const filtest = "/u/s/sahlmann/Documents/Fachpraktikum/FiltEST/filtest.sh"

# fluid properties
# viscosity in Pa.s
const 𝜇 = 4.79e-2
# density in kg/m^3
const 𝜌 = 866.2

# permeability of porous media in mm^2
const K_0 = 1e-4

# solid volume fraction in porous media
const 𝚽_0_ = [0.1, 0.2]

# minimum volume fraction for mixed voxel
const 𝜀 = 1e-3

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
const 𝜃 = [90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40]
const 𝜑 = [0, 10, 20, 30, 40, 45]

#
const methods = {
	0.4 => ["K_JJ", "K_KC", "K_xi", "fill", "remove", "midpoint"],
	0.1 => ["fill", "remove", "midpoint"]
}

# compress VTI-file?
const zipVTI = true

# do simulation or only write VTI-file?
const runFiltEST = true

# write simulation results into files?
const writetable = true
const writeCSV = true

# filenames, need to be consistent with flowSimInput.xml
# leave empty if no file should be written
const vtifilename = "results/housing.vti"
const solfilename = "results/savedK.sol"
const tablefilename = "results/results.txt"
const csvfilename = "results/results.csv"