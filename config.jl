# FiltEST location
const filtest="/u/s/sahlmann/Documents/Fachpraktikum/FiltEST/filtest.sh"

# fluid properties
const 𝜇 = 1.81e-5 # viscosity
const 𝜌 = 1.2 # density

# permeability of porous media in mm^-2
const K_0 = 7.0e-6

# solid volume fraction in porous media
const 𝚽_0_ = [0.1]#[0.1, 0.2]

# minimum volume fraction for mixed voxel
const eps = 1e-3

# dimension of a voxel
const dVoxel = [0.4]#, 0.1]

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
const 𝜃 = [80, 70, 60]#, 50]
const 𝜑 = [0, 10]#[0, 10, 20]

#
const methods = {
	0.4 => ["K_JJ", "K_xi", "fill", "remove"],
	0.1 => ["fill", "remove"]
}

# compress VTI-file?
const zipVTI = true

# do simulation or only write VTI-file?
const runFiltEST = false

# plot simulation results?
const plot = false

# write simulation results into files?
const writetable = false
const writeCSV = false

# filenames, need to be consistent with flowSimInput.xml
# leave empty if no file should be written
# FIXME: doppelung mit writetable/writeCSV
const vtifilename = "housing.vti"
const solfilename = "savedK.sol"
const tablefilename = "results.txt"
const csvfilename = "results.csv"