# FiltEST location
const filtest="/u/s/sahlmann/Documents/Fachpraktikum/FiltEST/filtest.sh"

# fluid properties
const mu = 1.81e-5 # viscosity
const rho = 1.2 # density

# permeability of porous media in mm^-2
const K_0 = 7.0e-6

# solid volume fraction in porous media
const Phi_0__ = [0.1]#[0.1, 0.2]

# minimum volume fraction for mixed voxel
const eps = 1e-3

# dimension of a voxel
const dVoxel_ = [0.4]#, 0.1]

# dimensions of housing
const dX = 26.4
const dY = 8.0
const dZ = 8.0
const dFilter = 1.6

# sizes in number of voxels
const NWall = 2
const NInlet = 4
const NOutlet = 4

# read input parameters
const theta_ = [80, 70, 60]#, 50]
const phi_ = [0, 10]#[0, 10, 20]

#
const method_ = {0.4 => ["K_JJ", "K_xi", "fill", "remove"],
				0.1 => ["fill", "remove"]
			}

# compress VTI-file
const zipVTI = true

#
const runFiltEST = false

#
const plot = false

#
const writetable = false
const writeCSV = false

# filenames, need to be consistent with flowSimInput.xml
# leave empty if no file should be written
const vtifilename = "housing.vti"
const solfilename = "savedK.sol"
const tablefilename = "results.txt"
const csvfilename = "results.csv"