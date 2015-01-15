# permeability of porous media in mm^-2
const K_0 = 7.0e-6

# solid volume fraction in porous media
const Phi_0_ = 0.1

# minimum volume fraction for mixed voxel
const eps = 1e-3

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

# read input parameters
const theta = [90, 80, 70, 60, 50]
const phi = [0, 10, 20]

#
const mode =["K_JJ", "K_xi", "fill", "remove"]

#
const plot = false

const table = true