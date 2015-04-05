import FiltEST_VTI: Mt

K_0 = 7.0e-6

include("../src/Permeability.jl")


# TEST permeability.jl
println("TEST Permeability")

# full voxel
for method in mixedvoxelmethod
	@eval @test $method[2](0, 0.1) == (Mt["Porous"], K_0)
end