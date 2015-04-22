import FiltEST_VTI: Mt

K_0 = 7.0e-6

include("../src/Permeability.jl")


# TEST permeability.jl
println("TEST Permeability")

# full voxel
for method in mixedvoxelmethod
	@eval @test $method[2](0, 0.1) == (Mt["Porous"], K_0)
end

# empty voxel
for method in mixedvoxelmethod
	@eval @test $method[2](1, 0.1) == (Mt["Fluid"], 1.0)
end

# half voxel
@test K_xi(0.5, 0.1) == (Mt["Porous"], 1.4e-5)
@test K_JJ(0.5, 0.1)[1] == Mt["Porous"]
@test_approx_eq_eps K_JJ(0.5, 0.1)[2] 2.10751e-5 1e-10
@test K_KC(0.5, 0.1)[1] == Mt["Porous"]
@test_approx_eq_eps K_KC(0.5, 0.1)[2] 3.29307e-5 1e-10
@test fill(0.5, 0.1) == (Mt["Porous"], K_0)
@test remove(0.5, 0.1) == (Mt["Fluid"], 1.0)
@test midpoint(0.5, 0.1) == (Mt["Porous"], K_0)

# 2/3 filled voxel
@test K_xi(1/3, 0.1) == (Mt["Porous"], 1.05e-5)
@test K_JJ(1/3, 0.1)[1] == Mt["Porous"]
@test_approx_eq_eps K_JJ(1/3, 0.1)[2] 1.36040e-5 1e-10
@test K_KC(1/3, 0.1)[1] == Mt["Porous"]
@test_approx_eq_eps K_KC(1/3, 0.1)[2] 1.75656e-5 1e-10
@test fill(1/3, 0.1) == (Mt["Porous"], K_0)
@test remove(1/3, 0.1) == (Mt["Fluid"], 1.0)
@test midpoint(1/3, 0.1) == (Mt["Porous"], K_0)

# 1/3 filled voxel
@test K_xi(2/3, 0.1) == (Mt["Porous"], 2.1e-5)
@test K_JJ(2/3, 0.1)[1] == Mt["Porous"]
@test_approx_eq_eps K_JJ(2/3, 0.1)[2] 3.78206e-5 1e-10
@test K_KC(2/3, 0.1)[1] == Mt["Porous"]
@test_approx_eq_eps K_KC(2/3, 0.1)[2] 7.80626e-5 1e-10
@test fill(2/3, 0.1) == (Mt["Porous"], K_0)
@test remove(2/3, 0.1) == (Mt["Fluid"], 1.0)
@test midpoint(2/3, 0.1) == (Mt["Fluid"], 1.0)
