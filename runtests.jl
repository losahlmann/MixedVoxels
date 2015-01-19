#using RunTests
using Base.Test

include("src/utils.jl")


tests = ["test/MixedVoxels_test.jl",
			"test/FiltEST_VTI_test.jl"]

println("Running tests:")
for test in tests
  include(test)
end

#@test hello("Julia") == "Hello, Julia"
#@test_approx_eq domath(2.0) 7.0

# TEST utils.jl
@test equalszero(2e-11)
@test !equalszero(2e-10)



#println(filter(x->x==Mt["Inflow"],material.data))
#material.data = 
#println(material)
#println(size(material.data))
#println(length(material.data))



#solfile = open("savedK.sol", "w")
#write(solfile, filter(x -> x < 1.0, permeability.data))
#close(solfile)


cubenormals = {[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]}
@test cubenormals[1] == [1.0,0.0,0.0]


# # cube presentation
# n_p = [-1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]
# dist = 1/(2*sqrt(3))
println("Finished tests!")
#exit(run_tests(["MixedVoxels_test.jl"]))
