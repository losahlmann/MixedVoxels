#using RunTests
using Base.Test

include("src/utils.jl")


tests = ["test/MixedVoxels_test.jl",
			"test/FiltEST_VTI_test.jl"]

# custom test handler
test_handler(r::Test.Success) = println("Success on $(r.expr)")
test_handler(r::Test.Failure) = println("Failure in $(r.expr)")
test_handler(r::Test.Error) = println("Error in $(r.expr)")

println("Running tests:")
Test.with_handler(test_handler) do
	for test in tests
		include(test)
	end


	@test string("Hello ", "Julia") == "Hello Julia"
	@test_approx_eq 7.000000000001 7.0
	@test_throws ErrorException error("An error!")

	# TEST utils.jl
	@test equalszero(2e-11)
	@test !equalszero(2e-10)


	cubenormals = {[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]}
	@test cubenormals[1] == [1.0,0.0,0.0]

	#println(filter(x->x==Mt["Inflow"],material.data))
	#material.data = 
	#println(material)
	#println(size(material.data))
	#println(length(material.data))



	#solfile = open("savedK.sol", "w")
	#write(solfile, filter(x -> x < 1.0, permeability.data))
	#close(solfile)



end
println("Finished tests!")

# # cube presentation
# n_p = [-1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]
# dist = 1/(2*sqrt(3))
#exit(run_tests(["MixedVoxels_test.jl"]))
