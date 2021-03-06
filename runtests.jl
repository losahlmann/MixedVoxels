#using RunTests
using Base.Test
using Lint

include("src/utils.jl")


tests = ["test/Permeability_test.jl",
			"test/MixedVoxels_test.jl",
			"test/FiltEST_VTI_test.jl",
			"test/rotated_filter_test.jl",
			"test/rel_vel_error_test.jl"]

lints = ["rotated_filter.jl",
			"src/FiltEST_VTI.jl",
			"src/MixedVoxels.jl",
			"src/Permeability.jl",
			"src/SaveK.jl"]

# custom test handler
test_handler(r::Test.Success) = nothing
test_handler(r::Test.Failure) = print_with_color(:red, "Failure in $(r.expr)\n")
test_handler(r::Test.Error) = print_with_color(:red, "Error in $(r.expr)\n")



function runtest()
	for test in tests
		include(test)
	end

	# TEST tests
	#@test string("Hello ", "Julia") == "Hello Julia1" # rises Failure
	#@test_approx_eq 7.0000000000001 7.0 # does not work: returns no succes message
	#@test_throws ErrorException error("An error!")

	# TEST utils.jl
	println("TEST equalszero")
	@test equalszero(2e-11)
	@test !equalszero(2e-10)

	# TEST miscellaneous
	println("TEST miscellaneous")
	cubenormals = {[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]}
	@test cubenormals[1] == [1.0,0.0,0.0]

end

print_with_color(:green, "Running tests:\n")
Test.with_handler(runtest, test_handler)
print_with_color(:green, "Finished tests!\n")

#print_with_color(:green, "Running lint:\n")
#map(lintfile, lints)
#print_with_color(:green, "Finished lint!\n")
