module SaveK

export saveK, readK

function saveK(permeability_data, filename)
	# create file
	solfile = open(filename, "w")

	# write only permeability of porous voxels
	porous = filter(perm -> perm < 1.0, permeability_data)

	# check if we selected the right amount of voxels
	#if length(porous) != length(filter(mat -> mat == Mt["Porous"], material.data))
	#	error("Porous voxel number in SOL-file is not correct!")
	#end

	for x in porous
		for i in 1:6
			write(solfile, x)
		end
	end

	# close file
	close(solfile)
end

function readK(filename, N)
	# open file
	solfile = open(filename, "r")

	permeability_data = Array(Float64, N)

	for n in 1:N
		permeability_data[n] = read(solfile, Float64)
		for i in 1:5
			read(solfile, Float64)
		end
	end

	# close file
	close(solfile)

	return permeability_data
end

end