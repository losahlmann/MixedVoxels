module SaveK

export saveK, readK

function saveK(permeability_data, filename)

	# create file
	solfile = open(filename, "w")

	# write only permeability of porous voxels
	porous = filter(perm -> perm < 1.0, permeability_data)

	# write each component of tensor (here: repeat scalar 6 times)
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

	# save permeability data for N voxels
	permeability_data = Array(Float64, N)

	for n in 1:N
		# read first scalar
		permeability_data[n] = read(solfile, Float64)
		# skip further 5 repetitions
		for i in 1:5
			read(solfile, Float64)
		end
	end

	# close file
	close(solfile)

	return permeability_data
end

end