using FiltEST_VTI

function relvelocityerror!(coarsesol, finesol)

	coarsesolvti = read_file(coarsesol)
	finesolvti = read_file(finesol)

	const h = coarsesolvti.spacing
	const epsilon = finesolvti.spacing
	const a = int(h/epsilon)

	coarsevelocity = coarsesolvti.voxeldata["Velocity_mm_per_sec"].data
	finevelocity = finesolvti.voxeldata["Velocity_mm_per_sec"].data

	Nx, Ny, Nz = size(coarsevelocity)

	relvelerror = DataArray(Float64, 1, Base.fill(0.0, Nx, Ny, Nz))

	for z in 1:Nz, y in 1:Ny, x in 1:Nx
		uibar = [0.0,0.0,0.0]
		for zz in 1:a, yy in 1:a, xx in 1:a
			j = ([x,y,z].-1) * a + [xx,yy,zz]
			uibar += finevelocity[j...]
		end
		uibar /= a^3

		relvelerror.data[x,y,z] = norm(coarsevelocity[x,y,z] - uibar) / norm(uibar)
	end

	add_data(coarsesolvti, "RelVelError", relvelerror)
	write_file(coarsesolvti, coarsesol, true, false)

end

function relvelocityerrorpath(path)
	files = readdir(path)
	finesols = filter(r".*_0\.1", files)
	for finesol in finesols
		finesetting = match(r"(.*)_0\.1", finesol).captures[1]
		coarsesols = filter(Regex("$(finesetting)_0\.4"), files)
		for coarsesol in coarsesols
			println("do "*coarsesol)
			relvelocityerror!(path*coarsesol, path*finesol)
		end
	end
end