include("../rel_vel_error.jl")

println("TEST relative error in velocity field")

coarsematerial = DataArray(Material, 1, Base.fill(Mt["Fluid"], 7, 3, 3))
coarsevelocity = DataArray(Float32, 3, Base.fill(float32([1.0,0.0,0.0]), 7, 3, 3))
finematerial = DataArray(Material, 1, Base.fill(Mt["Fluid"], 28, 12, 12))
finevelocity = DataArray(Float32, 3, Base.fill(float32([1.0,0.0,0.0]), 28, 12, 12))

coarsehousing = FiltEST_VTIFile()
coarsehousing.origin = [-0.8, -0.4, -0.4]
coarsehousing.spacing = 0.4

finehousing = FiltEST_VTIFile()
finehousing.origin = [-0.8, -0.4, -0.4]
finehousing.spacing = 0.1

add_data(coarsehousing, "Material", coarsematerial)
add_data(coarsehousing, "Velocity_mm_per_sec", coarsevelocity)
add_data(finehousing, "Material", finematerial)
add_data(finehousing, "Velocity_mm_per_sec", finevelocity)

write_file(coarsehousing, "test/coarsehousing.vti", true)
write_file(finehousing, "test/finehousing.vti", true)

relvelocityerror!("test/coarsehousing.vti", "test/finehousing.vti")

vticoarsefile = read_file("test/coarsehousing.vti")

@test vticoarsefile.voxeldata["RelVelError"].data == Base.fill(0.0, 7, 3, 3)

rm("test/coarsehousing.vti")
rm("test/finehousing.vti")