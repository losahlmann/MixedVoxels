using FiltEST_VTI
using Base.Test

#@test hello("Julia") == "Hello, Julia"
#@test_approx_eq domath(2.0) 7.0

# TEST type mt
println(typeof(Mt))
a = Mt["Fluid"]
@test typeof(a) == Int64
@test typeof(a) == Material
# TODO: @test typeof(a) == UInt8

# TEST type DataArray
material = DataArray(Material, 1, fill(Mt["Fluid"], 3, 3, 7))
material.data[:,:,1:3] = Mt["Inflow"]
println(filter(x->x==Mt["Inflow"],material.data))
#material.data = 
println(material)
println(size(material.data))
println(length(material.data))

# TEST creation of FiltEST-VTI-File
housing = FiltEST_VTIFile()
#housing.set_extent(3, 3, 7)
housing.origin = [-0.4, -0.4, -0.4]
housing.spacing = 0.4
println(typeof(housing))
println(typeof(material))
add_data(housing, "Material", material)
println(typeof(housing.voxeldata))
write_file(housing, "housing.vti", false)
