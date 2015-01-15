using FiltEST_VTI

# TEST type mt
println("TEST type Mt")
a = Mt["Fluid"]
@test typeof(a) == Uint16
@test typeof(a) == Material
# TODO: @test typeof(a) == UInt8

# TEST type DataArray
println("TEST type DataArray")
material = DataArray(Material, 1, fill(Mt["Solid"], 7, 3, 3))
material.data[2:6,2,2] = Mt["Fluid"]
material.data[2,2,2] = Mt["Inflow"]
material.data[6,2,2] = Mt["Outflow"]
material.data[4,2,2] = Mt["Porous"]

permeability = DataArray(Float64, 6, fill(1.0, 7, 3, 3))
permeability.data[4,2,2] = 7.0e-6

# TEST creation of FiltEST-VTI-File
println("TEST creation of FiltEST-VTI-File")
housing = FiltEST_VTIFile()
housing.origin = [-0.8, -0.4, -0.4]
housing.spacing = 0.4

#material = DataArray(Material, 1, [Mt["Inflow"], Mt["Porous"], Mt["Outflow"]])

# FIXME: Reihenfolge in Data
add_data(housing, "Material", material)
add_data(housing, "Permeability", permeability)

write_file(housing, "housing_test.vti", false)