using FiltEST_VTI
using Base.Test
import Zlib

# TEST type mt
# FiltEST expects UInt16 for material data (at least as data type in XML)
println("TEST type Mt")
a = Mt["Fluid"]
@test typeof(a) == Uint16
@test typeof(a) == Material
# TODO: @test typeof(a) == UInt8

# TEST type DataArray
println("TEST type DataArray")
material = DataArray(Material, 1, Base.fill(Mt["Solid"], 7, 3, 3))
material.data[2:6,2,2] = Mt["Fluid"]
material.data[2,2,2] = Mt["Inflow"]
material.data[6,2,2] = Mt["Outflow"]
material.data[4,2,2] = Mt["Porous"]

permeability = DataArray(Float64, 1, Base.fill(1.0, 7, 3, 3))
permeability.data[4,2,2] = 7.0e-6

# TEST creation of FiltEST-VTI-File
println("TEST creation of FiltEST-VTI-File")
housing = FiltEST_VTIFile()
housing.origin = [-0.8, -0.4, -0.4]
housing.spacing = 0.4

#material = DataArray(Material, 1, [Mt["Inflow"], Mt["Porous"], Mt["Outflow"], Mt["Fluid"]])

add_data(housing, "Material", material)
add_data(housing, "Permeability", permeability)

write_file(housing, "housing_test.vti", true)

#num2hex(Mt["Porous"])
#println(reinterpret(Uint8,[Mt["Inflow"], Mt["Porous"], Mt["Outflow"], Mt["Fluid"]]))

# uncompressed
#uncompressed = b"\x3c\x00\x0a\x00\x50\x00\x05\x00"
#println(uncompressed)
#println(uint16(uncompressed))
# compressed julia
#compressed_julia = b"\x78\x9c\xb3\xe1\x0a\x60\x05\x00\x01\xb7\x00\x9c"
#println("decompressed julia $(Zlib.decompress(compressed_julia))")
# == Uint8[60,10,80,5]

# compressed c++
#compressed_c = b"\x78\x9c\xb3\x61\xe0\x62\x08\x60\x60\x65\x00\x00\x03\x6e\x00\x9c"
#println("decompressed c $(Zlib.decompress(compressed_c))")
# == Uint8[60,0,10,0,80,0,5,0]