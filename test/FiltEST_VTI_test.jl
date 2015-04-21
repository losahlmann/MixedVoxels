using FiltEST_VTI
using SaveK
#import Zlib

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

write_file(housing, "test/housing.vti", true)
saveK(permeability.data, "test/savedK.sol")

vtireffile = open("test/housing_test_ref.vti")
vtitestfile = open("test/housing.vti")
solreffile = open("test/savedK_test_ref.sol")
soltestfile = open("test/savedK_test_ref.sol")
housing_ref = readall(vtireffile)
housing_test = readall(vtitestfile)
housing_ref = replace(housing_ref, r"<Date>(.+)<\/Date>", "<Date></Date>")
housing_test = replace(housing_test, r"<Date>(.+)<\/Date>", "<Date></Date>")
savedK_ref = readall(solreffile)
saveK_test = readall(soltestfile)

@test housing_ref == housing_test
@test savedK_ref == saveK_test

close(vtireffile)
close(vtitestfile)
close(solreffile)
close(soltestfile)

println("TEST parsing of XML tags")
@test readtagpair("""\t<Date>Di 21 Apr 11:35:57 2015</Date>""", "Date") == "Di 21 Apr 11:35:57 2015"
@test readtagpair("""\t\t<Date></Date>""", "Date") == ""
@test readtagpair("""\t\t<test></test>""", "Date") == false
@test readopentag("""\t <Material>""", "Material") == true
@test readopentag("""\t <Materials>""", "Material") == false
@test readopentag("""<Piece Extent="0 78 0 24 0 24">""", "Piece") == """Extent="0 78 0 24 0 24\""""
@test readopentag("""\t<DataArray type="UInt16" Name="Material" NumberOfComponents="1" format="appended" offset="0"/>""", "DataArray") == """type="UInt16" Name="Material" NumberOfComponents="1" format="appended" offset="0"/"""
@test readclosetag("""\t </Material>""", "Material") == true
@test readclosetag("""\t <Material>""", "Material") == false
@test readclosetag("""\t </Material>""", "Materials") == false
@test readclosetag("""\t </Date>""", "Material") == false
@test readentity("""type="UInt16" Name="Material" NumberOfComponents="1" format="appended" offset="0"/""", "Name") == "Material"
@test readentity("""type="UInt16" Name="Material" NumberOfComponents="1" format="appended" offset="0"/""", "type") == "UInt16"
@test readentity("""type="UInt16" Name="Material" NumberOfComponents="1" format="appended" offset="0"/""", "offset") == "0"
@test readentity("""type="UInt16" Name="Material" NumberOfComponents="1" format="appended" offset="0\"""", "offset") == "0"
@test readentity("""type="UInt16" Name="Material" NumberOfComponents="1" format="" offset="0"/""", "format") == ""
@test readentity("""type="UInt16" Name="Material" NumberOfComponents="1" format="appended" offset="0"/""", "Type") == false




println("TEST reading of FiltEST-VTI-File")

vtireffile = read_file("test/housing_test_ref.vti")
@test typeof(vtireffile) == FiltEST_VTIFile
@test vtireffile.traversedirection == "ZYX"
@test vtireffile.units == "mm"
@test vtireffile.dimension == 3
@test vtireffile.origin == [-0.8, -0.4, -0.4]
@test vtireffile.spacing == 0.4
@test vtireffile.dataorder == ["Material", "Permeability"]
@test haskey(vtireffile.voxeldata, "Material")
@test haskey(vtireffile.voxeldata, "Permeability")
mat = vtireffile.voxeldata["Material"]
perm = vtireffile.voxeldata["Permeability"]
@test typeof(mat) == DataArray
@test typeof(perm) == DataArray
#@test mat.offset == 0
#@test perm.offset == 40
@test mat.numberofblocks == 1
@test perm.numberofblocks == 1
@test mat.lastblocksize == 126
@test perm.lastblocksize == 504
@test mat.compressedblocksizes == [24]
@test perm.compressedblocksizes == [30]
@test size(mat.data) == (7,3,3)
@test size(perm.data) == (7,3,3)
@test mat.data[2,2,2] == Mt["Inflow"]
@test perm.data[4,2,2] == 7.0e-6


rm("test/housing.vti")
rm("test/savedK.sol")

println("TEST reading of Rotated-Filter-VTI-File")
housing_ref = read_file("test/housing_rotated_filter_ref.vti")
@test typeof(housing_ref) == FiltEST_VTIFile
@test housing_ref.traversedirection == "ZYX"
@test housing_ref.units == "mm"
@test housing_ref.dimension == 3
@test housing_ref.origin == [-2.4000000000000004, -0.8, -0.8]
@test housing_ref.spacing == 0.4
@test housing_ref.dataorder == ["Material", "Permeability"]
@test haskey(housing_ref.voxeldata, "Material")
@test haskey(housing_ref.voxeldata, "Permeability")
mat = housing_ref.voxeldata["Material"]
perm = housing_ref.voxeldata["Permeability"]
@test typeof(mat) == DataArray
@test typeof(perm) == DataArray
@test mat.numberofblocks == 3
@test perm.numberofblocks == 11
@test mat.lastblocksize == 24320
@test perm.lastblocksize == 31744
@test mat.compressedblocksizes == [322,380,235]
@test perm.compressedblocksizes == [215,1439,1477,1461,1435,1202,1469,1448,1493,1436,134]
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