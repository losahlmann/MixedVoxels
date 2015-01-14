using Base.Test
using FiltEST_VTI
using MixedVoxels

include("../src/utils.jl")



#@test hello("Julia") == "Hello, Julia"
#@test_approx_eq domath(2.0) 7.0

# TEST utils.jl
@test equalszero(2e-11)
@test !equalszero(2e-10)

# TEST type mt
a = Mt["Fluid"]
@test typeof(a) == Uint16
@test typeof(a) == Material
# TODO: @test typeof(a) == UInt8

# TEST type DataArray
#material = DataArray(Material, 1, fill(Mt["Solid"], 7, 3, 3))
#material.data[2:6,2,2] = Mt["Fluid"]
#material.data[2,2,2] = Mt["Inflow"]
#material.data[6,2,2] = Mt["Outflow"]
#material.data[4,2,2] = Mt["Porous"]

#permeability = DataArray(Float64, 6, fill(1.0, 7, 3, 3))
#permeability.data[4,2,2] = 7.0e-6


#println(filter(x->x==Mt["Inflow"],material.data))
#material.data = 
#println(material)
#println(size(material.data))
#println(length(material.data))

# TEST creation of FiltEST-VTI-File
#housing = FiltEST_VTIFile()
#housing.origin = [-0.8, -0.4, -0.4]
#housing.spacing = 0.4

# FIXME: Reihenfolge in Data
#add_data(housing, "Permeability", permeability)
#material = DataArray(Material, 1, [Mt["Inflow"], Mt["Porous"], Mt["Outflow"]])

#add_data(housing, "Material", material)

#write_file(housing, "housing.vti", false)

#solfile = open("savedK.sol", "w")
#write(solfile, filter(x -> x < 1.0, permeability.data))
#close(solfile)


cubenormals = {[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]}
@test cubenormals[1] == [1.0,0.0,0.0]


# Test 1: vol=0.5, 4 Schnittpunkte (sind Würfelecken),
# bzw. Schnitt sind zwei diagonal gegenüberliegende Würfelkanten
n_p = [-1/sqrt(2), 1/sqrt(2), 0]
d = 0.0
#s = sort_vertices(intersectionpoints, n_p)
#println(s)
#println(MixedVoxels.polygon_area(s))
#@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == 0.5
	
# Test 2: vol=0.25, 4 Schnittpunkte, davon zwei Würfelecken (eine Kante)
n_p = [-1/sqrt(5),2/sqrt(5),0.0]
d = 0.0
@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == 0.25

# Test 3: vol=0.5, 4 Schnittpunkte (sind Würfelecken),
# bzw. Schnitt sind zwei diagonal gegenüberliegende Würfelkanten
n_p = [-1/sqrt(2), 0.0, 1/sqrt(2)]
d = 0.0
@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == 0.5
#println(MixedVoxels.intersection_points(n_p, d))
#println(volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d))

# Test 4: vol=0.5, 4 Schnittpunkte
n_p = [-4.0/5.0, 0.0, 3.0/5.0]
d = -0.1
@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == 0.5
# println(volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d))

# Test 5: vol=0.125, 4 Schnittpunkte
n_p = [-1/sqrt(2), 1/sqrt(2), 0.0]
d = -1.0/(2.0*sqrt(2))
@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == 0.125
# println(volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d))

# Test 6: vol=0.25, 4 Schnittpunkte
n_p = [-1.0, 0.0, 0.0]
d = -0.75
@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == 0.25
# intersectionpoints = MixedVoxels.intersection_points(n_p, d)
# println(volumefraction(intersectionpoints, n_p, d))
# #println(intersectionpoints)
# #s = sort_vertices(intersectionpoints, n_p)
# #println(s)
# #@test s == {[0.75,0.0,0.0],[0.75,1.0,1.0],[0.75,0.0,1.0],[0.75,1.0,0.0]}
# #println(MixedVoxels.polygon_area(s))

# Test 7: vol=-1, kein Schnitt, Ebene liegt in Seitenfläche
n_p = [-1.0, 0.0, 0.0]
d = -1.0
@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == -1

# Test 8: vol=0.9375, 4 Schnittpunkte
n_p = [-1/sqrt(5), 2/sqrt(5), 0.0]
d = 3.0/(2.0*sqrt(5))
@test_approx_eq volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) 15/16
#println(volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d))

# Test 9: vol=-1, kein Schnitt, nur eine Kante liegt in Ebene
n_p = [-2/sqrt(5), 1/sqrt(5), 0.0]
d = -2.0/sqrt(5)
@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == -1
#println(volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d))

# Test 10: vol=-1, kein Schnitt, nur eine Kante liegt in Ebene
n_p = [-2/sqrt(5), 1/sqrt(5), 0.0]
d = 1.0/sqrt(5)
@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == -1
# println(volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d))

# Test 11: vol=0.00260416, 3 Schnittpunkte
n_p = [-1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]	
d = -3.0/(4.0*sqrt(3))
#println(volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d))
@test_approx_eq volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) 1/384

# Test 12: vol=0.3177, 6 Schnittpunkte
n_p = [-1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]	
d = 1.0/(4.0*sqrt(3))
@test_approx_eq_eps volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) 0.3177 1e-4
# println(volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d))

# Test 13:  vol=0.347222, 5 Schnittpunkte,
n_p = [-2/sqrt(17), 2/sqrt(17), 3/sqrt(17)]
d = 1/sqrt(17)
@test_approx_eq_eps volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) 0.347222 1e-6
# println(volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d))

# # cube presentation
# n_p = [-1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]
# dist = 1/(2*sqrt(3))

# 	# std::cout << "--------" << std::endl
# 	# stPoints intersectionPoints = intersectionPointsPlaneVoxel(n_p, dist]
# 	# std::cout << "N=" << intersectionPoints.N << std::endl
# 	# std::cout << "volumeFraction  = " << volumeFraction(intersectionPoints, n_p, dist)
# 	# 			<< std::endl
# 	# std::cout << "volumeFraction2 = " << volumeFraction2(n_p, dist) << std::endl