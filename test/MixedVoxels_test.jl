using MixedVoxels

import FiltEST_VTI: Mt, DataArray

include("../src/Permeability.jl")

K_0 = 7.0e-6
Phi_0_ = 0.1

# TEST volumefraction()
println("TEST volumefraction()")

planes = [
	{ # Test 1: vol=0.5, 4 Schnittpunkte (sind Würfelecken),
	# bzw. Schnitt sind zwei diagonal gegenüberliegende Würfelkanten
	"n_p" => [-1/sqrt(2), 1/sqrt(2), 0],
	"d" => 0.0,
	"vol" => 0.5,
	"approx_eps" => 0},
	{ # Test 2: vol=0.25, 4 Schnittpunkte, davon zwei Würfelecken (eine Kante)
	"n_p" => [-1/sqrt(5),2/sqrt(5),0.0],
	"d" => 0.0,
	"vol" => 0.25,
	"approx_eps" => 0},
	{ # Test 3: vol=0.5, 4 Schnittpunkte (sind Würfelecken),
	# bzw. Schnitt sind zwei diagonal gegenüberliegende Würfelkanten
	"n_p" => [-1/sqrt(2), 0.0, 1/sqrt(2)],
	"d" => 0.0,
	"vol" => 0.5,
	"approx_eps" => 0},
	{ # Test 4: vol=0.5, 4 Schnittpunkte
	"n_p" => n_p = [-4.0/5.0, 0.0, 3.0/5.0],
	"d" => -0.1,
	"vol" => 0.5,
	"approx_eps" => 0},
	{ # Test 5: vol=0.125, 4 Schnittpunkte
	"n_p" => [-1/sqrt(2), 1/sqrt(2), 0.0],
	"d" => -1.0/(2.0*sqrt(2)),
	"vol" => 0.125,
	"approx_eps" => 0},
	{ # Test 6: vol=0.25, 4 Schnittpunkte
	"n_p" => [-1.0, 0.0, 0.0],
	"d" => -0.75,
	"vol" => 0.25,
	"approx_eps" => 0},
	{ # Test 7: vol=-1, kein Schnitt, Ebene liegt in Seitenfläche
	"n_p" => [-1.0, 0.0, 0.0],
	"d" => -1.0,
	"vol" => -1,
	"approx_eps" => 0},
	{ # Test 8: vol=0.9375, 4 Schnittpunkte
	"n_p" => [-1/sqrt(5), 2/sqrt(5), 0.0],
	"d" => 3.0/(2.0*sqrt(5)),
	"vol" => 15/16,
	"approx_eps" => 1e-10},
	{ # Test 9: vol=-1, kein Schnitt, nur eine Kante liegt in Ebene
	"n_p" => [-2/sqrt(5), 1/sqrt(5), 0.0],
	"d" => -2.0/sqrt(5),
	"vol" => -1,
	"approx_eps" => 0},
	{ # Test 10: vol=-1, kein Schnitt, nur eine Kante liegt in Ebene
	"n_p" => [-2/sqrt(5), 1/sqrt(5), 0.0],
	"d" => 1.0/sqrt(5),
	"vol" => -1,
	"approx_eps" => 0},
	{ # Test 11: vol=0.00260416, 3 Schnittpunkte
	"n_p" => [-1/sqrt(3), 1/sqrt(3), 1/sqrt(3)],
	"d" => -3.0/(4.0*sqrt(3)),
	"vol" => 1/384,
	"approx_eps" => 1e-10},
	{ # Test 12: vol=0.3177, 6 Schnittpunkte
	"n_p" => [-1/sqrt(3), 1/sqrt(3), 1/sqrt(3)],
	"d" => 1.0/(4.0*sqrt(3)),
	"vol" => 0.3177,
	"approx_eps" => 1e-4},
	{ # Test 13:  vol=0.347222, 5 Schnittpunkte,
	"n_p" => [-2/sqrt(17), 2/sqrt(17), 3/sqrt(17)],
	"d" => 1/sqrt(17),
	"vol" => 0.347222,
	"approx_eps" => 1e-6}
]

for plane in planes
	n_p = plane["n_p"]
	d = plane["d"]
	vol = plane["vol"]
	eps = plane["approx_eps"]
	
	# do loop as macro to be able to see which test failed
	@eval begin
		if $eps == 0
			@test volumefraction(MixedVoxels.intersection_points($n_p, $d), $n_p, $d) == $vol
		elseif $eps > 0
			@test_approx_eq_eps volumefraction(MixedVoxels.intersection_points($n_p, $d), $n_p, $d) $vol $eps
		end
	end
end

#println("TEST Plane")

println("TEST Geometry")

plane = Plane([1,0,0], 0.5)

geom = Geometry()
geom.extent = [1, 1, 1]
geom.h = 1.0

add!(geom, plane)

m = DataArray(Uint16, 1, Base.fill(Mt["Fluid"], 1, 1, 1))
p = DataArray(Float64, 1, Base.fill(1.0, 1, 1, 1))

discretise(geom, x->K_xi(x,Phi_0_), m, p)

@test m.data[1,1,1] == Mt["Porous"]
@test p.data[1,1,1] == 1.4e-5


# Test 1: vol=0.5, 4 Schnittpunkte (sind Würfelecken),
# bzw. Schnitt sind zwei diagonal gegenüberliegende Würfelkanten
#n_p = [-1/sqrt(2), 1/sqrt(2), 0]
#d = 0.0
#@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == 0.5
	
# Test 2: vol=0.25, 4 Schnittpunkte, davon zwei Würfelecken (eine Kante)
#n_p = [-1/sqrt(5),2/sqrt(5),0.0]
#d = 0.0
#@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == 0.25

# Test 3: vol=0.5, 4 Schnittpunkte (sind Würfelecken),
# bzw. Schnitt sind zwei diagonal gegenüberliegende Würfelkanten
#n_p = [-1/sqrt(2), 0.0, 1/sqrt(2)]
#d = 0.0
#@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == 0.5

# Test 4: vol=0.5, 4 Schnittpunkte
#n_p = [-4.0/5.0, 0.0, 3.0/5.0]
#d = -0.1
#@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == 0.5

# Test 5: vol=0.125, 4 Schnittpunkte
#n_p = [-1/sqrt(2), 1/sqrt(2), 0.0]
#d = -1.0/(2.0*sqrt(2))
#@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == 0.125

# Test 6: vol=0.25, 4 Schnittpunkte
#n_p = [-1.0, 0.0, 0.0]
#d = -0.75
#@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == 0.25

# Test 7: vol=-1, kein Schnitt, Ebene liegt in Seitenfläche
#n_p = [-1.0, 0.0, 0.0]
#d = -1.0
#@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == -1

# Test 8: vol=0.9375, 4 Schnittpunkte
#n_p = [-1/sqrt(5), 2/sqrt(5), 0.0]
#d = 3.0/(2.0*sqrt(5))
#@test_approx_eq volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) 15/16

# Test 9: vol=-1, kein Schnitt, nur eine Kante liegt in Ebene
#n_p = [-2/sqrt(5), 1/sqrt(5), 0.0]
#d = -2.0/sqrt(5)
#@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == -1

# Test 10: vol=-1, kein Schnitt, nur eine Kante liegt in Ebene
#n_p = [-2/sqrt(5), 1/sqrt(5), 0.0]
#d = 1.0/sqrt(5)
#@test volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) == -1

# Test 11: vol=0.00260416, 3 Schnittpunkte
#n_p = [-1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]	
#d = -3.0/(4.0*sqrt(3))
#@test_approx_eq volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) 1/384

# Test 12: vol=0.3177, 6 Schnittpunkte
#n_p = [-1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]	
#d = 1.0/(4.0*sqrt(3))
#@test_approx_eq_eps volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) 0.3177 1e-4

# Test 13:  vol=0.347222, 5 Schnittpunkte,
#n_p = [-2/sqrt(17), 2/sqrt(17), 3/sqrt(17)]
#d = 1/sqrt(17)
#@test_approx_eq_eps volumefraction(MixedVoxels.intersection_points(n_p, d), n_p, d) 0.347222 1e-6
