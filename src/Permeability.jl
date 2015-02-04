#module Permeability

# TODO: method midpoint
#export K_xi, K_JJ, fill, remove, mixedvoxelmethod

function K_xi(xi, a...)
	K = K_0 / (1.0 - xi)
	return Mt["Porous"], K
end

function K_JJ(xi, Phi_0_, a...)
	K = K_0 / (1.0 - xi) * (log(1.0 - xi) / (log(Phi_0_) + 0.931) + 1.0)
	return Mt["Porous"], K
end

function fill(xi, a...)
	return Mt["Porous"], K_0
end

function remove(xi, a...)
	return Mt["Fluid"], 1.0
end

const mixedvoxelmethod = Dict{ASCIIString, Function}({
	"K_xi" => K_xi,
	"K_JJ" => K_JJ,
	"fill" => fill,
	"remove" => remove
})


#end