#module Permeability

# TODO: method midpoint
#export K_xi, K_JJ, fill, remove, mixedvoxelmethod

function K_xi(xi, Phi_0_)
	# FIXME: xi=1 ?
	K = K_0 / (1.0 - xi)
	return Mt["Porous"], K
end

function K_JJ(xi, Phi_0_)
	K = K_0 / (1.0 - xi) * (log(1.0 - xi) / (log(Phi_0_) + 0.931) + 1.0)
	return Mt["Porous"], K
end

function fill(xi, Phi_0_)
	return Mt["Porous"], K_0
end

function remove(xi, Phi_0_)
	return Mt["Fluid"], 1.0
end

function midpoint(xi, Phi_0_)
	if xi <= 0.5
		return Mt["Porous"], K_0
	else
		return Mt["Fluid"], 1.0
	end
end

# TODO: with symbols
const mixedvoxelmethod = Dict{ASCIIString, Function}({
	"K_xi" => K_xi,
	"K_JJ" => K_JJ,
	"fill" => fill,
	"remove" => remove,
	"midpoint" => midpoint
})


#end