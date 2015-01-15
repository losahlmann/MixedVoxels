module Permeability

export K_xi, K_JJ, fill, remove

const modes = Dict{ASCIIString, Function}({
	"K_xi" = K_xi,
	"K_JJ" = K_JJ,
	"fill" = fill,
	"remove" = remove
})

function K_xi(xi)
	return Mt["Porous"], K_0 / (1.0 - xi)
end

function K_JJ(xi)
	return Mt["Porous"], K_0 / (1.0 - xi) * (log(1.0 - xi) / (log(Phi_0_) + 0.931) + 1.0)
end

function fill(xi)
	return Mt["Porous"], K_0
end

function remove(xi)
	return Mt["Fluid"], 1.0
end


end