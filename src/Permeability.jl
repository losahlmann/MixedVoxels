function K_xi(xi, Phi_0_)
	if xi == 1
		return Mt["Fluid"], 1.0
	else
		K = K_0 / (1.0 - xi)
		return Mt["Porous"], K
	end
end

function K_JJ(xi, Phi_0_)
	if xi == 1
		return Mt["Fluid"], 1.0
	else
		K = K_0 / (1.0 - xi) * (log(1.0 - xi) / (log(Phi_0_) + 0.931) + 1.0)
		return Mt["Porous"], K
	end
end

function K_KC(xi, Phi_0_)
	if xi == 1
		return Mt["Fluid"], 1.0
	else
		K = K_0 / (1.0 - xi)^2 * (Phi_0_ / (1 - Phi_0_) * xi + 1.0)^3
		return Mt["Porous"], K
	end
end

function fill(xi, Phi_0_)
	if xi == 1
		return Mt["Fluid"], 1.0
	else
		return Mt["Porous"], K_0
	end
end

function remove(xi, Phi_0_)
	if xi == 0
		return Mt["Porous"], K_0
	else
		return Mt["Fluid"], 1.0
	end
end

function midpoint(xi, Phi_0_)
	if xi <= 0.5
		return Mt["Porous"], K_0
	else
		return Mt["Fluid"], 1.0
	end
end

# TODO: with symbols
mixedvoxelmethod = Dict{ASCIIString, Function}({
	"K_xi" => K_xi,
	"K_JJ" => K_JJ,
	"fill" => fill,
	"remove" => remove,
	"midpoint" => midpoint
})
