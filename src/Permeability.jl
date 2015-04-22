function K_xi(xi::Float64, Phi_0_::Float64)
	if xi == 1.0
		return Mt["Fluid"], 1.0
	else
		K = K_0 / (1.0 - xi)
		return Mt["Porous"], K
	end
end

function K_JJ(xi::Float64, Phi_0_::Float64)
	if xi == 1.0
		return Mt["Fluid"], 1.0
	else
		K = K_0 / (1.0 - xi) * (log(1.0 - xi) / (log(Phi_0_) + 0.931) + 1.0)
		return Mt["Porous"], K
	end
end

function K_KC(xi::Float64, Phi_0_::Float64)
	if xi == 1.0
		return Mt["Fluid"], 1.0
	else
		K = K_0 / (1.0 - xi)^2 * (Phi_0_ / (1.0 - Phi_0_) * xi + 1.0)^3
		return Mt["Porous"], K
	end
end

function fill(xi::Float64, Phi_0_::Float64)
	if xi == 1.0
		return Mt["Fluid"], 1.0
	else
		return Mt["Porous"], K_0
	end
end

function remove(xi::Float64, Phi_0_::Float64)
	if xi == 0.0
		return Mt["Porous"], K_0
	else
		return Mt["Fluid"], 1.0
	end
end

function midpoint(xi::Float64, Phi_0_::Float64)
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
	"K_KC" => K_KC,
	"fill" => fill,
	"remove" => remove,
	"midpoint" => midpoint
})
