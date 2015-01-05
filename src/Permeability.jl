module Permeability

export K_xi, K_JJ

function K_xi(xi)
	return K_0 / (1.0 - xi)
end

function K_JJ(xi)
	return K_0 / (1.0 - xi) * (log(1.0 - xi) / (log(Phi_0_) + 0.931) + 1.0)
end

end