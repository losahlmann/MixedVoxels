# TODO: equalszero als macro
function equalszero(x, epss=1e-10)
	return (-epss <= x <= epss) ? true : false
end