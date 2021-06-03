struct Domain{FT<:Real}
    α :: FT
    β :: FT
    g :: FT
    H :: FT
end

function DomainProperties(
	FT = Float64;
	α :: Real = 0.2,
	β :: Real = 0.5,
	g :: Real = 1,
	H :: Real = 1
)
	return Domain{FT}(α,β,g,H)

end

function show(io::IO, D::Domain{FT}) where FT <: Real
    print(
		io,
		"The Domain{$FT} properties are as follows:\n",
		"     Rayleigh Damping (α) : ", D.α, '\n',
		"    Coriolis Gradient (β) : ", D.β, '\n',
		"              Gravity (g) : ", D.g, '\n',
		"       Surface Height (H) : ", D.H, '\n'
	)
end
