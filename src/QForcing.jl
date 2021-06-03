struct QStructs{FT<:Real}
	A  :: FT
	Lx :: FT
	Qx :: FT
	Ly :: FT
	Qy :: FT
end

function QfieldProperties(
	FT = Float64;
	A  :: Real = 1,
	Lx :: Real = 2,
	Qx :: Real = 0,
	Ly :: Real = 2,
	Qy :: Real = 0
)

	return QStructs{FT}(A,Lx,Qx,Ly,Qy)

end

function createQ(
	QParams::Array{QStructs{FT}},
	g::Grid{FT},
) where FT <: Real

	xc = g.xc
	yc = g.yc
	nQ = length(QParams)
	Q  = zeros(FT,g.nx,g.ny)

	for Qs in QParams, jj = 1 : g.ny, ii = 1 : g.nx
		if abs(xc[ii]) < Qs.Lx
			Q[ii,jj] += Qs.A * cos(pi*(xc[ii]-Qs.Qx)*0.5/Qs.Lx) *
							   exp(-((yc[jj]-Qs.Qy)/Qs.Ly)^2)
		end
	end

	return Q

end

function show(io::IO, Q::QStructs{FT}) where FT <: Real
    print(
		io,
		"The Q{$FT}-field forcing is defined as follows:\n",
		"           Amplitude : ", Q.A, '\n',
		"        Size (Lx,Ly) : ", (Q.Lx,Q.Ly), '\n',
		"    Position (Qx,Qy) : ", (Q.Qx,Q.Qy), '\n'
	)
end
