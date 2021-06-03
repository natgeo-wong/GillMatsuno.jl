struct Simulation{FT<:Real}
     δt :: FT
     tt :: FT
	 ft :: FT
    fnc :: AbstractString
end

function CreateSimulation(
	FT = Float64;
	δt  :: Real,
	tt  :: Real,
	ft  :: Real,
	fnc :: AbstractString
)

	return Simulation{FT}(δt,tt,ft,fnc)

end

function show(io::IO, S::Simulation{FT}) where FT <: Real
    print(
		io,
		"Gill-Matsuno Simulation{$FT} Setup:\n",
		"           Time Step (δt) : ", S.δt, '\n',
		"           Stop Time (tt) : ", S.tt, '\n',
		"    Output Frequency (ft) : ", S.ft, '\n',
		"    NetCDF Filepath (fnc) : ", S.fnc, '\n'
	)
end
