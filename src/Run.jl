function runGillMatsuno(
    S  :: Simulation{FT},
    G  :: Grid{FT},
    Qv :: Vector{QStructs{FT}},
    D  :: Domain{FT},
) where FT <: Real

	@info "$(now()) - Initializing the Gill-Matsuno Model ..."

    nx,ny = G.nx,G.ny

    Q  = createQ(Qv,G)
    rt = floor(Int32,S.tt/S.δt) + 1	# Number of running steps
    nt = floor(Int32,S.tt/S.ft) + 1	# Number of output steps
    st = 0.                 	# Simulation time
    ot = 0.                  	# Indicator as to whether to save output
    mt = 1                  	# Index in output matrix

	@info "$(now()) - The simulation will run for $(rt-1) timesteps, and output the fields for $(nt-1) timesteps"

	@info "$(now()) - Preallocating arrays for the wind and potential fields ..."

    ϕ  = zeros(FT,nx,ny);    u  = zeros(FT,nx,ny);    v  = zeros(FT,nx,ny+1);
    ϕn = zeros(FT,nx,ny);    un = zeros(FT,nx,ny);    vn = zeros(FT,nx,ny+1);
    ϕf = zeros(FT,nx,ny,nt); uf = zeros(FT,nx,ny,nt); vf = zeros(FT,nx,ny+1,nt);
	t  = zeros(nt)

    for it = 1 : rt

        ϕfield!(ϕn,ϕ,u,v,G,D,S,Q)
        ufield!(un,u,ϕ,v,G,D,S)
        vfield!(vn,v,ϕ,u,G,D,S)
        ϕ .= ϕn
		u .= un
		v .= vn

		if isnan(sum(ϕ)) || isnan(sum(u)) || isnan(sum(v))
			error("There are NaN values in the fields, indicating that the model is unstable.  Please reduce the timestep")
		end

        st += S.δt
        if st >= ot
            ot += S.ft
            if !isone(it)
            	mt += 1
				@info "$(now()) - Extracting output at $(@sprintf("%06.2f",st)) model seconds ..."
                ϕf[:,:,mt] .= ϕ
                uf[:,:,mt] .= u
                vf[:,:,mt] .= v
				t[mt] = st
            end
        end

    end

	savefields(uf,vf,ϕf,t,G,S)

end

function savefields(
	u :: Array{FT},
	v :: Array{FT},
	ϕ :: Array{FT},
	t :: Array{<:Real},
	G :: Grid{FT},
	S :: Simulation{FT}
) where FT <: Real

	if isfile(S.fnc)
        @info "$(now()) - Stale NetCDF file $(S.fnc) detected.  Overwriting ..."
        rm(S.fnc);
    end
	ds = NCDataset(S.fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created as output by the GillMatsuno.jl package on $(Dates.now())"
    ));

	ds.dim["xC"] = G.nx
    ds.dim["yC"] = G.ny
    ds.dim["xF"] = G.nx
    ds.dim["yF"] = G.ny + 1
    ds.dim["time"] = length(t)

	ncxC = defVar(ds,"xC",FT,("xC",),attrib = Dict("long_name" => "X center",))
	ncyC = defVar(ds,"yC",FT,("yC",),attrib = Dict("long_name" => "Y center",))
	ncxF = defVar(ds,"xF",FT,("xF",),attrib = Dict("long_name" => "X face",))
	ncyF = defVar(ds,"yF",FT,("yF",),attrib = Dict("long_name" => "Y face",))
	nct  = defVar(ds,"time",Float64,("time",),attrib = Dict("long_name" => "time",))

	ncu = defVar(ds,"u",FT,("xF","yC","time"),attrib = Dict(
		"long_name" => "u_component_of_wind",
		"full_name" => "Nondimensional u-Wind",
		"units"		=> "L/T",
	))

	ncv = defVar(ds,"v",FT,("xC","yF","time"),attrib = Dict(
		"long_name" => "v_component_of_wind",
		"full_name" => "Nondimensional v-Wind",
		"units"		=> "L/T",
	))

	ncϕ = defVar(ds,"ϕ",FT,("xC","yC","time"),attrib = Dict(
		"full_name" => "geopotential",
		"full_name" => "Nondimensional Height",
		"units"		=> "H",
	))

	ncxC[:] = G.xc
	ncyC[:] = G.yc
	ncxF[:] = G.xf
	ncyF[:] = G.yf

	nct[:] = t
	ncu[:] = u
	ncv[:] = v
	ncϕ[:] = ϕ

	close(ds)

	@info "$(now()) - Output fields from the Gill-Matsuno run has been saved."

end

function ϕfield!(
    ϕn::Array{FT}, ϕ::Array{FT}, u::Array{FT}, v::Array{FT},
    G::Grid{FT}, D::Domain{FT}, S::Simulation{FT}, Q::Array{FT}
) where FT <: Real

    for jj = 1 : G.ny, ii = 1 : G.nx

		vy = (v[ii,jj+1] - v[ii,jj]) / G.δy

        if ii != 1
			ux = (u[ii,jj] - u[ii-1,jj]) / G.δx
		else
			ux = (u[1,jj] - u[G.nx,jj]) / G.δx
        end

        ϕn[ii,jj] = ϕ[ii,jj] - S.δt * (D.α * ϕ[ii,jj] + D.g * D.H * (ux+vy) + Q[ii,jj])

    end

    return

end

function ufield!(
    un::Array{FT}, u::Array{FT}, ϕ::Array{FT}, v::Array{FT},
    G::Grid{FT}, D::Domain{FT}, S::Simulation{FT}
) where FT <: Real

    for jj = 1 : G.ny, ii = 1 : G.nx

        if ii != G.nx
			ϕx = (ϕ[ii+1,jj] - ϕ[ii,jj]) / G.δx
			vi = 0.25 * (v[ii,jj]+v[ii+1,jj]+v[ii,jj+1]+v[ii+1,jj+1])
		else
			ϕx = (ϕ[1,jj] - ϕ[G.nx,jj]) / G.δx
			vi = 0.25 * (v[G.nx,jj]+v[1,jj]+v[G.nx,jj+1]+v[1,jj+1])
        end

        un[ii,jj] = u[ii,jj] - S.δt * (D.α * u[ii,jj] - D.β * G.yc[jj] * vi + ϕx)

    end

    return

end

function vfield!(
    vn::Array{FT}, v::Array{FT}, ϕ::Array{FT}, u::Array{FT},
    G::Grid{FT}, D::Domain{FT}, S::Simulation{FT}
) where FT <: Real

    for jj = 2 : G.ny, ii = 1 : G.nx

		ϕy = (ϕ[ii,jj] - ϕ[ii,jj-1]) / G.δy

        if ii != 1
			ui = 0.25 * (u[ii,jj]+u[ii-1,jj]+u[ii,jj-1]+u[ii-1,jj-1])
		else
			ui = 0.25 * (u[1,jj]+u[G.nx,jj]+u[1,jj-1]+u[G.nx,jj-1])
        end

        vn[ii,jj] = v[ii,jj] - S.δt * (D.α * v[ii,jj] + D.β * G.yf[jj] * ui + ϕy)

    end

    return

end
