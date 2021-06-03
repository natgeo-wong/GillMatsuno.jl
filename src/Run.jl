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

    for it = 1 : rt

        ϕfield!(ϕn,ϕ,u,v,G,D,S,Q)
        ufield!(un,u,ϕ,v,G,D,S)
        vfield!(vn,v,ϕ,u,G,D,S)
        ϕ .= ϕn
		u .= un
		v .= vn

        st += S.δt
        if st >= ot
            ot += S.ft
            if !isone(it)
            	mt += 1
				@info "$(now()) - Saving output at $(@sprintf("%06.2f",st)) model seconds ..."
                ϕf[:,:,mt] .= ϕ
                uf[:,:,mt] .= u
                vf[:,:,mt] .= v
            end
        end

    end

	return ϕf

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
