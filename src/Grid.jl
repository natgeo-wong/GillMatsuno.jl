struct Grid{FT<:Real}
	xmin :: FT
	xmax :: FT
	nx   :: Int
	δx   :: FT
	xc   :: Array{FT}
	xf   :: Array{FT}
	ymin :: FT
	ymax :: FT
	ny   :: Int
	δy   :: FT
	yc   :: Array{FT}
	yf   :: Array{FT}
end

function GenerateGrid(
	FT = Float64;
	size, x, y,
)

	nx,ny = size
	xmin,xmax = x
	ymin,ymax = y

	xv = collect(range(xmin,xmax,length=nx+1)); δx = (xv[end] - xv[1]) / nx
	xf = xv[2:end]
	xc = (xv[1:(end-1)] .+ xv[2:end]) / 2
    yf = collect(range(ymin,ymax,length=ny+1)); δy = (yf[end] - yf[1]) / ny
	yc = (yf[1:(end-1)] .+ yf[2:end]) / 2

	return Grid{FT}(
		xmin,xmax,nx,δx,xc,xf,
		ymin,ymax,ny,δy,yc,yf
	)

end

function show(io::IO, g::Grid{FT}) where FT <: Real
    print(
		io,
		"The Grid{$FT} is defined as follows:\n",
		"                  domain : x ∈ [$(g.xmin), $(g.xmax)], y ∈ [$(g.ymin), $(g.ymin)]\n",
		"      resolution (nx,ny) : ", (g.nx,g.ny), '\n',
		"    grid spacing (δx,δy) : ", (g.δx,g.δy), '\n'
	)
end
