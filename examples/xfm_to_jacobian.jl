using Minc2 # for reading MINC2 files
using Interpolations
using ArgParse
using StaticArrays
using LinearAlgebra

function calculate_jacobian(
            out_vol::Array{T,3}, 
            v2w::Minc2.AffineTransform{C},
            tfm::Vector{Minc2.AnyTransform{C,C}};
            interp::I=BSpline(Quadratic(Line(OnCell()))),
            fill=0.0,
            ftol=1.0/80,
            max_iter=10) where {C,T,I}

    # calculate scaling matrix from the voxel to world matrix
    f = svd(v2w.rot)
    dir_cos = f.U * f.Vt
    sc = Base.inv(v2w.rot * Base.inv(dir_cos))
    @info "Scaling matrix" sc

    # First step: generate vector field of transformations
    vector_field = Array{T}(undef, 3, size(out_vol)...)

    @simd for c in CartesianIndices(out_vol)
        orig = Minc2.transform_point(v2w, c )
        dst  = Minc2.transform_point(tfm, orig; ftol, max_iter )
        
        @inbounds vector_field[:,c] .= dst # .- orig
    end

    # Second step: calculate jacobian determinant 
    vector_field_itp = extrapolate( interpolate( vector_field, (NoInterp(),interp,interp,interp)), Flat())

    @simd for c in CartesianIndices(out_vol)
        grad = hcat([ Interpolations.gradient( vector_field_itp, i, Tuple(c)...) for i in 1:3 ]...)
        out_vol[c] = det(grad'*sc)
    end

    out_vol
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "transform"
            help = "Input transform file"
            required = true
        "out"
            help = "Output jacobian minc file"
            required = true
        "--like"
            help = "Reference "
        "--ftol"
            help = "Tolerance for inverse nonlinear transform"
            default=1.0/80
            arg_type=Float64
        "--max_iter"
            help = "Maximum number of iterations per voxel for inverse transform"
            default=10
            arg_type=Int
        "--order"
            help = "Interpolation order [0-3]"
            arg_type = Int
            default = 2
            range_tester = x->0<=x<4
        "--fill"
            help = "Fill value"
            arg_type = Float64
            default = 0.0
        "--invert"
            help = "Invert transform"
            action = :store_true

    end
    parse_args(ARGS, s)
end

#####
begin 
    args = parse_commandline()

    tfm=Minc2.load_transforms(args["transform"])

    if args["invert"]
        tfm=Minc2.inv(tfm)
    end

    @info "Input transform" tfm



    if !isnothing(args["like"])
        # search for reference volume in .xfm file
        out_vol,out_hdr,out_store_hdr = Minc2.empty_like_minc_volume_std(args["like"],Float64)
    else
        idx_grid=findfirst(x->(x isa Minc2.GridTransform || x isa Minc2.InverseGridTransform), tfm)
        if isnothing(idx_grid)
            @error "Can't find grid to use as reference, provide --like "
            exit 
        end

        out_vol = Array{Float64}(undef, size(tfm[idx_grid].vector_field)[2:4]...)
        out_hdr = Minc2.create_header_from_v2w(size(out_vol), tfm[idx_grid].voxel_to_world,vector_dim=false)
        out_store_hdr = out_hdr
    end

    v2w=Minc2.voxel_to_world(out_hdr)

    if args["order"] == 0     # nearest
        #in_vol_itp = extrapolate( interpolate( in_vol, BSpline(Constant())),args["fill"])
        calculate_jacobian(out_vol, v2w,tfm; interp=BSpline(Constant()),fill=args["fill"])
    elseif args["order"] == 1 # linear
        calculate_jacobian(out_vol, v2w,tfm; interp=BSpline(Linear()),fill=args["fill"])
    elseif args["order"] == 2 # quadratic
        calculate_jacobian(out_vol, v2w,tfm; interp=BSpline(Quadratic(Line(OnCell()))),fill=args["fill"])
    elseif args["order"] == 3 # cubic
        calculate_jacobian(out_vol, v2w,tfm; interp=BSpline(Cubic(Line(OnCell()))),fill=args["fill"])
    end

    Minc2.write_minc_volume_std(args["out"], UInt16, out_store_hdr, out_vol)
end
