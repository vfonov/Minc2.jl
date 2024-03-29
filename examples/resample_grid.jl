using Minc2 # for reading MINC2 files
using Interpolations
using ArgParse
using StaticArrays

function resample_grid_volume(in_vol::Array{T,4}, 
            out_vol::Array{T,4}, 
            v2w::Minc2.AffineTransform{C}, 
            w2v::Minc2.AffineTransform{C}, 
            itfm::Vector{Minc2.AnyTransform{C,C}};
            interp::I=BSpline(Quadratic(Line(OnCell()))),
            fill=0.0,
            ftol=1.0/80,
            max_iter=10) where {C,T,I}

    # NEED to interpolate only over spatial dimensions
    in_vol_itp = extrapolate( interpolate( in_vol, (NoInterp(),interp,interp,interp)),fill)

    @simd for c in CartesianIndices(view(out_vol,1,:,:,:))
        orig = Minc2.transform_point(v2w, c )
        dst  = Minc2.transform_point(itfm, orig; ftol, max_iter )
        dst_v= Minc2.transform_point(w2v, dst ) .+ 1.0
        
        for i in eachindex(view(out_vol,:,1,1,1))
            @inbounds out_vol[i,c] = in_vol_itp(i, dst_v...)
        end
    end
    out_vol
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "in"
            help = "Input grid minc file"
            required = true
        "transform"
            help = "Transform"
            required = true
        "out"
            help = "Output grid minc file"
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
    end
    parse_args(ARGS, s)
end

args = parse_commandline()


in_vol,in_hdr,in_store_hdr = Minc2.read_minc_volume_std(args["in"], Float64)

if !isnothing(args["like"])
    out_vol,out_hdr,out_store_hdr = Minc2.empty_like_minc_volume_std(args["like"],Float64)
else
    out_vol = Array{Float64}(undef, size(in_vol)...)
    out_hdr = in_hdr
    out_store_hdr = in_store_hdr
end

v2w=Minc2.voxel_to_world(out_hdr)
w2v=Minc2.world_to_voxel(in_hdr)

tfm=Minc2.load_transforms(args["transform"])
itfm=Minc2.inv(tfm)

if args["order"] == 0     # nearest
    resample_grid_volume(in_vol, out_vol,v2w,w2v,itfm; interp=BSpline(Constant()),fill=args["fill"])
elseif args["order"] == 1 # linear
    resample_grid_volume(in_vol,out_vol,v2w,w2v,itfm; interp=BSpline(Linear()),fill=args["fill"])
elseif args["order"] == 2 # quadratic
    resample_grid_volume(in_vol,out_vol,v2w,w2v,itfm; interp=BSpline(Quadratic(Line(OnCell()))),fill=args["fill"])
elseif args["order"] == 3 # cubic
    resample_grid_volume(in_vol,out_vol,v2w,w2v,itfm; interp=BSpline(Cubic(Line(OnCell()))),fill=args["fill"])
end


Minc2.write_minc_volume_std(args["out"], UInt16, out_store_hdr, out_vol)
