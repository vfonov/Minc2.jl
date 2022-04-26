using Minc2 # for reading MINC2 files
using Interpolations
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "in"
            help = "Input minc file"
            required = true
        "transform"
            help = "Transform"
            required = true
        "out"
            help = "Output minc file"
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

# xfm = Minc2.open_xfm_file(args["transform"])

# for i in 1:Minc2.get_n_concat(xfm)
#     @info "Transform:$(i) type:$(Minc2.get_n_type(xfm;n=i-1))"
#     if Minc2.get_n_type(xfm;n=i-1)==Minc2.MINC2_XFM_LINEAR
#         @info "matrix:",Minc2.get_linear_transform(xfm,n=i-1)
#     end
# end



in_vol,in_hdr,in_store_hdr = Minc2.read_minc_volume_std(args["in"], Float64)

@info "in_vol:",size(in_vol)
@info "in_hdr:",in_hdr

if args["order"] == 0     # nearest
    in_vol_itp = extrapolate( interpolate( in_vol, BSpline(Constant())),args["fill"])
elseif args["order"] == 1 # linear
    in_vol_itp = extrapolate( interpolate( in_vol, BSpline(Linear())),args["fill"])
elseif args["order"] == 2 # quadratic
    in_vol_itp = extrapolate( interpolate( in_vol, BSpline(Quadratic(Line(OnCell())))), args["fill"])
elseif args["order"] == 3 # cubic
    in_vol_itp = extrapolate( interpolate( in_vol, BSpline(Cubic(Line(OnCell())))), args["fill"])
end

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

@info "tfm:",tfm , "ixfm:",itfm

Threads.@threads for c in CartesianIndices(out_vol)
    orig = Minc2.transform_point(v2w, c )
    dst  = Minc2.transform_point(itfm, orig; ftol=args["ftol"], max_iter=args["max_iter"] )
    dst_v= Minc2.transform_point(w2v, dst ) .+ 1.0
    
    out_vol[c] = in_vol_itp( dst_v... )
    #out_vol[c] = sqrt(sum((orig - dst).^2))
end

Minc2.write_minc_volume_std(args["out"], UInt16, out_store_hdr, out_vol)
