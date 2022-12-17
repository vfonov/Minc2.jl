using Minc2 # for reading MINC2 files
using ArgParse


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "in"
            help = "Input minc file"
            required = true
        "--transform"
            help = "Transform"
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
        "--invert"
            help = "Invert transform"
            action = :store_true

    end
    parse_args(ARGS, s)
end

args = parse_commandline()


in_vol=Minc2.read_volume(args["in"],store=Float64)

if !isnothing(args["like"])
    out_vol = Minc2.empty_volume_like(args["like"], store=Float64)
else
    # out_vol = Array{Float64}(undef, size(in_vol)...)
    # out_hdr = in_hdr
    # out_store_hdr = in_store_hdr
    out_vol = Minc2.empty_volume_like(in_vol)
end

if !isnothing(args["transform"])
    tfm=Minc2.load_transforms(args["transform"])
    if args["invert"]
        tfm=Minc2.inv(tfm)
    end
else
    tfm=nothing
end

@info "Transform:" tfm

Minc2.resample_volume!(in_vol,out_vol;tfm,order=args["order"],
    fill=args["fill"],ftol=args["ftol"],max_iter=args["max_iter"])


#Minc2.write_minc_volume_std(args["out"], UInt16, out_store_hdr, out_vol)
Minc2.save_volume(args["out"],out_vol,store=UInt16, history=Minc2.format_history(ARGS))
