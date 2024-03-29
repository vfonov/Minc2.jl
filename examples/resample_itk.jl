using Minc2 # for reading MINC2 files
using ArgParse


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "in"
            help = "Input minc file"
            required = true
        "transform"
            help = "Transform"
            required = false
            nargs = '*'
        "--out"
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

in_vol=Minc2.read_nifti_volume(args["in"], store=Float64)

if !isnothing(args["like"])
    out_vol = Minc2.empty_volume_like(args["like"],Float64)
else
    out_vol = Minc2.empty_volume_like(in_vol)
end

itfm = Minc2.read_ants_transform.(args["transform"])

Minc2.resample_volume!(out_vol, in_vol; itfm, order=args["order"], 
    fill=args["fill"], ftol=args["ftol"], max_iter=args["max_iter"])

Minc2.save_nifti_volume(args["out"],out_vol, store=Float32, history=Minc2.format_history(ARGS))
