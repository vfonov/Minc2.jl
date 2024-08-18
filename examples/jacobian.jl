using Minc2 # for reading MINC2 files
using ArgParse
using Interpolations


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "transform"
            help = "Transform"
            required = true
        "ref"
            help = "Reference minc file"
            required = true
        "out"
            help = "Output minc file"
            required = true
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

jac = Minc2.empty_volume_like(args["ref"], store=Float64)

tfm=Minc2.load_transforms(args["transform"])
if args["invert"]
    tfm=Minc2.inv(tfm)
end

if args["order"]==1
    interp=BSpline(Linear())
elseif args["order"]==2  # default
    interp=BSpline(Quadratic(Line(OnCell())))
elseif args["order"]==3
    interp=BSpline(Cubic(Line(OnCell())))
else
    @error "Unsupported order" args["order"]
end

Minc2.calculate_jacobian!(tfm, jac;ftol=args["ftol"],max_iter=args["max_iter"],interp=interp)
Minc2.save_volume(args["out"],jac, store=UInt16, history=Minc2.format_history(ARGS))
