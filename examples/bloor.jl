# Apply gaussian blurring to a minc file
# Equivalent to mincblur -fwhm <n> in.mnc out.mnc

using Minc2 
using ArgParse
using ImageFiltering

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "in"
            help = "Input minc file"
            required = true
        "--fwhm"
            default=1.0
            arg_type=Float64
            help = "Blooring kernel"
        "out"
            help = "Output minc file"
            required = true
    end
    parse_args(ARGS, s)
end

args = parse_commandline()


in_vol = Minc2.read_volume(args["in"],store=Float64)
out_vol = Minc2.empty_volume_like(in_vol)
kernel = KernelFactors.gaussian((args["fwhm"]/2.355, args["fwhm"]/2.355, args["fwhm"]/2.355))

imfilter!(Minc2.array(out_vol), Minc2.array(in_vol), kernel)


Minc2.save_volume(args["out"], out_vol, store=UInt16, history=Minc2.format_history(ARGS))
