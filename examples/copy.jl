# Create a new minc file with same dimensions
# And copy voxel data

using Minc2 # for reading MINC2 files
using ArgParse


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "in"
            help = "Input minc file"
            required = true
        "out"
            help = "Output file"
            required = true
    end
    parse_args(ARGS, s)
end

args = parse_commandline()

a = Minc2.open_minc_file(args["in"])
mri, store_hdr = Minc2.read_minc_volume_raw(a, Float32)
Minc2.close_minc_file(a)

b = Minc2.define_minc_file(store_hdr, Float32, Float64)
Minc2.create_minc_file(b,args["out"])
Minc2.write_minc_history(b,"from Julia")
Minc2.write_minc_volume_raw(b, convert(Array{Float64,4}, mri))
Minc2.close_minc_file(b)
