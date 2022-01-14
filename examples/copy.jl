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
mri, hdr, store_hdr=Minc2.read_minc_volume(a,Float32)
history = Minc2.read_minc_history(a)
Minc2.close_minc_file(a)

println("History:", history)
println("store_hdr:",store_hdr)
println("hdr:",hdr)

b = Minc2.define_minc_file(store_hdr, UInt16, Float32)
Minc2.create_minc_file(b,args["out"])
Minc2.write_minc_history(b,"from Julia")
Minc2.write_minc_volume(b,mri)
Minc2.close_minc_file(b)
