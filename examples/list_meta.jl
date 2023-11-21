# List all metadata stored with minc file, similar to mincheader

using Minc2 # for reading MINC2 files
using ArgParse


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "in"
            help = "Input minc file"
            required = true
    end
    parse_args(ARGS, s)
end

args = parse_commandline()

a = Minc2.open_minc_file(args["in"])

for g in Minc2.groups(a)
    for i in Minc2.attributes(a,g)
        v = Minc2.read_attribute(a,g,i)
        @info "$(g):$(i) $(v)"
    end
end
Minc2.close_minc_file(a)