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

groups = Minc2.groups(a)

for g in groups
    attrs = Minc2.attributes(a,g)
    for i in attrs
        v = Minc2.read_attribute(a,g,i)
        @info "$(g):$(i) $(v)"
    end
end
Minc2.close_minc_file(a)