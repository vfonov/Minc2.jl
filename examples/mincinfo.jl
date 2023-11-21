# Show information about minc files, similar to mincinfo
# And store the information in a CSV file

using Minc2 # for reading MINC2 files
using ArgParse
using DataFrames
using CSV

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "in"
            help = "Input minc file"
            required = true
            nargs = '*'
        "-v","--verbose"
            help = "Verbose screen output"
            action = :store_true
        "--csv"
            help = "Output CSV file"
    end
    parse_args(ARGS, s)
end

function get_minc_info(minc_file)
    local file_info=Dict()
    a = Minc2.open_minc_file(minc_file)
    axis_to_names=Dict(Minc2.DIM_X=>"x",Minc2.DIM_Y=>"y",Minc2.DIM_Z=>"z",Minc2.DIM_VEC=>"v",Minc2.DIM_TIME=>"t")
    hdr = Minc2.store_header(a)
    file_info["file"]=minc_file

    for i = 1:length(hdr.start)
        file_info["$(axis_to_names[hdr.axis[i]])step"] = hdr.step[i]
        file_info["$(axis_to_names[hdr.axis[i]])"] = hdr.dims[i]
    end

    file_info["series_description"]= Minc2.get_attribute(a,"acquisition","series_description")
    file_info["mr_acq_type"]       = Minc2.get_attribute(a,"acquisition","mr_acq_type")
    file_info["flip_angle"]        = Minc2.get_attribute(a,"acquisition","flip_angle")
    file_info["echo_time"]         = Minc2.get_attribute(a,"acquisition","echo_time")
    file_info["slice_thickness"]   = Minc2.get_attribute(a,"acquisition","slice_thickness")
    file_info["repetition_time"]   = Minc2.get_attribute(a,"acquisition","repetition_time")
    file_info["acquisition_id"]    = Minc2.get_attribute(a,"acquisition","acquisition_id")

    if !isnothing(file_info["acquisition_id"])
        file_info["acquisition_id"]=parse(Int,file_info["acquisition_id"])
    else
        file_info["acquisition_id"]=0
    end

    Minc2.close_minc_file(a)
    return file_info
end

args = parse_commandline()
#@info args["in"], length(args["in"])

q=DataFrame(file=String[],
    x=Int[],y=Int[],z=Int[],t=Int[],
    xstep=Float64[],ystep=Float64[],zstep=Float64[],tstep=Float64[],
    series_description=String[],
    mr_acq_type=String[],
    flip_angle=Float64[],
    echo_time=Float64[],
    slice_thickness=Float64[],
    repetition_time=Float64[],
    acquisition_id=Int64[]
    )

for i in args["in"]
    try
      push!(q, get_minc_info(i);cols=:union)
    catch e
      @error "Error loading $(i)"
      @error e
    end
end

if !isnothing(args["csv"])
    CSV.write(args["csv"], q)
else
    if !args["verbose"]
        sort!(q, [:acquisition_id])

        transform!(q, AsTable([:x,:y,:z,:t])=>ByRow(q->join(skipmissing(q),"x")) => :size)
        transform!(q, AsTable([:xstep,:ystep,:zstep,:tstep])=>ByRow(q->join(map(x->round(x,digits=2),skipmissing(q)),"x")) => :step)

        show(IOContext(stdout, :compact => false), q[:,["file","series_description","size","step"]], allrows=true, summary=false, truncate=65535,show_row_number=false)
    else
        show(IOContext(stdout, :compact => false), q[:,["file","series_description","x","y","z","t"]],allrows=true,truncate=65535,summary=false,show_row_number=false)
    end
end
