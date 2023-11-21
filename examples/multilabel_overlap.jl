using Minc2 # for reading MINC2 files
using ArgParse
using StatsBase
using SplitApplyCombine

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "in"
            help = "Input minc files"
            required = true
            nargs = '*'
        "-v","--verbose"
            help = "Verbose screen output"
            action = :store_true
        "-b","--bg"
            help = "Include background (zero label)"
            action = :store_true
        "--maj"
            help = "Output majority file"
        "--ovl"
            help = "Output overlap file"
    end
    parse_args(ARGS, s)
end

args = parse_commandline()
if length(args["in"])<1
    @error "need input files"
    exit(1)
end

ref=Minc2.read_volume(args["in"][1], store=UInt8)
lbs=sort(unique(vec(ref.vol)))
nlab=length(lbs)
remap_fw=Dict(lbs[i]=>UInt8(i-1) for i in 1:length(lbs))
remap_bw=Dict(UInt8(i-1)=>lbs[i] for i in 1:length(lbs))

seg=combinedims([Minc2.read_volume(i, store=UInt8).vol for i in args["in"] ])

cseg=map(x->remap_fw[x],seg)


function gen_ovl(arr,nlab_)
    # discrete version of the formula from http://dx.doi.org/10.1109/TMI.2006.880587
    N=length(arr)
    cts=[count(==(element),arr) for element in 0:nlab_ ]
    U=sum(j->N*(N-1)/2-(N-j)*(N-j-1)/2, cts)
    I=sum(j->j*(j-1)/2, cts)
    L=argmax(cts)-1

    return I,U,L
end

function gen_ovl_nz(arr,nlab_)
    # discrete version of the formula from http://dx.doi.org/10.1109/TMI.2006.880587
    N=length(arr)
    cts=[count(==(element),arr) for element in 1:nlab_ ]
    U=sum(j->N*(N-1)/2-(N-j)*(N-j-1)/2, cts)
    I=sum(j->j*(j-1)/2, cts)

    L=U>0 ? argmax(cts) : 0 # deal with voxels where evything is bg 

    return I,U,L
end


if args["bg"]
    I,U,maj=invert(gen_ovl.(splitdimsview(cseg,(1,2,3)),nlab))
else
    I,U,maj=invert(gen_ovl_nz.(splitdimsview(cseg,(1,2,3)),nlab))
end

if !isnothing(args["maj"])
    omaj=map(x->remap_bw[x],maj)
    Minc2.save_volume(args["maj"],Minc2.Volume3D(omaj,Minc2.voxel_to_world(ref)),store=UInt8)
end


if !isnothing(args["ovl"])
    ovl=I./U
    ovl[isinf.(ovl) .|| isnan.(ovl)] .= 0.0 # happens when we have background only and no bg is included

    Minc2.save_volume(args["ovl"],Minc2.Volume3D(ovl,Minc2.voxel_to_world(ref)),store=Float32)
end

println( sum(I)/sum(U))
