using Minc2 # for reading MINC2 files
using Interpolations
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "in"
            help = "Input minc file"
            required = true
        "transform"
            help = "Transform"
            required = true
        "out"
            help = "Output minc file"
            required = true
        "--like"
            help = "Reference "
    end
    parse_args(ARGS, s)
end

args = parse_commandline()

xfm = Minc2.open_xfm_file(args["transform"])

for i in 1:Minc2.get_n_concat(xfm)
    @info "Transform:$(i) type:$(Minc2.get_n_type(xfm;n=i-1))"
    if Minc2.get_n_type(xfm;n=i-1)==Minc2.MINC2_XFM_LINEAR
        @info "matrix:",Minc2.get_linear_transform(xfm,n=i-1)
    end
end

in_vol,in_hdr,in_store_hdr=Minc2.read_minc_volume_std(args["in"], Float64)

@info "in_vol:",size(in_vol)
@info "in_hdr:",in_hdr


in_vol_itp = extrapolate( interpolate( in_vol, BSpline(Linear())),0.0)

out_vol=zeros(size(in_vol))
# simple case 
v2w=Minc2.voxel_to_world(in_hdr)
w2v=Minc2.world_to_voxel(in_hdr)

tfm=Minc2.AffineTransform()
tfm.mat=Minc2.get_linear_transform(xfm,n=0)

itfm=Minc2.AffineTransform()
itfm.mat=inv(tfm.mat)

@info "v2w:",v2w
@info "w2v:",w2v
@info "ixfm:",itfm


for c in CartesianIndices(out_vol)

    dst = Minc2.transform_point(w2v, Minc2.transform_point(itfm, Minc2.transform_point(v2w, Float64[Tuple(c)...] .- 1.0 ))) .+ 1.0

    out_vol[c] = in_vol_itp( dst[1], dst[2], dst[3] )
end


Minc2.write_minc_volume_std(args["out"], UInt16, in_store_hdr, out_vol)

