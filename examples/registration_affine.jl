using Minc2 # for reading MINC2 files
using Interpolations
using ArgParse
using StaticArrays
using Optim
using ImageFiltering # for smoothing

function save_volume(fn, vol;store::Type{T}=Float32) where {T}
    Minc2.write_minc_volume_std(fn, store, 
        Minc2.create_header_from_v2w(size(vol.vol), vol.v2w,vector_dim=(length(size(vol.vol))==4)), vol.vol)
end
  
  
function read_volume(fn)
    in_vol,in_hdr,in_store_hdr = Minc2.read_minc_volume_std(fn, Float64)
    v2w=Minc2.voxel_to_world(in_hdr)

    return (vol=in_vol,v2w=v2w)
end
  

function resample_volume!(out,
            in;
            itfm=Minc2.IdentityTransform(),
            interp=BSpline(Linear()),
            fill=0.0 )

    out_vol=out.vol
    v2w=out.v2w

    in_vol=in.vol
    w2v=Minc2.inv(in.v2w)

    in_vol_itp = extrapolate( interpolate( in_vol, interp),fill)

    @simd for c in CartesianIndices(out_vol)
        orig = Minc2.transform_point(v2w, c )
        dst  = Minc2.transform_point(itfm, orig)
        dst_v= Minc2.transform_point(w2v, dst ) .+ 1.0
        
        @inbounds out_vol[c] = in_vol_itp( dst_v... )
    end
    out
end

function smooth_downsample(in; factor=2, smooth=factor*2)
    #@info size(in.vol),in.v2w
    # smooth
    vol_smooth=imfilter(in.vol, Kernel.gaussian( (smooth,smooth,smooth) ))
    out_vol = similar(in.vol, size(in.vol) .÷ factor)
    #mag = size(in.vol) ./ size(out_vol)
    out_v2w = Minc2.AffineTransform(
        SA_F64[ factor 0 0 ;0 factor 0 ;0 0 factor ] * in.v2w.rot,
        in.v2w.shift + (SA_F64[factor/2.0-0.5,factor/2.0-0.5,factor/2.0-0.5]'*in.v2w.rot)'
     )

    out=(vol=out_vol,v2w=out_v2w)
    # downsample
    resample_volume!(out,in)
    return out
end



function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "in1"
            help = "Input minc file 1"
            required = true
        "in2"
            help = "Input minc file 2"
            required = true
        "out"
            help = "Output transform file"
            required = true
        "--transform"
            help = "Initial transform"
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


#in_vol,in_hdr,in_store_hdr = Minc2.read_minc_volume_std(args["in"], Float64)
in1=read_volume(args["in1"])
in2=read_volume(args["in2"])

if !isnothing(args["transform"])
    @info "loading",args["transform"]
    itfm=Minc2.load_transforms(args["transform"])[1]
else
    itfm=Minc2.AffineTransform()
end


@info "itfm:",itfm 

function SSD(v1,v2)
    mapreduce((x,y)->(x-y)^2,(+),v1.vol,v2.vol)
end

function loss!(temp, v1, v2, param)
    itfm=Minc2.AffineTransform(reshape(param[1:9],(3,3)),param[10:12])
    resample_volume!(temp, v1, itfm=itfm)
    SSD(temp,v2)
end

#exit(1)
temp_=(vol=similar(in2.vol), v2w=in2.v2w)
resample_volume!(temp_, in1, itfm=itfm)
save_volume("init.mnc",temp_)

for s in [8,4,2] # scales
    @info "Level:", s
    # downsample and smooth input volume
    t_in1=smooth_downsample(in1, factor=s,smooth=s*2)
    t_in2=smooth_downsample(in2, factor=s,smooth=s*2)

    #save_volume("tmp_in1_$(s).mnc",t_in1)
    #save_volume("tmp_in2_$(s).mnc",t_in2)

    # 
    temp=(vol=similar(t_in2.vol), v2w=t_in2.v2w)
    # 
    init_par=vcat(reshape(itfm.rot,9),itfm.shift)
    result=optimize( x->loss!(temp,t_in1,t_in2,x), init_par,LBFGS(),
        Optim.Options(
                    store_trace = false,
                    show_trace = false)
                    )
    @info result
    param=Optim.minimizer(result)
    @info param
    global itfm = Minc2.AffineTransform(reshape(param[1:9],(3,3)),param[10:12])

    resample_volume!(temp_, in1, itfm=itfm)
    save_volume("res_$(s).mnc",temp_)
end

Minc2.save_transforms(args["out"], Minc2.AnyTransform{Float64,Float64}[tfm])
