using Minc2
using BenchmarkTools
using Profile


if false
    @info "Nifti io"
    b=@benchmark Minc2.read_nifti_volume("template0.nii.gz")
    show(IOContext(stdout), MIME("text/plain"), b)
    println()

    @info "MINC io"
    b=@benchmark Minc2.read_volume("std_sub-1053744_ses-2_t1w.mnc")
    show(IOContext(stdout), MIME("text/plain"), b)
    println()


    @info "XFM io"
    b=@benchmark Minc2.load_transforms("MRT05___nl_0_NL.xfm")
    show(IOContext(stdout), MIME("text/plain"), b)
    println()
end


v=Minc2.read_volume("std_sub-1053744_ses-2_t1w.mnc")

if false
@info "identity resample"
b=@benchmark Minc2.resample_volume($v;like=$v)
show(IOContext(stdout), MIME("text/plain"), b)
println()
end

#itfm=Minc2.load_transforms("cplx.xfm")
itfm=Minc2.load_transforms("MRT05___nl_0_NL.xfm")
@info "nl resample, naive"
b=@benchmark Minc2.resample_volume($v;like=$v,itfm=$itfm[1])
show(IOContext(stdout), MIME("text/plain"), b)
println()

# itfm=Minc2.load_transforms("cplx.xfm")
# @info "nl resample, naive"
# b=@benchmark Minc2.resample_volume_bcast($v;like=$v,itfm=$itfm)
# show(IOContext(stdout), MIME("text/plain"), b)
# println()
