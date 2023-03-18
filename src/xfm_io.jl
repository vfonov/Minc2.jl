"""
Low level: transformation types stored in .xfm file
"""
@enum XFM begin
    MINC2_XFM_LINEAR                 = Cint(minc2_simple.MINC2_XFM_LINEAR)
    MINC2_XFM_THIN_PLATE_SPLINE      = Cint(minc2_simple.MINC2_XFM_THIN_PLATE_SPLINE)
    MINC2_XFM_USER_TRANSFORM         = Cint(minc2_simple.MINC2_XFM_USER_TRANSFORM)
    MINC2_XFM_CONCATENATED_TRANSFORM = Cint(minc2_simple.MINC2_XFM_CONCATENATED_TRANSFORM)
    MINC2_XFM_GRID_TRANSFORM         = Cint(minc2_simple.MINC2_XFM_GRID_TRANSFORM)
end


"""
Low level: minc2_simple XFM transform handle
"""
mutable struct TransformHandle
    x::Ref
    function TransformHandle()
    ret = new( Ref(minc2_simple.minc2_xfm_allocate0()) )

    finalizer(x -> c"minc2_simple.minc2_xfm_destroy"(x[]), ret.x)
    return ret
    end
end


"""
Low level: Open transform xfm file, return handle
"""
function open_xfm_file(fname::String)::TransformHandle
    h = TransformHandle()
    @minc2_check minc2_simple.minc2_xfm_open(h.x[], fname)
    return h
end


"""
Low level: Save information into file from an open handle
"""
function save_xfm_file(h::TransformHandle, path::String)
    @minc2_check  minc2_simple.minc2_xfm_save(h.x[],path)
end

"""
Low level: Transform point x,y,z
"""
function transform_point(h::TransformHandle, xyz::Vector{Float64})::Vector{Float64}
    xyz_out=zeros(Float64,3)
    @minc2_check  minc2_simple.minc2_xfm_transform_point(h.x[],xyz,xyz_out)
    return xyz_out
end

"""
Low level: Inverse transform point x,y,z 
"""
function inverse_transform_point(h::TransformHandle, xyz::Vector{Float64})::Vector{Float64}
    xyz_out=zeros(Float64,3)
    @minc2_check  minc2_simple.minc2_xfm_inverse_transform_point(h.x[],xyz,xyz_out)
    return xyz_out
end

"""
Low level: Invert transform
"""
function invert_transform(h::TransformHandle)
    @minc2_check  minc2_simple.minc2_xfm_invert(h.x[])
end

"""
Low level: Get number of transformations
"""
function get_n_concat(h::TransformHandle)::Int64
    n = Ref{Int}(0)
    @minc2_check minc2_simple.minc2_xfm_get_n_concat( h.x[], n )
    return n[]
end

"""
Low level: Get transform type 
"""
function get_n_type(h::TransformHandle;n::Int64=0)::XFM
    t = Ref{Int}(0)
    @minc2_check minc2_simple.minc2_xfm_get_n_type( h.x[], n, t )
    return XFM(t[])
end

"""
Low level: extract reference to a grid file from open handle
"""
function get_grid_transform(h::TransformHandle;n::Int64=0)
    c_file=Ref{c"char *"}()
    inv=Ref{Int}(0)

    @minc2_check minc2_simple.minc2_xfm_get_grid_transform(h.x[],n,inv,c_file)
    r=unsafe_string(c_file[])

    Libc.free(c_file[])

    return (r,inv[]!=0)
end

"""
Low level: extract AffineTransform{Float64} from open handle
"""
function get_linear_transform(h::TransformHandle;n::Int64=0)::AffineTransform{Float64}
    mat=zeros(Float64,4,4)
    @minc2_check minc2_simple.minc2_xfm_get_linear_transform(h.x[], n, Base.unsafe_convert(Ptr{Cdouble},mat))
    return AffineTransform(mat')
end

"""
Low level: extacto transformation parameters from affine transform
"""
function get_linear_transform_param(h::TransformHandle;n::Int64=0,center::Union{Nothing,Vector{Float64}}=nothing)
    if isnothing(center)
        center=zeros(Float64,3)
    end

    translations=zeros(Float64,3)
    scales=zeros(Float64,3)
    scales=zeros(Float64,3)
    shears=zeros(Float64,3)
    rotations=zeros(Float64,3)
    @minc2_check minc2_simple.minc2_xfm_extract_linear_param(h.x[],n,center,translations,scales,shears,rotations)

    return (center=center,translations=translations,scales=scales,shears=shears,rotations=rotations)
end

"""
Low level: Append affine transform to an open transformation handle
"""
function append_linear_transform(h::TransformHandle, lin::AffineTransform)
    mat = Matrix{Float64}( Float64[lin.rot lin.shift;0 0 0 1]')
    @minc2_check minc2_simple.minc2_xfm_append_linear_transform(h.x[],
        Base.unsafe_convert(Ptr{Cdouble}, mat) )
end

"""
Low level: Append grid transform  to an open transformation handle
"""
function append_grid_transform(h::TransformHandle, grid_file::String; inv::Bool=false)
    @minc2_check minc2_simple.minc2_xfm_append_grid_transform(h.x[], grid_file, inv)
end

"""
Low level: concatenate two transfomations
"""
function concat_xfm(h::TransformHandle, i::TransformHandle)
    @minc2_check minc2_simple.minc2_xfm_concat_xfm(h.x[], i.x[])
end

"""
Low level: Load all transforms from open .XFM handle
"""
function load_transforms(h::TransformHandle)::Vector{AnyTransform}
    r=Vector{AnyTransform}()
    for i in 1:get_n_concat(h)
        t=get_n_type(h,n=i-1)
        if t==MINC2_XFM_LINEAR
            push!(r,get_linear_transform(h,n=i-1))
        elseif t==MINC2_XFM_GRID_TRANSFORM
            grid_fname, inv_grid=get_grid_transform(h,n=i-1)

            grid_vol, grid_hdr, grid_store_hdr = Minc2.read_minc_volume_std(grid_fname, Float64)

            if inv_grid
                push!(r,InverseGridTransform(voxel_to_world(grid_hdr), grid_vol))
            else
                push!(r,       GridTransform(voxel_to_world(grid_hdr), grid_vol))
            end
        else
            # unsupported type 
            throw( Minc2Error("Unsupported transform type: $t")) 
        end
    end
    return r
end


"""
Load transformations from .xfm file
"""
function load_transforms(fname::String)::Vector{AnyTransform} 
    h = Minc2.open_xfm_file(fname)
    load_transforms(h)
end

"""
Save transformations into .xfm file
"""
function save_transforms(fname::String, 
        xfm::Union{Vector{XFM}, XFM};
        grid_store::Type{T}=Float32 ) where {T, XFM<:AnyTransform}
    h = TransformHandle()
    grid_ctr=0# count grid files

    if typeof(xfm) <: Vector
        _xfm = xfm
    else
        _xfm = AnyTransform[xfm]
    end

    for x in xfm
        if x isa AffineTransform
            append_linear_transform(h,x)
        elseif x isa GridTransform
            grid_file_name=replace(fname,r"\.xfm$"=>"") * "_grid_$(grid_ctr).mnc"

            write_minc_volume_std(grid_file_name, grid_store, 
                create_header_from_v2w(size(x.vector_field), x.voxel_to_world,vector_dim=true), 
                    x.vector_field)
            
            append_grid_transform(h, grid_file_name)

            grid_ctr+=1
        elseif x isa InverseGridTransform
            grid_file_name=replace(fname,r"\.xfm$"=>"") * "_grid_$(grid_ctr).mnc"

            write_minc_volume_std(grid_file_name, grid_store, 
                create_header_from_v2w(size(x.vector_field), x.voxel_to_world,vector_dim=true), 
                    x.vector_field)
            
            append_grid_transform(h,grid_file_name; inv=true)
        elseif x isa IdentityTransform
            # do nothing, just skip iver
        else
            throw( Minc2Error("Unsupported transform: $x"))
        end 
    end
    save_xfm_file(h,fname)
end 
