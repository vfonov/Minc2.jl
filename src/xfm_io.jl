"""
Axis typed from MINC volume, TODO: make this compatible with NIFTI ?
"""
@enum XFM begin
    MINC2_XFM_LINEAR                 = Cint(minc2_simple.MINC2_XFM_LINEAR)
    MINC2_XFM_THIN_PLATE_SPLINE      = Cint(minc2_simple.MINC2_XFM_THIN_PLATE_SPLINE)
    MINC2_XFM_USER_TRANSFORM         = Cint(minc2_simple.MINC2_XFM_USER_TRANSFORM)
    MINC2_XFM_CONCATENATED_TRANSFORM = Cint(minc2_simple.MINC2_XFM_CONCATENATED_TRANSFORM)
    MINC2_XFM_GRID_TRANSFORM         = Cint(minc2_simple.MINC2_XFM_GRID_TRANSFORM)
end


"""
minc2_simple XFM transform handle
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
Open transform xfm file, return handle
"""
function open_xfm_file(fname::String)::TransformHandle
    h = TransformHandle()
    @minc2_check minc2_simple.minc2_xfm_open(h.x[], fname)
    return h
end


"""
Save information into file
"""
function save_xfm_file(h::TransformHandle, path::String)
    @minc2_check  minc2_simple.minc2_xfm_save(h.x[],path)
end

"""
Transform point x,y,z using libminc
"""
function transform_point(h::TransformHandle, xyz::Vector{Float64})::Vector{Float64}
    xyz_out=zeros(Float64,3)
    @minc2_check  minc2_simple.minc2_xfm_transform_point(h.x[],xyz,xyz_out)
    return xyz_out
end

"""
Inverse transform point x,y,z  using libminc
"""
function inverse_transform_point(h::TransformHandle, xyz::Vector{Float64})::Vector{Float64}
    xyz_out=zeros(Float64,3)
    @minc2_check  minc2_simple.minc2_xfm_inverse_transform_point(h.x[],xyz,xyz_out)
    return xyz_out
end

"""
Invert transform
"""
function invert_transform(h::TransformHandle)
    @minc2_check  minc2_simple.minc2_xfm_invert(h.x[])
end

"""
Get number of transformations
"""
function get_n_concat(h::TransformHandle)::Int64
    n = Ref{Int}(0)
    @minc2_check minc2_simple.minc2_xfm_get_n_concat( h.x[], n )
    return n[]
end

"""
Get transform type 
"""
function get_n_type(h::TransformHandle;n::Int64=0)::XFM
    t = Ref{Int}(0)
    @minc2_check minc2_simple.minc2_xfm_get_n_type( h.x[], n, t )
    return XFM(t[])
end

function get_grid_transform(h::TransformHandle;n::Int64=0)
    c_file=Ref{c"char *"}()
    inv=Ref{Int}(0)

    @minc2_check minc2_simple.minc2_xfm_get_grid_transform(h.x[],n,inv,c_file)
    r=unsafe_string(c_file[])

    Libc.free(c_file[])

    return (r,inv[]!=0)
end

function get_linear_transform(h::TransformHandle;n::Int64=0)::AffineTransform{Float64}
    mat=zeros(Float64,4,4)
    @minc2_check minc2_simple.minc2_xfm_get_linear_transform(h.x[], n, Base.unsafe_convert(Ptr{Cdouble},mat))
    return AffineTransform(mat')
end

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
Append affine transform
"""
function append_linear_transform(h::TransformHandle, lin::AffineTransform)
    
    mat = Matrix{Float64}( Float64[lin.rot lin.shift;0 0 0 1]')
    @minc2_check minc2_simple.minc2_xfm_append_linear_transform(h.x[],
        Base.unsafe_convert(Ptr{Cdouble}, mat) )
end

"""
Append grid transform 
"""
function append_grid_transform(h::TransformHandle, grid_file::String; inv::Bool=false)
    @minc2_check minc2_simple.minc2_xfm_append_grid_transform(h.x[], grid_file, inv)
end

"""

"""
function concat_xfm(h::TransformHandle, i::TransformHandle)
    @minc2_check minc2_simple.minc2_xfm_concat_xfm(h.x[], i.x[])
end

"""
High level interface to load complete transform into memory 
"""
function load_transforms(h::TransformHandle)::Vector{AnyTransform{Float64,Float64}}
    r=Vector{AnyTransform{Float64,Float64}}()
    for i in 1:get_n_concat(h)
        t=get_n_type(h,n=i-1)
        if t==MINC2_XFM_LINEAR
            push!(r,get_linear_transform(h,n=i-1))
        elseif t==MINC2_XFM_GRID_TRANSFORM
            grid_fname, inv_grid=get_grid_transform(h,n=i-1)

            grid_vol, grid_hdr, grid_store_hdr = Minc2.read_minc_volume_std(grid_fname, Float64)

            if inv_grid
                push!(r,InverseGridTransform{Float64,Float64}(voxel_to_world(grid_hdr),grid_vol))
            else
                push!(r,       GridTransform{Float64,Float64}(voxel_to_world(grid_hdr),grid_vol))
            end
        else
            # unsupported type 
            throw( Minc2Error("Unsupported transform type: $t")) 
        end
    end
    return r
end


"""
High level interface to load complete transform into memory 
"""
function load_transforms(fname::String)::Vector{AnyTransform{Float64,Float64}} 
    h = Minc2.open_xfm_file(fname)
    load_transforms(h)
end


function save_transforms(fname::String, xfm::Vector{AnyTransform{Float64,Float64}};
        grid_store::Type{T}=Float32 ) where {T}
    h = TransformHandle()
    grid_ctr=0# count grid files
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
        else
            throw( Minc2Error("Unsupported transform: $x"))
        end 
    end
    save_xfm_file(h,fname)
end 