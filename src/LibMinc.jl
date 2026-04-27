module LibMinc

using libminc_jll
export libminc_jll

using CEnum: CEnum, @cenum

@cenum minc2_dimensions::UInt32 begin
    MINC2_DIM_UNKNOWN = 0
    MINC2_DIM_X = 1
    MINC2_DIM_Y = 2
    MINC2_DIM_Z = 3
    MINC2_DIM_TIME = 4
    MINC2_DIM_VEC = 5
    MINC2_DIM_END = 255
end

@cenum minc2_type::Int32 begin
    MINC2_ORIGINAL = 0
    MINC2_BYTE = 1
    MINC2_SHORT = 3
    MINC2_INT = 4
    MINC2_FLOAT = 5
    MINC2_DOUBLE = 6
    MINC2_STRING = 7
    MINC2_UBYTE = 100
    MINC2_USHORT = 101
    MINC2_UINT = 102
    MINC2_SCOMPLEX = 1000
    MINC2_ICOMPLEX = 1001
    MINC2_FCOMPLEX = 1002
    MINC2_DCOMPLEX = 1003
    MINC2_MAX_TYPE_ID = 1004
    MINC2_UNKNOWN = -1
end

@cenum minc2_xfm::UInt32 begin
    MINC2_XFM_LINEAR = 1
    MINC2_XFM_THIN_PLATE_SPLINE = 2
    MINC2_XFM_USER_TRANSFORM = 3
    MINC2_XFM_CONCATENATED_TRANSFORM = 4
    MINC2_XFM_GRID_TRANSFORM = 5
    MINC2_XFM_END = 6
end

struct minc2_dimension
    id::Cint
    length::Cint
    irregular::Cint
    step::Cdouble
    start::Cdouble
    have_dir_cos::Cint
    dir_cos::NTuple{3, Cdouble}
end

mutable struct minc2_info_iterator end

const minc2_info_iterator_handle = Ptr{minc2_info_iterator}

mutable struct minc2_file_iterator end

const minc2_file_iterator_handle = Ptr{minc2_file_iterator}

@cenum var"##Ctag#277"::Int32 begin
    MINC2_SUCCESS = 0
    MINC2_ERROR = -1
end

mutable struct minc2_file end

const minc2_file_handle = Ptr{minc2_file}

mutable struct minc2_xfm_file end

const minc2_xfm_file_handle = Ptr{minc2_xfm_file}

struct minc2_tags
    n_volumes::Cint
    n_tag_points::Cint
    tags_volume1::Ptr{Cdouble}
    tags_volume2::Ptr{Cdouble}
    weights::Ptr{Cdouble}
    structure_ids::Ptr{Cint}
    patient_ids::Ptr{Cint}
    labels::Ptr{Ptr{Cchar}}
end

const minc2_tags_handle = Ptr{minc2_tags}

function minc2_allocate(h)
    ccall((:minc2_allocate, libminc2_simple), Cint, (Ptr{minc2_file_handle},), h)
end

function minc2_allocate0()
    ccall((:minc2_allocate0, libminc2_simple), minc2_file_handle, ())
end

function minc2_init(h)
    ccall((:minc2_init, libminc2_simple), Cint, (minc2_file_handle,), h)
end

function minc2_free(h)
    ccall((:minc2_free, libminc2_simple), Cint, (minc2_file_handle,), h)
end

function minc2_destroy(h)
    ccall((:minc2_destroy, libminc2_simple), Cint, (minc2_file_handle,), h)
end

function minc2_open(h, path)
    ccall((:minc2_open, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cchar}), h, path)
end

function minc2_open_rdwr(h, path)
    ccall((:minc2_open_rdwr, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cchar}), h, path)
end

function minc2_define(h, store_dims, store_data_type, data_type)
    ccall((:minc2_define, libminc2_simple), Cint, (minc2_file_handle, Ptr{minc2_dimension}, Cint, Cint), h, store_dims, store_data_type, data_type)
end

function minc2_create(h, path)
    ccall((:minc2_create, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cchar}), h, path)
end

function minc2_close(h)
    ccall((:minc2_close, libminc2_simple), Cint, (minc2_file_handle,), h)
end

function minc2_ndim(h, ndim)
    ccall((:minc2_ndim, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cint}), h, ndim)
end

function minc2_nelement(h, nelement)
    ccall((:minc2_nelement, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cint}), h, nelement)
end

function minc2_data_type(h, _type)
    ccall((:minc2_data_type, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cint}), h, _type)
end

function minc2_storage_data_type(h, _type)
    ccall((:minc2_storage_data_type, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cint}), h, _type)
end

function minc2_slice_ndim(h, slice_ndim)
    ccall((:minc2_slice_ndim, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cint}), h, slice_ndim)
end

function minc2_setup_standard_order(h)
    ccall((:minc2_setup_standard_order, libminc2_simple), Cint, (minc2_file_handle,), h)
end

function minc2_get_representation_dimensions(h, dims)
    ccall((:minc2_get_representation_dimensions, libminc2_simple), Cint, (minc2_file_handle, Ptr{Ptr{minc2_dimension}}), h, dims)
end

function minc2_get_store_dimensions(h, dims)
    ccall((:minc2_get_store_dimensions, libminc2_simple), Cint, (minc2_file_handle, Ptr{Ptr{minc2_dimension}}), h, dims)
end

function minc2_compare_voxel_dimensions(one, two)
    ccall((:minc2_compare_voxel_dimensions, libminc2_simple), Cint, (Ptr{minc2_dimension}, Ptr{minc2_dimension}), one, two)
end

function minc2_compare_dimensions(one, two)
    ccall((:minc2_compare_dimensions, libminc2_simple), Cint, (Ptr{minc2_dimension}, Ptr{minc2_dimension}), one, two)
end

function minc2_load_complete_volume(h, buffer, representation_type)
    ccall((:minc2_load_complete_volume, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cvoid}, Cint), h, buffer, representation_type)
end

function minc2_save_complete_volume(h, buffer, representation_type)
    ccall((:minc2_save_complete_volume, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cvoid}, Cint), h, buffer, representation_type)
end

function minc2_set_scaling(h, use_global_scaling, use_slice_scaling)
    ccall((:minc2_set_scaling, libminc2_simple), Cint, (minc2_file_handle, Cint, Cint), h, use_global_scaling, use_slice_scaling)
end

function minc2_set_volume_range(h, value_min, value_max)
    ccall((:minc2_set_volume_range, libminc2_simple), Cint, (minc2_file_handle, Cdouble, Cdouble), h, value_min, value_max)
end

function minc2_set_slice_range(h, start, value_min, value_max)
    ccall((:minc2_set_slice_range, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cint}, Cdouble, Cdouble), h, start, value_min, value_max)
end

function minc2_world_to_voxel(h, world, voxel)
    ccall((:minc2_world_to_voxel, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cdouble}, Ptr{Cdouble}), h, world, voxel)
end

function minc2_world_to_voxel_vec(h, n, stride, world, voxel)
    ccall((:minc2_world_to_voxel_vec, libminc2_simple), Cint, (minc2_file_handle, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), h, n, stride, world, voxel)
end

function minc2_voxel_to_world(h, voxel, world)
    ccall((:minc2_voxel_to_world, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cdouble}, Ptr{Cdouble}), h, voxel, world)
end

function minc2_voxel_to_world_vec(h, n, stride, voxel, world)
    ccall((:minc2_voxel_to_world_vec, libminc2_simple), Cint, (minc2_file_handle, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), h, n, stride, voxel, world)
end

function minc2_copy_metadata(src, dst)
    ccall((:minc2_copy_metadata, libminc2_simple), Cint, (minc2_file_handle, minc2_file_handle), src, dst)
end

function minc2_write_hyperslab(h, start, count, buffer, representation_type)
    ccall((:minc2_write_hyperslab, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid}, Cint), h, start, count, buffer, representation_type)
end

function minc2_read_hyperslab(h, start, count, buffer, representation_type)
    ccall((:minc2_read_hyperslab, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid}, Cint), h, start, count, buffer, representation_type)
end

function minc2_data_type_name(minc2_type_id)
    ccall((:minc2_data_type_name, libminc2_simple), Ptr{Cchar}, (Cint,), minc2_type_id)
end

function minc2_dim_type_name(minc2_dim_id)
    ccall((:minc2_dim_type_name, libminc2_simple), Ptr{Cchar}, (Cint,), minc2_dim_id)
end

function minc2_get_attribute_type(h, group, attr, minc2_type_)
    ccall((:minc2_get_attribute_type, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}), h, group, attr, minc2_type_)
end

function minc2_get_attribute_length(h, group, attr, attr_length)
    ccall((:minc2_get_attribute_length, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}), h, group, attr, attr_length)
end

function minc2_read_attribute(h, group, attr, buf, buf_size)
    ccall((:minc2_read_attribute, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cvoid}, Cint), h, group, attr, buf, buf_size)
end

function minc2_write_attribute(h, group, attr, buf, buf_size, minc2_type_)
    ccall((:minc2_write_attribute, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cvoid}, Cint, Cint), h, group, attr, buf, buf_size, minc2_type_)
end

function minc2_delete_attribute(h, group, attr)
    ccall((:minc2_delete_attribute, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cchar}, Ptr{Cchar}), h, group, attr)
end

function minc2_delete_group(h, group)
    ccall((:minc2_delete_group, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cchar}), h, group)
end

function minc2_allocate_info_iterator()
    ccall((:minc2_allocate_info_iterator, libminc2_simple), minc2_info_iterator_handle, ())
end

function minc2_free_info_iterator(it)
    ccall((:minc2_free_info_iterator, libminc2_simple), Cint, (minc2_info_iterator_handle,), it)
end

function minc2_stop_info_iterator(it)
    ccall((:minc2_stop_info_iterator, libminc2_simple), Cint, (minc2_info_iterator_handle,), it)
end

function minc2_start_group_iterator(h, group_it)
    ccall((:minc2_start_group_iterator, libminc2_simple), Cint, (minc2_file_handle, minc2_info_iterator_handle), h, group_it)
end

function minc2_start_attribute_iterator(h, group, it)
    ccall((:minc2_start_attribute_iterator, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cchar}, minc2_info_iterator_handle), h, group, it)
end

function minc2_iterator_group_next(it)
    ccall((:minc2_iterator_group_next, libminc2_simple), Cint, (minc2_info_iterator_handle,), it)
end

function minc2_iterator_attribute_next(it)
    ccall((:minc2_iterator_attribute_next, libminc2_simple), Cint, (minc2_info_iterator_handle,), it)
end

function minc2_iterator_group_name(it)
    ccall((:minc2_iterator_group_name, libminc2_simple), Ptr{Cchar}, (minc2_info_iterator_handle,), it)
end

function minc2_iterator_attribute_name(it)
    ccall((:minc2_iterator_attribute_name, libminc2_simple), Ptr{Cchar}, (minc2_info_iterator_handle,), it)
end

function minc2_timestamp(argc, argv)
    ccall((:minc2_timestamp, libminc2_simple), Ptr{Cchar}, (Cint, Ptr{Ptr{Cchar}}), argc, argv)
end

function minc2_xfm_allocate(h)
    ccall((:minc2_xfm_allocate, libminc2_simple), Cint, (Ptr{minc2_xfm_file_handle},), h)
end

function minc2_xfm_allocate0()
    ccall((:minc2_xfm_allocate0, libminc2_simple), minc2_xfm_file_handle, ())
end

function minc2_xfm_init(h)
    ccall((:minc2_xfm_init, libminc2_simple), Cint, (minc2_xfm_file_handle,), h)
end

function minc2_xfm_free(h)
    ccall((:minc2_xfm_free, libminc2_simple), Cint, (minc2_xfm_file_handle,), h)
end

function minc2_xfm_destroy(h)
    ccall((:minc2_xfm_destroy, libminc2_simple), Cint, (minc2_xfm_file_handle,), h)
end

function minc2_xfm_open(h, path)
    ccall((:minc2_xfm_open, libminc2_simple), Cint, (minc2_xfm_file_handle, Ptr{Cchar}), h, path)
end

function minc2_xfm_save(h, path)
    ccall((:minc2_xfm_save, libminc2_simple), Cint, (minc2_xfm_file_handle, Ptr{Cchar}), h, path)
end

function minc2_xfm_transform_point(h, in, out)
    ccall((:minc2_xfm_transform_point, libminc2_simple), Cint, (minc2_xfm_file_handle, Ptr{Cdouble}, Ptr{Cdouble}), h, in, out)
end

function minc2_xfm_inverse_transform_point(h, in, out)
    ccall((:minc2_xfm_inverse_transform_point, libminc2_simple), Cint, (minc2_xfm_file_handle, Ptr{Cdouble}, Ptr{Cdouble}), h, in, out)
end

function minc2_xfm_transform_point_vec(h, n, stride, in, out)
    ccall((:minc2_xfm_transform_point_vec, libminc2_simple), Cint, (minc2_xfm_file_handle, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), h, n, stride, in, out)
end

function minc2_xfm_inverse_transform_point_vec(h, n, stride, in, out)
    ccall((:minc2_xfm_inverse_transform_point_vec, libminc2_simple), Cint, (minc2_xfm_file_handle, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), h, n, stride, in, out)
end

function minc2_xfm_invert(h)
    ccall((:minc2_xfm_invert, libminc2_simple), Cint, (minc2_xfm_file_handle,), h)
end

function minc2_xfm_get_n_concat(h, n)
    ccall((:minc2_xfm_get_n_concat, libminc2_simple), Cint, (minc2_xfm_file_handle, Ptr{Cint}), h, n)
end

function minc2_xfm_get_n_type(h, n, xfm_type)
    ccall((:minc2_xfm_get_n_type, libminc2_simple), Cint, (minc2_xfm_file_handle, Cint, Ptr{Cint}), h, n, xfm_type)
end

function minc2_xfm_get_linear_transform(h, n, matrix)
    ccall((:minc2_xfm_get_linear_transform, libminc2_simple), Cint, (minc2_xfm_file_handle, Cint, Ptr{Cdouble}), h, n, matrix)
end

function minc2_xfm_get_grid_transform(h, n, inverted, grid_file)
    ccall((:minc2_xfm_get_grid_transform, libminc2_simple), Cint, (minc2_xfm_file_handle, Cint, Ptr{Cint}, Ptr{Ptr{Cchar}}), h, n, inverted, grid_file)
end

function minc2_xfm_append_linear_transform(h, matrix)
    ccall((:minc2_xfm_append_linear_transform, libminc2_simple), Cint, (minc2_xfm_file_handle, Ptr{Cdouble}), h, matrix)
end

function minc2_xfm_append_grid_transform(h, grid_path, inv)
    ccall((:minc2_xfm_append_grid_transform, libminc2_simple), Cint, (minc2_xfm_file_handle, Ptr{Cchar}, Cint), h, grid_path, inv)
end

function minc2_xfm_concat_xfm(h, o)
    ccall((:minc2_xfm_concat_xfm, libminc2_simple), Cint, (minc2_xfm_file_handle, minc2_xfm_file_handle), h, o)
end

function minc2_xfm_append_linear_param(h, center, translations, scales, shears, rotations)
    ccall((:minc2_xfm_append_linear_param, libminc2_simple), Cint, (minc2_xfm_file_handle, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), h, center, translations, scales, shears, rotations)
end

function minc2_xfm_extract_linear_param(h, n, center, translations, scales, shears, rotations)
    ccall((:minc2_xfm_extract_linear_param, libminc2_simple), Cint, (minc2_xfm_file_handle, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), h, n, center, translations, scales, shears, rotations)
end

function minc2_iterator_allocate0()
    ccall((:minc2_iterator_allocate0, libminc2_simple), minc2_file_iterator_handle, ())
end

function minc2_iterator_free(h)
    ccall((:minc2_iterator_free, libminc2_simple), Cint, (minc2_file_iterator_handle,), h)
end

function minc2_iterator_input_start(h, m, data_type)
    ccall((:minc2_iterator_input_start, libminc2_simple), Cint, (minc2_file_iterator_handle, minc2_file_handle, Cint), h, m, data_type)
end

function minc2_iterator_output_start(h, m, data_type)
    ccall((:minc2_iterator_output_start, libminc2_simple), Cint, (minc2_file_iterator_handle, minc2_file_handle, Cint), h, m, data_type)
end

function minc2_multi_iterator_input_start(h, m, data_type, fnum)
    ccall((:minc2_multi_iterator_input_start, libminc2_simple), Cint, (minc2_file_iterator_handle, Ptr{minc2_file_handle}, Cint, Cint), h, m, data_type, fnum)
end

function minc2_multi_iterator_output_start(h, m, data_type, fnum)
    ccall((:minc2_multi_iterator_output_start, libminc2_simple), Cint, (minc2_file_iterator_handle, Ptr{minc2_file_handle}, Cint, Cint), h, m, data_type, fnum)
end

function minc2_iterator_next(h)
    ccall((:minc2_iterator_next, libminc2_simple), Cint, (minc2_file_iterator_handle,), h)
end

function minc2_iterator_get_values(h, val)
    ccall((:minc2_iterator_get_values, libminc2_simple), Cint, (minc2_file_iterator_handle, Ptr{Cvoid}), h, val)
end

function minc2_iterator_put_values(h, val)
    ccall((:minc2_iterator_put_values, libminc2_simple), Cint, (minc2_file_iterator_handle, Ptr{Cvoid}), h, val)
end

function minc2_tags_allocate0()
    ccall((:minc2_tags_allocate0, libminc2_simple), minc2_tags_handle, ())
end

function minc2_tags_free(tags)
    ccall((:minc2_tags_free, libminc2_simple), Cint, (minc2_tags_handle,), tags)
end

function minc2_tags_load(tags, file)
    ccall((:minc2_tags_load, libminc2_simple), Cint, (minc2_tags_handle, Ptr{Cchar}), tags, file)
end

function minc2_tags_save(tags, file)
    ccall((:minc2_tags_save, libminc2_simple), Cint, (minc2_tags_handle, Ptr{Cchar}), tags, file)
end

function minc2_tags_init(tags, n_tag_points, n_volumes, have_weights, have_strucure_ids, have_patient_ids, have_labels)
    ccall((:minc2_tags_init, libminc2_simple), Cint, (minc2_tags_handle, Cint, Cint, Cint, Cint, Cint, Cint), tags, n_tag_points, n_volumes, have_weights, have_strucure_ids, have_patient_ids, have_labels)
end

function minc2_get_variable_ndims(h, path, name, ndims)
    ccall((:minc2_get_variable_ndims, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}), h, path, name, ndims)
end

function minc2_get_variable_dims(h, path, name, dims)
    ccall((:minc2_get_variable_dims, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}), h, path, name, dims)
end

function minc2_get_variable_type(h, path, name, minc2_type_)
    ccall((:minc2_get_variable_type, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}), h, path, name, minc2_type_)
end

function minc2_read_variable_raw(h, path, name, representation_type, start, count, buffer)
    ccall((:minc2_read_variable_raw, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cchar}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid}), h, path, name, representation_type, start, count, buffer)
end

function minc2_write_variable_raw(h, path, name, representation_type, start, count, buffer)
    ccall((:minc2_write_variable_raw, libminc2_simple), Cint, (minc2_file_handle, Ptr{Cchar}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid}), h, path, name, representation_type, start, count, buffer)
end

# exports
const PREFIXES = ["minc2_", "MINC2_", "miget_", "miset_"]
for name in names(@__MODULE__; all=true), prefix in PREFIXES
    if startswith(string(name), prefix)
        @eval export $name
    end
end

end # module
