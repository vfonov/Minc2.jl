using Test, Minc2
using StatsBase
using LinearAlgebra
using StaticArrays

@testset "Reading 3D volumes" begin
    for i in [
            "input/t1_z+_byte_cor_2.mnc",
            "input/t1_z+_byte_cor_3.mnc",
            "input/t1_z+_byte_cor.mnc",
            "input/t1_z-_byte_cor.mnc",
            "input/t1_z+_byte_sag.mnc",
            "input/t1_z-_byte_sag.mnc",
            "input/t1_z+_byte_trans.mnc",
            "input/t1_z-_byte_trans.mnc",
            
            "input/t1_z+_ubyte_cor.mnc",
            "input/t1_z-_ubyte_cor.mnc",
            "input/t1_z+_ubyte_sag.mnc",
            "input/t1_z-_ubyte_sag.mnc",
            "input/t1_z+_ubyte_trans.mnc",
            "input/t1_z-_ubyte_trans.mnc",
            "input/t1_z+_ubyte_yxz_nonorm.mnc",
            "input/t1_z+_ubyte_yxz_nonorm_noParams.mnc",
                
            "input/t1_z+_float_cor.mnc",
            "input/t1_z-_float_cor.mnc",
            "input/t1_z+_float_sag.mnc",
            "input/t1_z-_float_sag.mnc",
            "input/t1_z+_float_trans.mnc",
            "input/t1_z-_float_trans.mnc",
            "input/t1_z+_float_yxz_nonorm.mnc",
            "input/t1_z+_float_yxz_nonorm_noParams.mnc",
            "input/t1_z+_float_yxz_norm.mnc",

            "input/t1_z+_double_cor.mnc",
            "input/t1_z-_double_cor.mnc",
            "input/t1_z+_double_sag.mnc",
            "input/t1_z-_double_sag.mnc",
            "input/t1_z+_double_trans.mnc",
            "input/t1_z-_double_trans.mnc",

            "input/t1_z+_long_cor.mnc",
            "input/t1_z-_long_cor.mnc",
            "input/t1_z+_long_sag.mnc",
            "input/t1_z-_long_sag.mnc",
            "input/t1_z+_long_trans.mnc",
            "input/t1_z-_long_trans.mnc",
            "input/t1_z+_ulong_cor.mnc",
            "input/t1_z-_ulong_cor.mnc",
            "input/t1_z+_ulong_sag.mnc",
            "input/t1_z-_ulong_sag.mnc",
            "input/t1_z+_ulong_trans.mnc",
            "input/t1_z-_ulong_trans.mnc",

            "input/t1_z+_short_cor.mnc",
            "input/t1_z-_short_cor.mnc",
            "input/t1_z+_short_sag.mnc",
            "input/t1_z-_short_sag.mnc",
            "input/t1_z+_short_trans.mnc",
            "input/t1_z-_short_trans.mnc",
            "input/t1_z+_ushort_cor.mnc",
            "input/t1_z-_ushort_cor.mnc",
            "input/t1_z+_ushort_sag.mnc",
            "input/t1_z-_ushort_sag.mnc",
            "input/t1_z+_ushort_trans.mnc",
            "input/t1_z-_ushort_trans.mnc",
            ]

        @testset "Testing $(i)" begin
            vol,hdr,store_hdr = Minc2.read_minc_volume_std(i,Float64)
            @test hdr.dims == [30,40,10]
            @test size(vol) == (30,40,10)
            # TODO: test hdr too
            if contains(i,"_byte_") || contains(i,"_ubyte_")
                @test mean(vol) ≈ 35.63256359 atol=0.5
            elseif contains(i,"_short_") || contains(i,"_ushort_")
                @test mean(vol) ≈ 35.63256359 atol=0.1
            elseif contains(i,"_long_") || contains(i,"_ulong_")
                @test mean(vol) ≈ 35.63256359 atol=0.1
            elseif contains(i,"_float_") 
                @test mean(vol) ≈ 35.63256359 atol=0.01
            else
                @test mean(vol) ≈ 35.63256359 atol=0.000001
            end
        end
    end
end

@testset "Reading 4D volume, raw" begin
    vol,hdr = Minc2.read_minc_volume_raw("input/dti_sample.mnc",Float64)

    @test hdr.dims == [8,8,8,32]

    @test size(vol) == (8,8,8,32)

    expected_means=[ # calculated with mincstats
        62.40276136
        22.04154596
        18.43924784
        18.13223724
        19.55196869
        17.26168143
        18.46567965
        24.67909026
        16.74767032
        22.37696711
        18.43018655
        18.98854761
        25.19646741
        23.33803361
        25.03639432
        21.22607346
        28.81376303
        24.05946873
        25.20108933
        22.81017135
        21.76358845
        23.59570726
        30.45337502
        27.16038751
        27.85480669
        27.18367963
        18.87198478
        25.29296543
        17.92635387
        28.95309226
        18.43434116
        21.82324256
    ]
    means=vec(mean(vol,dims=(1,2,3)))

    @test length(means) == length(expected_means)

    @test means ≈ expected_means atol=0.001 # ??
end

@testset "Writing 3D volume in Float64" begin
    mktempdir() do tmp
        in_vol,in_hdr,in_stor_hdr=Minc2.read_minc_volume_std("input/t1_z+_double_cor.mnc", Float64)

        Minc2.write_minc_volume_std(joinpath(tmp,"test1.mnc"), Float64, in_stor_hdr, in_vol)

        out_vol,out_hdr,out_stor_hdr = Minc2.read_minc_volume_std(joinpath(tmp,"test1.mnc"), Float64)

        @test in_hdr.dims == out_hdr.dims
        @test in_hdr.start == out_hdr.start
        @test in_hdr.step == out_hdr.step
        @test in_hdr.dir_cos == out_hdr.dir_cos
        @test in_vol ≈ out_vol
    end
end

@testset "Writing 3D volume in short" begin
    mktempdir() do tmp
        in_vol,in_hdr,in_stor_hdr=Minc2.read_minc_volume_std("input/t1_z+_double_cor.mnc", Float64)

        Minc2.write_minc_volume_std(joinpath(tmp,"test1.mnc"), Int16, in_stor_hdr, in_vol)

        out_vol,out_hdr,out_stor_hdr = Minc2.read_minc_volume_std(joinpath(tmp,"test1.mnc"), Float64)

        @test in_hdr.dims == out_hdr.dims
        @test in_hdr.start == out_hdr.start
        @test in_hdr.step == out_hdr.step
        @test in_hdr.dir_cos == out_hdr.dir_cos
        @test in_vol ≈ out_vol atol=0.1
    end
end


@testset "Re-generate header" begin
    in_vol,in_hdr,in_stor_hdr=Minc2.read_minc_volume_std("input/t1_z+_double_cor.mnc", Float64)

    v2w=Minc2.voxel_to_world(in_hdr)
    new_hdr=Minc2.create_header_from_v2w(size(in_vol),v2w)

    @test in_hdr.start ≈ new_hdr.start atol=1e-4 
    @test in_hdr.step ≈ new_hdr.step atol=1e-6 
    @test in_hdr.dir_cos ≈ new_hdr.dir_cos atol=1e-6 
    @test in_hdr.dims == new_hdr.dims
end


@testset "read and write ITK transforms" begin
    mktempdir() do tmp
        # linear transform
        xfm = Minc2.read_itk_txt_transform("input/ants_linear.txt")
        @test xfm isa Minc2.AffineTransform

        Minc2.save_itk_txt_transform(joinpath(tmp,"test_ants_linear.txt"),xfm)
        xfm2 = Minc2.read_itk_txt_transform(joinpath(tmp,"test_ants_linear.txt"))
        @test xfm2 isa Minc2.AffineTransform
        @test xfm.rot ≈ xfm2.rot
        @test xfm.shift ≈ xfm2.shift

        # nonlinear warp
        grid_xfm = Minc2.read_itk_nifti_transform("input/ADNI_fixed_MNI-ICBM152_moving_setting_is_fastfortesting1Warp.nii.gz")
        @test grid_xfm isa Minc2.GridTransform
        Minc2.save_itk_nifti_transform(joinpath(tmp,"test_ants_warp.nii.gz"),grid_xfm)
        grid2_xfm = Minc2.read_itk_nifti_transform(joinpath(tmp,"test_ants_warp.nii.gz"))
        @test grid2_xfm isa Minc2.GridTransform

        @test Minc2.voxel_to_world(grid_xfm).rot   ≈  Minc2.voxel_to_world(grid2_xfm).rot
        @test Minc2.voxel_to_world(grid_xfm).shift ≈  Minc2.voxel_to_world(grid2_xfm).shift
        @test grid_xfm.vector_field  ≈  grid2_xfm.vector_field
    end
end



@testset "read and write minc .xfm transforms" begin
    mktempdir() do tmp
        # linear transform
        xfms = Minc2.load_transforms("input/linear.xfm")
        @test length(xfms)==1
        @test xfms[1] isa Minc2.AffineTransform

        Minc2.save_transforms(joinpath(tmp,"test_linear.xfm"),xfms)
        xfms2 = Minc2.load_transforms(joinpath(tmp,"test_linear.xfm"))
        @test length(xfms2)==1
        @test xfms[1] isa Minc2.AffineTransform

        @test xfms[1].rot ≈ xfms2[1].rot
        @test xfms[1].shift ≈ xfms2[1].shift

        # nonlinear warp
        grid_xfms = Minc2.load_transforms("input/test_nonlinear.xfm")
        @test length(grid_xfms)==1
        @test grid_xfms[1] isa Minc2.GridTransform

        Minc2.save_transforms(joinpath(tmp,"test_nonlinear.xfm"),grid_xfms)
        grid2_xfms = Minc2.load_transforms(joinpath(tmp,"test_nonlinear.xfm"))
        @test length(grid2_xfms)==1
        @test grid2_xfms[1] isa Minc2.GridTransform

        @test Minc2.voxel_to_world(grid_xfms[1]).rot   ≈  Minc2.voxel_to_world(grid2_xfms[1]).rot
        @test Minc2.voxel_to_world(grid_xfms[1]).shift ≈  Minc2.voxel_to_world(grid2_xfms[1]).shift
        @test grid_xfms[1].vector_field  ≈  grid2_xfms[1].vector_field
    end
end



@testset "Test Transforms" begin
    ### TODO: check that  all trasformation works as expected

    @testset "IdentityTransform" begin
        xfm = Minc2.IdentityTransform()
        @test Minc2.transform_point(xfm, SA_F64[1, 2, 3]) ≈ SA_F64[1, 2, 3]
        @test Minc2.transform_point(Minc2.inv(xfm),SA_F64[1, 2, 3]) ≈ SA_F64[1, 2, 3]
    end

    @testset "AffineTransform" begin
        xfm = Minc2.AffineTransform( [1.0 0 0;0 1 0;0 0 1], [1.0 1.0 1.0])
        @test Minc2.transform_point(xfm, SA_F64[1, 2, 3]) ≈ SA_F64[2, 3, 4]
        @test Minc2.transform_point(Minc2.inv(xfm),SA_F64[1, 2, 3]) ≈ SA_F64[0, 1, 2]

        xfm = Minc2.AffineTransform( [0.0 1.0 0;1.0 0.0 0;0 0 1], [0.0 0.0 0.0])
        @test Minc2.transform_point(xfm, SA_F64[1, 2, 3]) ≈ SA_F64[2, 1, 3]
        @test Minc2.transform_point(Minc2.inv(xfm),SA_F64[1, 2, 3]) ≈ SA_F64[2, 1, 3]
    end

    # TODO: come up with test for nonlinear transforms
end