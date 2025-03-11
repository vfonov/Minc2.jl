using Test, Minc2
using StatsBase
using LinearAlgebra
using StaticArrays
using DelimitedFiles
using Tables

@testset "Low level io" begin
    @testset "Try open missing file" begin
        @test_throws Minc2.Minc2Error Minc2.open_minc_file("input/missing.mnc")
    end

    @testset "Reading attributes" begin
        h=Minc2.open_minc_file("input/t1_z+_byte_cor.mnc")
        @test length(Minc2.read_attribute(h,"","history")) == 777
        @test length(Minc2.read_attribute(h,"","ident")) == 45
        ### TODO: add tests for more attributes
        Minc2.close_minc_file(h)
    end

    @testset "Testing coronal byte" begin
        h=Minc2.open_minc_file("input/t1_z+_byte_cor.mnc")

        @test Minc2.store_type(h) == Type{Int8}
        @test Minc2.representation_type(h) == Type{Float32}
        @test Minc2.ndim(h) == 3

        sh = Minc2.store_header(h)

        @test sh.dims == [30, 10, 40]
        @test sh.start == [-18.000000000000018, -9.999999999999998, 59.999999999999986]
        @test sh.step == [1.0,1.2,1.0]
        @test isapprox(sh.dir_cos, [ 0.998021 0.052304 -0.034899; 0.03576 -0.015602 0.999239;-0.05172 0.998509 0.017442],rtol=1e-6)
        @test sh.axis == [Minc2.DIM_X, Minc2.DIM_Z, Minc2.DIM_Y]

        # instruct minc that we want to read in "standard" order
        Minc2.setup_standard_order( h )
        rh = Minc2.representation_header(h)

        @test rh.dims == [30, 40, 10]
        @test rh.start == [-18.000000000000018, 59.999999999999986, -9.999999999999998]
        @test rh.step == [1.0,1.0,1.2 ]
        @test isapprox(rh.dir_cos, [ 0.998021 0.052304 -0.034899; -0.05172 0.998509 0.017442; 0.03576 -0.015602 0.999239],rtol=1e-6)
        @test rh.axis == [Minc2.DIM_X, Minc2.DIM_Y, Minc2.DIM_Z]

        Minc2.close_minc_file(h)
    end

    @testset "Testing sagittal long" begin
        h = Minc2.open_minc_file("input/t1_z-_long_sag.mnc")
        
        @test Minc2.store_type(h) == Type{Int32}
        @test Minc2.representation_type(h) == Type{Float32}
        @test Minc2.ndim(h) == 3

        sh = Minc2.store_header(h)

        @test sh.dims == [40, 10, 30]
        @test sh.start == [59.999999999999986, 0.8000000000000007, -18.000000000000018]
        @test sh.step == [1.0,-1.2,1.0]
        @test isapprox(sh.dir_cos, [-0.05172 0.998509 0.017442; 0.03576 -0.015602 0.999239;0.998021 0.052304 -0.034899],rtol=1e-6)
        @test sh.axis == [Minc2.DIM_Y, Minc2.DIM_Z, Minc2.DIM_X]

        # instruct minc that we want to read in "standard" order, should be same as before
        Minc2.setup_standard_order( h )
        rh = Minc2.representation_header(h)

        @test rh.dims == [30, 40, 10]
        @test rh.start == [-18.000000000000018, 59.999999999999986, -9.999999999999998]
        @test rh.step == [1.0,1.0,1.2 ]
        @test isapprox(rh.dir_cos, [ 0.998021 0.052304 -0.034899; -0.05172 0.998509 0.017442; 0.03576 -0.015602 0.999239],rtol=1e-6)
        @test rh.axis == [Minc2.DIM_X, Minc2.DIM_Y, Minc2.DIM_Z]

        Minc2.close_minc_file(h)
    end

    @testset "Testing volume types, Char" begin
        h = Minc2.open_minc_file("input/3DCharImage.mnc")
        @test Minc2.store_type(h) == Type{Int8}
        @test Minc2.representation_type(h) == Type{Int8}
        Minc2.close_minc_file(h)
    end

    @testset "Testing volume types, UShort" begin
        h = Minc2.open_minc_file("input/3DUShortImage.mnc")
        @test Minc2.store_type(h) == Type{UInt16}
        @test Minc2.representation_type(h) == Type{UInt16}
        Minc2.close_minc_file(h)
    end

    @testset "Testing volume types, Int" begin
        h = Minc2.open_minc_file("input/3DIntImage.mnc")
        @test Minc2.store_type(h) == Type{Int32}
        @test Minc2.representation_type(h) == Type{Int32}
        Minc2.close_minc_file(h)
    end
end

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

@testset "Testing reading wrong voxel type" begin
    vol,hdr,store_hdr = Minc2.read_minc_volume_std("input/t1_z-_long_sag.mnc",Int16)
    @test eltype(vol) == Int16
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


@testset "Test Transforms" begin
    ### TODO: check that  all transformations works as expected

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


@testset "Read and Write Linear ITK transforms" begin
    mktempdir() do tmp
        # linear transform
        xfm = Minc2.read_itk_txt_transform("input/ants_linear.txt")
        @test xfm isa Minc2.AffineTransform

        # check if transformation is correct:
        # antsApplyTransformsToPoints -d 3 --precision 1 --input input/ants_points_input.csv --output input/ants_points_linear_output.csv --transform input/ants_linear.txt
        input_points=readdlm("input/ants_points_input.csv",',',header=true)[1] # x,y,z,t,label,comment 
        ref_output_points=readdlm("input/ants_points_linear_output.csv",',',header=true)[1] # x,y,z,t,label,comment 

        # check all points
        for i in axes(input_points,1)
            @test Minc2.transform_point(xfm, SVector{3,Float64}(input_points[i,1:3])) ≈ SVector{3,Float64}(ref_output_points[i,1:3]) atol=1e-6
        end

        Minc2.save_itk_txt_transform(joinpath(tmp,"test_ants_linear.txt"),xfm)
        xfm2 = Minc2.read_itk_txt_transform(joinpath(tmp,"test_ants_linear.txt"))
        @test xfm2 isa Minc2.AffineTransform
        @test xfm.rot ≈ xfm2.rot
        @test xfm.shift ≈ xfm2.shift
    
    end
end

@testset "Read and Write NonLinear ITK transforms" begin
    mktempdir() do tmp
        # nonlinear warp
        grid_xfm = Minc2.read_itk_nifti_transform("input/ADNI_fixed_MNI-ICBM152_moving_setting_is_fastfortesting1Warp.nii.gz")
        @test grid_xfm isa Minc2.GridTransform

        # antsApplyTransformsToPoints -d 3 --precision 1 --input input/ants_points_input.csv --output input/ants_points_nonlinear_output.csv --transform input/ADNI_fixed_MNI-ICBM152_moving_setting_is_fastfortesting1Warp.nii.gz
        input_points=readdlm("input/ants_points_input.csv",',',header=true)[1] # x,y,z,t,label,comment 
        ref_output_points=readdlm("input/ants_points_nonlinear_output.csv",',',header=true)[1] # x,y,z,t,label,comment 

        for i in axes(input_points,1)
            @test Minc2.transform_point(grid_xfm, SVector{3,Float64}(input_points[i,1:3])) ≈ SVector{3,Float64}(ref_output_points[i,1:3]) atol=1e-6
        end

        Minc2.save_itk_nifti_transform(joinpath(tmp,"test_ants_warp.nii.gz"),grid_xfm)
        grid2_xfm = Minc2.read_itk_nifti_transform(joinpath(tmp,"test_ants_warp.nii.gz"))
        @test grid2_xfm isa Minc2.GridTransform

        @test Minc2.voxel_to_world(grid_xfm).rot   ≈  Minc2.voxel_to_world(grid2_xfm).rot
        @test Minc2.voxel_to_world(grid_xfm).shift ≈  Minc2.voxel_to_world(grid2_xfm).shift
        @test grid_xfm.vector_field  ≈  grid2_xfm.vector_field
    end
end


@testset "Read and Write MINC .xfm transforms" begin
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

@testset "Test Spatial coordinates correctness" begin
    ### TODO: check that  all trasformation works as expected
    @testset "Low level interface" begin
        h = Minc2.open_minc_file("input/t1_z+_float_cor.mnc")
        @test Minc2.voxel_to_world(h, [0.0,0.0,0.0]) ≈ [-21.425182690867980995,59.125105562514228552,-8.3176836593094876093] atol = 1e-9
        @test Minc2.voxel_to_world(h, [10.0,10.0,10.0]) ≈ [-11.533050841537136222,69.446014720535146125,3.4986096455089885637] atol = 1e-9
        Minc2.close_minc_file(h)
    end

    @testset "Oblique volume" begin
        vol = Minc2.read_volume("input/t1_z+_float_cor.mnc", store=Float64)
        v2w = Minc2.voxel_to_world(vol)

        @test Minc2.transform_point(v2w, SA_F64[0,0,0] ) ≈ SA_F64[-21.425182690867980995, 59.125105562514228552, -8.3176836593094876093] atol=1e-6
        @test Minc2.transform_point(v2w, SA_F64[10,10,10] ) ≈ SA_F64[-11.533050841537136222,69.446014720535146125,3.4986096455089885637] atol=1e-6

        start, step, dir_cos = Minc2.decompose(v2w)
        @test start ≈ SA_F64[-18, 60, -10] atol=1e-4
        @test  step ≈ SA_F64[  1,  1, 1.2] atol=1e-5

        @test dir_cos[:,1] ≈ [0.998021217040696, 0.0523040113746069, -0.0348990075895229] atol=1e-6
        @test dir_cos[:,2] ≈ [-0.0517200153907152, 0.998509297133945, 0.0174420051903491] atol=1e-6
        @test dir_cos[:,3] ≈ [0.0357599860692531, -0.0156019939220494, 0.999238610734185] atol=1e-6
        
        # calculate COM in voxel coordinates
        arr = Minc2.array(vol)
        v_com = @SVector [reduce((x,y)->x+arr[y]*(y[i]-1), collect(CartesianIndices(arr));init=0.0)/sum(arr) for i in 1:3]
        @test v_com ≈ SA_F64[13.72049578,  12.72096456,  4.525147114] atol=1e-3
        
        w_com = Minc2.transform_point(v2w, v_com)
        @test w_com ≈ SA_F64[-8.195582241, 72.46002233, -3.148594157] atol=1e-3

    end

    @testset "AffineTransform" begin
        xfms = Minc2.load_transforms("input/linear.xfm")
        @test length(xfms)==1
        xfm=xfms[1]
        @test xfm isa Minc2.AffineTransform
        @test Minc2.transform_point(xfm, SA_F64[0,0,0]) ≈ SA_F64[-0.176377685791353,-1.07878519423042,-0.546123974572532] atol=1e-9
        @test Minc2.transform_point(xfm, SA_F64[0,10,0]) ≈ SA_F64[-0.072901424071292,8.88340936772789,-0.121799365233295] atol=1e-9
        @test Minc2.transform_point(xfm, SA_F64[0,0,10]) ≈ SA_F64[-0.129016039943794,-1.02947230651711,9.51151610111107] atol=1e-9
        @test Minc2.transform_point(xfm, SA_F64[10,10,10]) ≈ SA_F64[9.96594152594497,8.75495294839407,9.92949690428783]   atol=1e-9
    end

    @testset "NonLinearTransform" begin
        xfms = Minc2.load_transforms("input/test_nonlinear.xfm")
        @test length(xfms)==1
        xfm=xfms[1]
        @test xfm isa Minc2.GridTransform
        # grid interpolation algorithms in minc and in Julia are different :(
        @test Minc2.transform_point(xfm, SA_F64[0,73.4967727661133,-32.4086036682129])    ≈ SA_F64[0.151066607496551,80.9965964661819,-36.0780171315063] atol=0.2
        @test Minc2.transform_point(xfm, SA_F64[-42.2117652893066,48.4688186645508,52.8470573425293])   ≈ SA_F64[-45.8093159136667,51.4849761254315,57.8832412837079] atol=0.2
        @test Minc2.transform_point(xfm, SA_F64[43.294116973877,48.4688186645508,52.8470573425293]) ≈ SA_F64[46.4595841478604,50.8305390103988,58.2716356097526]   atol=0.2

        # test invert
        ixfm=Minc2.inv(xfm)
        @test Minc2.transform_point(ixfm, SA_F64[0,73.4967727661133,-32.4086036682129])    ≈ SA_F64[-0.18140356051401,66.6279125797094,-29.0028658475105] atol=0.3
        @test Minc2.transform_point(ixfm, SA_F64[-42.2117652893066,48.4688186645508,52.8470573425293])   ≈ SA_F64[-38.3819344872913, 45.9948061533514, 48.218509626188] atol=0.4
        @test Minc2.transform_point(ixfm, SA_F64[43.294116973877,48.4688186645508,52.8470573425293]) ≈ SA_F64[39.8429712902809,46.7448946781746,47.7003233054328]   atol=0.3

    end
end


@testset "Test Jacobian calculation" begin
    @testset "Test identity tfm to grid" begin
        xfm = Minc2.AffineTransform( [1.0 0 0;0 1 0;0 0 1], [0.0 0.0 0.0])
        dummy = Minc2.Volume3D(zeros(Float64, (20,20,20)), Minc2.AffineTransform( [0.0 1.0 0;1.0 0.0 0;0 0 1], [-10.0 -10.0 -10.0]))
        grid = Minc2.tfm_to_grid(xfm, dummy)
        @test size(Minc2.array(grid)) == (3,20,20,20)
        @test all( Minc2.array(grid) .≈ 0.0)
    end

    @testset "Test y-shift tfm to grid" begin
        xfm = Minc2.AffineTransform( [1.0 0 0;0 1 0;0 0 1], [0.0 1.0 0.0])
        dummy = Minc2.Volume3D(zeros(Float64, (20,20,20)), Minc2.AffineTransform( [0.0 1.0 0;1.0 0.0 0;0 0 1], [-10.0 -10.0 -10.0]))
        grid = Minc2.tfm_to_grid(xfm, dummy)
        @test size(Minc2.array(grid)) == (3,20,20,20)

        @test all( Minc2.array(grid)[1,:,:,:] .≈ 0.0)
        @test all( Minc2.array(grid)[2,:,:,:] .≈ 1.0)
        @test all( Minc2.array(grid)[3,:,:,:] .≈ 0.0)
    end

    @testset "Test x-shift tfm to grid" begin
        xfm = Minc2.AffineTransform( [1.0 0 0;0 1 0;0 0 1], [1.0 0.0 0.0])
        dummy = Minc2.Volume3D(zeros(Float64, (20,20,20)), Minc2.AffineTransform( [0.0 1.0 0;1.0 0.0 0;0 0 1], [-10.0 -10.0 -10.0]))
        grid = Minc2.tfm_to_grid(xfm, dummy)
        @test size(Minc2.array(grid)) == (3,20,20,20)

        @test all( Minc2.array(grid)[1,:,:,:] .≈ 1.0)
        @test all( Minc2.array(grid)[2,:,:,:] .≈ 0.0)
        @test all( Minc2.array(grid)[3,:,:,:] .≈ 0.0)
    end

    @testset "Test z-shift tfm to grid" begin
        xfm = Minc2.AffineTransform( [1.0 0 0;0 1 0;0 0 1], [0.0 0.0 1.0])
        dummy = Minc2.Volume3D(zeros(Float64, (20,20,20)), Minc2.AffineTransform( [0.0 1.0 0;1.0 0.0 0;0 0 1], [-10.0 -10.0 -10.0]))
        grid = Minc2.tfm_to_grid(xfm, dummy)
        @test size(Minc2.array(grid)) == (3,20,20,20)

        @test all( Minc2.array(grid)[1,:,:,:] .≈ 0.0)
        @test all( Minc2.array(grid)[2,:,:,:] .≈ 0.0)
        @test all( Minc2.array(grid)[3,:,:,:] .≈ 1.0)
    end


    @testset "Identity transform, unit jacobian" begin
        # identity transform
        xfm = Minc2.AffineTransform( [1.0 0 0;0 1 0;0 0 1], [0.0 0.0 0.0])
        # volume centered at origin
        jac = Minc2.Volume3D(zeros(Float64, (20,20,20)), Minc2.AffineTransform( [1.0 0.0 0;0.0 1.0 0;0 0 1], [-10.0 -10.0 -10.0]))

        Minc2.calculate_jacobian!(xfm, jac)
        @test mean(Minc2.array(jac)) ≈ 1.0 atol=1e-6
    end

    @testset "Uniform scaling jacobian" begin
        # identity transform
        xfm = Minc2.AffineTransform( [0.95 0 0;0 0.95 0;0 0 0.95], [0.0 0.0 0.0])

        # volume centered at origin
        jac = Minc2.Volume3D(zeros(Float64, (20,20,20)), Minc2.AffineTransform( [1.0 0.0 0;0.0 1.0 0;0 0 1], [-10.0 -10.0 -10.0]))
        Minc2.calculate_jacobian!(xfm, jac)
        @test mean(Minc2.array(jac)) ≈ det(xfm.rot) atol=1e-6
    end


    @testset "Z-scaling jacobian" begin
        # identity transform
        xfm = Minc2.AffineTransform( [1.0 0 0;0 1.0 0;0 0 0.9], [0.0 0.0 0.0])
        # volume centered at origin
        jac = Minc2.Volume3D(zeros(Float64, (20,20,20)), Minc2.AffineTransform( [1.0 0.0 0;0.0 1.0 0;0 0 1], [-10.0 -10.0 -10.0]))
        Minc2.calculate_jacobian!(xfm, jac)
        @test mean(Minc2.array(jac)) ≈ det(xfm.rot) atol=1e-6
    end

    @testset "X-scaling jacobian" begin
        # identity transform
        xfm = Minc2.AffineTransform( [0.9 0 0;0 1.0 0;0 0 1.0], [0.0 0.0 0.0])
        Minc2.save_transforms("test_scale_x.xfm",[xfm])

        # volume centered at origin
        jac = Minc2.Volume3D(zeros(Float64, (20,20,20)), Minc2.AffineTransform( [1.0 0.0 0;0.0 1.0 0;0 0 1], [-10.0 -10.0 -10.0]))

        Minc2.calculate_jacobian!(xfm, jac)
        @test mean(Minc2.array(jac)) ≈ det(xfm.rot) atol=1e-6
    end

    @testset "Pure rotation transform, unit jacobian" begin
        #xfm = Minc2.load_transforms("r30.xfm")[1]
        xfm = Minc2.AffineTransform( [
         1 0 0;
         0 0.866025388240814 -0.5;
         0 0.5 0.866025388240814;
        ], [0.0 0.0 0.0])
        # volume centered at origin
        jac = Minc2.Volume3D(zeros(Float64, (20,25,30)), Minc2.AffineTransform( [1.0 0.0 0;0.0 1.0 0;0 0 1], [-10.0 -12.5 -15.0]))

        Minc2.calculate_jacobian!(xfm, jac)
        @test mean(Minc2.array(jac)) ≈ det(xfm.rot) atol=1e-6
    end
end
