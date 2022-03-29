using Test, Minc2
using StatsBase
#using Random

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

@testset "Writing 3D volumes" begin
    mktempdir() do tmp
        in_vol,in_hdr,in_stor_hdr=Minc2.read_minc_volume_std("input/t1_z+_double_cor.mnc", Float64)

        Minc2.write_minc_volume_std(joinpath(tmp,"test1.mnc"),Float64,in_stor_hdr,in_vol)

        out_vol,out_hdr,out_stor_hdr=Minc2.read_minc_volume_std(joinpath(tmp,"test1.mnc"), Float64)

        @test in_hdr.dims == out_hdr.dims
        @test in_hdr.start == out_hdr.start
        @test in_hdr.step == out_hdr.step
        @test in_hdr.dir_cos == out_hdr.dir_cos
        @test in_vol ≈ out_vol
    end
end
# TODO: write volume

# TODO: test 4D and 5D volume

# TODO: check xfm files
