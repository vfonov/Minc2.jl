using Test, Minc2
using StatsBase

@testset "Reading volumes" begin
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


# TODO: write volume

# TODO: test 4D and 5D volume

# TODO: check xfm files