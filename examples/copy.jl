using Minc2 # for reading MINC2 files


a = Minc2.open_minc_file("/home/vfonov/mni/icbm152_model_09c/mni_icbm152_t1_tal_nlin_asym_09c.mnc")
mri, hdr, store_hdr=Minc2.read_minc_volume(a,Float32)
Minc2.close_minc_file(a)

println("store_hdr:",store_hdr)
println("hdr:",hdr)

b = Minc2.define_minc_file(store_hdr, UInt16, Float32)
Minc2.create_minc_file(b,"test.mnc")
Minc2.write_minc_volume(b,mri)
Minc2.close_minc_file(b)
