using Minc2 # for reading MINC2 files


in=Minc2.read_volume("mni_icbm152_t1_tal_nlin_sym_09c.mnc")
out=Minc2.reshape_volume(in; dimrange=[[10,-20],[20,-10],[-15,15]],fill_val=10.0)
Minc2.save_volume("test_reshape.mnc",out,store=Int16)
