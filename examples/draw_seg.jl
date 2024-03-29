# Show MNI-ICBM152 template contours with tissue masks overlays

using CairoMakie
using Colors

using Minc2 # for reading MINC2 files


anat=Minc2.read_volume("/data/data01/vfonov/models/icbm152_model_09c/mni_icbm152_t1_tal_nlin_sym_09c.mnc",store=Float64)
seg=Minc2.read_volume("/data/data01/vfonov/models/icbm152_model_09c/mni_icbm152_cls_tal_nlin_sym_09c.mnc",store=UInt8)

fig = Figure()
gc = fig[1, 1] = GridLayout()

Minc2.draw_outline_with_labels(gc, anat, seg, 
    labels=Dict([1=>"CSF",2=>"GM",3=>"WM"]),
    nslices = 5)

resize_to_layout!(fig)

save("mni_icbm152_segmentation.png", fig, px_per_unit = 1) 


