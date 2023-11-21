# Show MNI-ICBM152 template contours with GM proability map

using CairoMakie
using Colors
using Minc2 


anat=Minc2.read_volume("/data/data01/vfonov/models/icbm152_model_09c/mni_icbm152_t1_tal_nlin_sym_09c.mnc",store=Float64)
gm=Minc2.read_volume("/data/data01/vfonov/models/icbm152_model_09c/mni_icbm152_gm_tal_nlin_sym_09c.mnc",store=Float64)

fig = Figure()
gc = fig[1, 1] = GridLayout()

# set background to transparent
Minc2.array(gm)[Minc2.array(gm) .< 1e-6] .= NaN

Minc2.draw_outline_with_heatmap(gc, anat, gm, 
    heat_limits=(0.0,1.0),cmap=:turbo,
    nslices = 5)

resize_to_layout!(fig)

save("mni_icbm152_gm.png", fig, px_per_unit = 1)
