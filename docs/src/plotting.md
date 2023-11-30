```@meta
CurrentModule = Minc2

DocTestSetup  = quote
    using Minc2
end
```

# High level functions for plotting 3D volumes

These are helper functions for plotting `Volume3D` objects using Makie. To make them available one must first load `Makie` or equivalent (i.e `CairoMakie`)

### Show MNI-ICBM152 template contours with tissue masks overlays

```julia
using CairoMakie
using Minc2

# read T1w scan
icbm=Minc2.read_volume("mni_icbm152_t1_tal_nlin_sym_09c.mnc", store=Float64)
# read label mask, represent it as array of bytes 
lab=Minc2.read_volume("mni_icbm152_cls_tal_nlin_sym_09c.mnc", store=UInt8)


fig = Figure()
gc = fig[1, 1] = GridLayout()

Minc2.draw_outline_with_labels(gc, icbm, lab, 
    labels=Dict([1=>"CSF",2=>"GM",3=>"WM"]),
    nslices = 5)

resize_to_layout!(fig)

save("mni_icbm152_segmentation.png", fig, px_per_unit = 1)
```

Will produce
![MNI-ICBM152 Segmentation](https://github.com/vfonov/Minc2.jl/blob/main/docs/src/assets/mni_icbm152_segmentation.png?raw=true)

```@docs
Minc2.draw_outline_with_labels
```

### Show MNI-ICBM152 template contours with GM proability map

```julia
using CairoMakie
using Minc2

# read T1w scan
anat=Minc2.read_volume("mni_icbm152_t1_tal_nlin_sym_09c.mnc",store=Float64)
# read field
gm=Minc2.read_volume("mni_icbm152_gm_tal_nlin_sym_09c.mnc",store=Float64)

fig = Figure()
gc = fig[1, 1] = GridLayout()

# set background to transparent
Minc2.array(gm)[Minc2.array(gm) .< 1e-6] .= NaN

Minc2.draw_outline_with_heatmap(gc, anat, gm, 
    heat_limits=(0.0,1.0),cmap=:turbo,
    nslices = 5)

resize_to_layout!(fig)

save("mni_icbm152_gm.png", fig, px_per_unit = 1) 
```

Will produce
![MNI-ICBM152 GM](https://github.com/vfonov/Minc2.jl/blob/main/docs/src/assets/mni_icbm152_gm.png?raw=true)

```@docs
Minc2.draw_outline_with_heatmap
```

### Show two-tailed t-statistics on DBM style analysis



```julia
using CairoMakie
using Minc2

# read T1w scan and crop black area around ROI
anat=Minc2.crop_volume(
    Minc2.read_volume("mni_icbm152_t1_tal_nlin_sym_09c.mnc",store=Float64),
    [[15,15],[20,20],[20,20]])

# Read brain ROI, should be the same sampling as anat
mask=Minc2.crop_volume(
    Minc2.read_volume("mni_icbm152_t1_tal_nlin_sym_09c_mask.mnc",store=UInt8),
    [[15,15],[20,20],[20,20]])

# read statistics map, does not have to be the same resolution as anat
vol_t=Minc2.read_volume("t_AD_vs_CN_2mm.mnc",store=Float64)

# statistical limits (two-tailed)
stat_limits  = [2.9 , 5]

# start new figure
fig = Figure()

g1  = f[1,1] = GridLayout()
Label(f[1,1,Top()],"AD vs CN-")

draw_outline_with_heatmap_symmetric(g1, anat, mask, vol_stat;
    heat_limits=stat_limits,
    show_colorbar=true,
    cmap_pos=:YlOrRd_9,cmap_neg=:YlGnBu_9,
    nslices=5,
    colorbar_label="t-statistics")

resize_to_layout!(fig)

save("t-statistics.png", fig, px_per_unit = 1) 
```

![T-statistics](https://github.com/vfonov/Minc2.jl/blob/main/docs/src/assets/t-statistics.png?raw=true)

```@docs
Minc2.draw_outline_with_heatmap_symmetric
```
