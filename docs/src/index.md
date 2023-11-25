# MINC2 for Julia

[![CI](https://github.com/vfonov/Minc2.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vfonov/Minc2.jl/actions/workflows/CI.yml)

Read and write MINC2 files from Julia, reading and writing .xfm transformation files, 
Manipulations with volumes

## Examples

### Read a 3D minc volume, calculate mean value

Minc command:

```Shell
> mincstats -mean mni_icbm152_t1_tal_nlin_sym_09c.mnc
Mean:              29.61005195
```

Julia version:

```julia
using Minc2
using StatsBase

# read volume from minc file, represent it as array of doubles (Float64 in julia world)
icbm=Minc2.read_volume("mni_icbm152_t1_tal_nlin_sym_09c.mnc", store=Float64)
@show mean(Minc2.array(icbm))
```

Output `mean(Minc2.array(icbm)) = 29.61005194874031`

```@docs
Minc2.read_volume
Minc2.array(vol::Minc2.Volume3D)
```


### Read a 3D volume with real values, and another one with mask labels and show statistics per label

```julia
using Minc2
using StatsBase

# read T1w scan
icbm=Minc2.read_volume("mni_icbm152_t1_tal_nlin_sym_09c.mnc", store=Float64)
# read label mask, represent it as array of bytes 
lab=Minc2.read_volume("mni_icbm152_cls_tal_nlin_sym_09c.mnc", store=UInt8)

# calculate image statistics per label
for (i,l) in [1=>"CSF",2=>"GM",3=>"WM"]
    println("$(l):$(mean(Minc2.array(icbm)[ Minc2.array(lab).==i ]))")
end
```

Output:

```
CSF:31.263195070792978
GM:64.88399427312729
WM:84.74593912912998
```


### Apply transformation stored in .xfm file to a minc volume, save result

Create transformation with minc command:`param2xfm -rotations 30 0 0 rotate.xfm`

```julia
using Minc2

# read T1w scan
icbm=Minc2.read_volume("mni_icbm152_t1_tal_nlin_sym_09c.mnc", store=Float64)

# read transformation
tfm=Minc2.load_transforms("rotate.xfm")

# apply transformation to the volume
transformed_icbm=Minc2.resample_volume(icbm; tfm, order=2, fill=0.0)

# write the resulting volume into file, append history
# the resulting volume will be stored with unsigned short minc data type
# but will preserve floating point vales used fixe-point arithmetics
Minc2.save_volume("transformed_icbm.mnc",transformed_icbm, store=UInt16, history="Julia example")
```

```@docs
Minc2.load_transforms(::String)
Minc2.resample_volume
Minc2.save_volume
```


### Integrate jacobians per ROI, based on a transformation in .xfm file

```julia
using Minc2

# read ROI labels
rois = Minc2.read_volume("rois.mnc", store=UInt8)

# calculate voxel volume 
start, step, dir_cos = Minc2.decompose(Minc2.voxel_to_world(rois))
voxel_volume = abs(reduce(*,step,init=1))

# read transformation
tfm = Minc2.load_transforms("nonlinear.xfm")

# define sampling parameters for jacobian field
jacobians = Minc2.empty_volume_like(lab, store=Float64)

# calculate per-voxel jacobians
Minc2.calculate_jacobian!(tfm, jacobians)

# Integrate per ROI
for (i,l) in [1=>"ROI-1",2=>"ROI-2",3=>"ROI-3"]
    println("$(l):$( sum((Minc2.array(rois) .== i) .* Minc2.array(jac) )*voxel_volume)")
end
```

```@docs
Minc2.decompose
Minc2.empty_volume_like
Minc2.calculate_jacobian!
```


## Graphic examples, using Makie

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
![MNI-ICBM152](https://github.com/vfonov/Minc2.jl/blob/main/docs/src/assets/mni_icbm152_segmentation.png?raw=true)


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
![MNI-ICBM152](https://github.com/vfonov/Minc2.jl/blob/main/docs/src/assets/mni_icbm152_gm.png?raw=true)

```@docs
Minc2.draw_outline_with_heatmap
```


## More examples

See `examples` directory for more examples:

- `bloor.jl` - Apply Blooring kernel to minc volume, equivalent to mincbloor
- `copy.jl` - Create a minc volume with the same dimensions as another one, copy voxels
- `draw_heatmap.jl` - Show MNI-ICBM152 template contours with GM proability map, as above
- `draw_seg.jl` - Show MNI-ICBM152 template contours with tissue masks overlays, as above
- `list_meta.jl` - List all metadata stored with minc file, similar to `mincheader`
- `mriview_dual.jl` - Show two minc volumes side-by-side, interactively
- `mriview.jl` -  Show minc volume, interactively
- `multilabel_overlap.jl` - calculate voxel-wise generalized overlap coeffecient , method from http://dx.doi.org/10.1109/TMI.2006.880587
- `resample_grid.jl` - apply transformation to _grid file
- `resample_itk.jl` - resample nifti file with transformations produces by ANTS (multiple transformations can be applied in one go), equivalent to `antsApplyTransforms`
- `resample.jl`  - apply transformation to minc volume, equivalent to `mincresample`
- `xfm_to_jacobian.jl` - calculate jacobian determinant field, for a given .xfm file
- `mincinfo.jl` -  show information about mincfile(s) and store it in csv file, equivalent to `mincinfo`


## Documentation for MINC

Documentations for underlying libminc and minc tools is available at https://en.wikibooks.org/wiki/MINC


