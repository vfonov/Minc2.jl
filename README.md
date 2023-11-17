# MINC2 for Julia

Read and write MINC2 files from Julia, reading and writing .xfm transformation files, 
Manipulations with volumes


## Examples


### Read a 3D minc volume, calculate mean value

Minc command: 
``` > mincstats -mean mni_icbm152_t1_tal_nlin_sym_09c.mnc
Mean:              29.61005195
```

Julia version:
```
using Minc2
using StatsBase

# read volume from minc file, represent it as array of doubles (Float64 in julia world)
icbm=Minc2.read_volume("mni_icbm152_t1_tal_nlin_sym_09c.mnc", store=Float64)
@show mean(Minc2.array(icbm))
```

Output `mean(Minc2.array(icbm)) = 29.61005194874031`


### Read a 3D volume with real values, and another one with mask labels and show statistics per label

```
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

to create transformation:`param2xfm -rotations 30 0 0 rotate.xfm`
```
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

### Integrate jacobians per ROI, based on a transformation in .xfm file
```
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


### More examples 

See `examples` directory for more examples