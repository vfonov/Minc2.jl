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

icbm=Minc2.read_volume("mni_icbm152_t1_tal_nlin_sym_09c.mnc", store=Float64)
@show mean(Minc2.array(icbm))
```
Should print `mean(Minc2.array(icbm)) = 29.61005194874031`


### Read a 3D volume with real values, and another one with mask labels and show statistics per label

```
using Minc2
using StatsBase

icbm=Minc2.read_volume("mni_icbm152_t1_tal_nlin_sym_09c.mnc", store=Float64)
lab=Minc2.read_volume("mni_icbm152_t1_tal_nlin_sym_09c_cls.mnc", store=Uint8)

@show mean(Minc2.array(icbm))
```


See `examples` directory for more examples