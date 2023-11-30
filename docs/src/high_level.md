```@meta
CurrentModule = Minc2

DocTestSetup  = quote
    using Minc2
end
```

# High level functions for operations on MRI scans

## 3D Volume information

```@docs
Minc2.voxel_to_world
Minc2.world_to_voxel
Minc2.array
Minc2.history

```

## File IO functions

### Loading 3D volumes 

```@docs
Minc2.read_volume
Minc2.read_nifti_volume
```

### Saving 3D volumes

```@docs
Minc2.save_volume
Minc2.save_nifti_volume
```

### Loading transformations

```@docs
Minc2.read_transforms
Minc2.read_itk_nifti_transform
Minc2.read_ants_transform
```

### Saving transformations

```@docs
Minc2.save_transforms
Minc2.save_itk_nifti_transform
Minc2.save_itk_txt_transform
```

## 3D Volume creation

```@docs
Minc2.Volume3D
Minc2.empty_volume_like
Minc2.full_volume_like
```

## 3D Volume manipulations

```@docs
Minc2.resample_volume
Minc2.resample_volume!
Minc2.resample_grid
Minc2.resample_grid!
Minc2.crop_volume
```

## 3D Volume and geometrical transformation link

```@docs
Minc2.calculate_jacobian
Minc2.calculate_jacobian!
Minc2.GridTransform
Minc2.InverseGridTransform
Minc2.tfm_to_grid
Minc2.normalize_tfm
```
