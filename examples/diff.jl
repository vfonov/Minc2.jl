using NIfTI
i1=niread("test_warp_ants.nii.gz")
i2=niread("std_warp_ants.nii.gz") #  "std_sub-1053744_ses-2_t1w.nii.gz"
o=NIVolume(i1.header,i1.extensions,i1.raw-i2.raw)
niwrite("diff_warp.nii.gz",o)
@info sum(abs.(i1.raw-i2.raw))/length(i1.raw)
