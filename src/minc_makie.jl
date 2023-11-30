##### visualization helpers

"""
    draw_outline_with_labels(
            layout, 
            anat, 
            seg; 
            labels, 
            cmap = :rainbow,
            levels = [20,30],
            show_colorbar = true,
            nslices = 4)

Draws a segmentation `seg` on top of `anat` in a grid layout `layout`

* `layout` is a grid layout
* `anat` is a Volume3D with anatomical image
* `seg` is a Volume3D with discrete segemntations
* `labels` is a dictionary of labels to show
* `cmap` is a color map to use
* `levels` is a list of levels of `anat` to use for contours
* `show_colorbar` is a boolean flag to show colorbar
* `nslices` is a number of slices to show
"""
function draw_outline_with_labels(
        layout, 
        anat, 
        seg; 
        labels, 
        cmap = :rainbow,
        levels = [20,30],
        show_colorbar = true,
        nslices = 4)
    
    sz = size(array(anat))
    w = max(sz...)

    relabel = Dict( collect(keys(labels)) .=>  1:length(labels)  )
    # relabel segmentation
    relabel_seg = map(x->x in keys(labels) ? relabel[x] : 0 , array(seg))
    
    # TODO: generate labels
    # TODO: add implicit background label (0)
    colors = Makie.cgrad(cmap, length(labels), categorical = true)

    # for showing ROIs
    opts_brain=(colorrange = [1, length(keys(labels))], colormap=colors,  
        interpolate=false, lowclip=Makie.colorant"transparent", 
        alpha=0.7,transparency=true)

    opt_c=(levels=levels, color=:black, alpha=0.7, transparency=true)

    # width=w, height=w, 
    aopt = (aspect=Makie.DataAspect(), 
        titlevisible =false, topspinevisible=false,yticksvisible = false,
        yticklabelsvisible = false, bottomspinevisible=false, leftspinevisible=false, 
        rightspinevisible = false, xgridvisible=false,ygridvisible=false)

    gc = layout
    baxs = []

    # show brain with ROIs
    
    for ax in 1:3, s in 1:nslices
        mask_s = selectdim(relabel_seg, ax, fld1(sz[ax]*s,nslices+1))
        anat_s = selectdim(array(anat), ax, fld1(sz[ax]*s,nslices+1)) 
        (w,h)  = size(anat_s)
        b1a = Makie.Axis(gc[ax,s]; width=w, height=h, aopt...)
        Makie.hidedecorations!(b1a, grid = false)

        push!(baxs,b1a)

        Makie.heatmap!(b1a, mask_s;opts_brain...)
        Makie.contour!(b1a, anat_s;opt_c...)
    end

    Makie.tightlimits!.(baxs)
    #colsize!(gc, 1, Aspect(1, 1.0))
    #colsize!(gc, 2, Aspect(1, 1.0))
    # add color bar
    if show_colorbar
        Makie.Colorbar( gc[:,nslices+1], colormap=colors, 
                limits = (1, length(labels)+1),
                ticks=   (collect(1:length(labels)) .+ 0.5, collect(values(labels))))
    end
    Makie.colgap!(gc, 0)
    Makie.rowgap!(gc, 0)
    return gc
end

"""
    draw_outline_with_heatmap(
        layout, anat, heat; 
        cmap=:rainbow, levels=[20, 40],
        heat_limits=nothing, 
        show_colorbar=true,
        nslices=4)

Draws a heatmap of `heat` on top of `anat` in a grid layout `layout`.

* `layout` is a grid layout
* `anat` is a Volume3D with anatomical image
* `heat` is a Volume3D with overlay heatmap
* `heat_limits` limits for the heatmap
* `cmap` is a color map to use
* `levels` is a list of levels of `anat` to use for contours
* `show_colorbar` is a boolean flag to show colorbar
* `nslices` is a number of slices to show
"""
function draw_outline_with_heatmap(
    layout, anat, heat; 
    cmap=:rainbow, levels=[20, 40],
    heat_limits=nothing, 
    show_colorbar=true,
    nslices=4)

    sz = size(Minc2.array(anat))
    w = max(sz...)

    if isnothing(heat_limits)
        heat_limits = extrema(Minc2.array(heat))
    end
    
    # for showing ROIs
    opts_brain=(colorrange = heat_limits, colormap=cmap,  
        interpolate=false, lowclip=Makie.colorant"transparent", 
        alpha=0.7,transparency=true)

    opt_c=(levels=levels, color=:black, alpha=0.7, transparency=true)

    # width=w, height=w, 
    aopt=(aspect=Makie.DataAspect(), 
        titlevisible =false, topspinevisible=false,yticksvisible = false,
        yticklabelsvisible = false, bottomspinevisible=false, leftspinevisible=false, 
        rightspinevisible = false, xgridvisible=false,ygridvisible=false)

    gc = layout
    baxs = []

    # show brain with ROIs
   
    for ax in 1:3, s in 1:nslices
        mask_s = selectdim(Minc2.array(heat), ax, fld1(sz[ax]*s,nslices+1))
        anat_s = selectdim(Minc2.array(anat), ax, fld1(sz[ax]*s,nslices+1)) 
        (w,h)  = size(anat_s)
        b1a = Makie.Axis(gc[ax,s]; width=w, height=h, aopt...)
        Makie.hidedecorations!(b1a, grid = false)

        push!(baxs,b1a)

        Makie.heatmap!(b1a, mask_s;opts_brain...)
        Makie.contour!(b1a, anat_s;opt_c...)
    end

    Makie.tightlimits!.(baxs)
    # add color bar
    if show_colorbar
        Makie.Colorbar( gc[:,nslices+1], colormap=cmap, 
                limits = heat_limits) # highclip = :cyan, lowclip = :red
    end
    Makie.colgap!(gc, 0)
    Makie.rowgap!(gc, 0)
    return gc
end

"""
    draw_outline_with_heatmap_symmetric(
        layout, anat, mask, heat; 
        cmap_pos=:YlOrRd_9,
        cmap_neg=:YlGnBu_9,
        levels=[20, 40],
        heat_limits=nothing,
        show_colorbar=true,
        under=nothing, 
        over=nothing,
        nslices=4,
        alpha=0.7,
        colorbar_label="",
        slices_gap=5,
        colorbar_gap=10)

Draws a two-sided (symmetric positive and negative) heatmap of `heat` on top of `anat` in a grid layout `layout`.

* `layout` is a grid layout
* `anat` is a Volume3D with anatomical image
* `mask` is a Volume3D with ROI
* `heat` is a Volume3D with overlay heatmap, can be arbitrary sampling
* `heat_limits` limits for the heatmap
* `cmap_pos` is a color map to use for positive values, if nothing then don't show positive values
* `cmap_neg` is a color map to use for negative values, if nothing then don't show negative values
* `levels` is a list of levels of `anat` to use for contours
* `show_colorbar` is a boolean flag to show colorbar
* `colorbar_label` is a label for the colorbar
* `nslices` is a number of slices to show
"""
function draw_outline_with_heatmap_symmetric(
    layout, anat, mask, heat; 
    cmap_pos=:YlOrRd_9,
    cmap_neg=:YlGnBu_9,
    levels=[20, 40],
    heat_limits=nothing,
    show_colorbar=true,
    under=nothing, 
    over=nothing,
    nslices=4,
    alpha=0.7,
    colorbar_label="",
    slices_gap=5,
    colorbar_gap=10)

    sz = size(Minc2.array(anat))
    w = max(sz...)

    if !isnothing(cmap_neg)
        cmap_neg_rev = Reverse(cmap_neg)
    end

    # resample heat to match sampling of mask
    heat_r = Minc2.empty_volume_like(mask; store=Float64)
    Minc2.resample_volume!(heat, heat_r)

    if isnothing(heat_limits)
        vol_sel = Minc2.array(heat_r)[Minc2.array(mask)>0]
        vol_pos = vol_sel[vol_sel .> 0]
        vol_neg = - vol_sel[vol_sel .< 0]
        heat_limits = extrema(vcat(vol_pos, vol_neg))
    end
    
    # for showing ROIs
    opts_pos=(colorrange = heat_limits, 
        colormap=cmap_pos,
        interpolate=false, lowclip=nothing, highclip=over, 
        alpha, transparency=true)

    opts_neg=(colorrange = (-heat_limits[2],-heat_limits[1]), 
        colormap=cmap_neg_rev,
        interpolate=false, highclip=nothing, lowclip=under, 
        alpha, transparency=true)

    opt_c=(levels=levels, color=:black, alpha=0.7, transparency=true)

    # width=w, height=w, 
    aopt=(aspect=Makie.DataAspect(), 
        titlevisible =false, topspinevisible=false,yticksvisible = false,
        yticklabelsvisible = false, bottomspinevisible=false, leftspinevisible=false, 
        rightspinevisible = false, xgridvisible=false,ygridvisible=false)

    gc = layout
    baxs = []

    # show brain with ROIs
   
    for ax in 1:3, s in 1:nslices
        heat_pos = deepcopy(selectdim(Minc2.array(heat_r), ax, fld1(sz[ax]*s,nslices+1)))
        heat_neg = deepcopy(heat_pos)

        # make transparent
        heat_pos[heat_pos .<  heat_limits[1]] .= NaN
        heat_neg[heat_pos .> -heat_limits[1]] .= NaN

        mask_s = selectdim(Minc2.array(mask), ax, fld1(sz[ax]*s,nslices+1))
        anat_s = selectdim(Minc2.array(anat), ax, fld1(sz[ax]*s,nslices+1))

        # make background transparent
        heat_pos[mask_s .== 0] .= NaN
        heat_neg[mask_s .== 0] .= NaN

        (w,h)  = size(anat_s)
        b1a = Makie.Axis(gc[ax,s]; width=w, height=h, aopt...)
        Makie.hidedecorations!(b1a, grid = false)

        push!(baxs,b1a)

        if !isnothing(cmap_pos)
            Makie.heatmap!(b1a, heat_pos;opts_pos...)
        end

        if !isnothing(cmap_neg)
            Makie.heatmap!(b1a, heat_neg;opts_neg...)
        end
        Makie.contour!(b1a, anat_s;opt_c...)
    end

    Makie.tightlimits!.(baxs)

    # add color bars
    if show_colorbar
        gb = gc[:,nslices+1] = GridLayout(default_rowgap=colorbar_gap, default_colgap=0)
        # separate from the main heatmaps
        if !isnothing(cmap_pos)
            Makie.Colorbar( gb[1,1], colormap=cmap_pos, 
                    limits = heat_limits, 
                    highclip = over, 
                    lowclip = nothing)
        end
        if !isnothing(cmap_neg)
            Makie.Colorbar( gb[2,1], colormap=cmap_neg_rev, 
                    limits = (-heat_limits[2],-heat_limits[1]), 
                    highclip = nothing, lowclip = under)
        end
        Makie.Label(gb[:,2], colorbar_label , rotation = pi/2)
    end
    Makie.colgap!(gc, slices_gap)
    Makie.rowgap!(gc, slices_gap)

    if show_colorbar
        Makie.colgap!(gc, nslices, colorbar_gap)
    end
    return gc
end