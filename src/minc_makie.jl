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