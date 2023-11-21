# Simple MRI viewer

# This script displays a MINC volume using ImageView.
using Minc2 # for reading MINC2 files

using ImageView
using Gtk.ShortNames
using ImageCore
using Images
using IndirectArrays, Colors
using Statistics
using ArgParse

function mriview(path,path2,zoom::Int=5)
    a = Minc2.open_minc_file(path)
    scan, = Minc2.read_minc_volume(a,Float32)
    b = Minc2.open_minc_file(path2)
    mask, = Minc2.read_minc_volume(b,UInt8)

    println("Mask stats: mean=$(mean(mask))")
    

    #_max = maximum(scan)
    scaler = scaleminmax(minimum(scan),maximum(scan))

    #mriseg = colorview(RGB, mask)
    #mriseg = N0f8.( mask )
    scan_rgb = RGB.( scan ./ 200 )

    colors = distinguishable_colors(6)
    #mask_rgb = IndirectArray(mask, colors)

    scan_rgb[mask .> 0] .= mask_rgb[mask .> 0]
    #println(scan_rgb)
    println(typeof(scan_rgb))

    gui1 = imshow(scan_rgb;axes=(2,1),name="Axial")
    #imshow(mriseg, nothing, gui1["roi"]["zoomregion"], gui1["roi"]["slicedata"])

    gui2 = imshow(scan_rgb;axes=(3,1),name="Sagittal")
    gui3 = imshow(scan_rgb;axes=(3,2),name="Coronal")

    if (!isinteractive())
        # Create a condition object
        c = Condition()

        signal_connect(gui1["gui"]["window"], :destroy) do widget
            notify(c)
        end
        signal_connect(gui2["gui"]["window"], :destroy) do widget
            notify(c)
        end
        signal_connect(gui3["gui"]["window"], :destroy) do widget
            notify(c)
        end

        # wait for quit
        wait(c)
    end
    println("Done")

end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "in1"
            help = "Input minc file"
            required = true
        "in2"
            help = "Input minc file"
            required = true
    end
    parse_args(ARGS, s)
end


args = parse_commandline()

mriview(args["in1"],args["in1"])


