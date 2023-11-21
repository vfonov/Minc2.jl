# Simple MRI viewer

# This script displays a MINC volume using ImageView.
using Minc2 # for reading MINC2 files

using ImageView
using Gtk.ShortNames
using ImageCore
using Statistics
using ArgParse

function mriview(path, zoom::Int=5)
    a = Minc2.open_minc_file(path)
    mri, = Minc2.read_minc_volume(a,Float32)

    println(typeof(mri))

    scaler = scaleminmax(minimum(mri),maximum(mri))

    gui3 = imshow(mri;axes=(3,2),name="Coronal",scalei=scaler);
    gui1 = imshow(mri;axes=(2,1),name="Axial",scalei=scaler);
    gui2 = imshow(mri;axes=(3,1),name="Sagittal",scalei=scaler);

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
        "in"
            help = "Input minc file"
            required = true
    end
    parse_args(ARGS, s)
end

args = parse_commandline()

mriview(args["in"])
