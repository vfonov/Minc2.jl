# Simple MRI viewer

# This script displays a MINC volume using ImageView.
using Minc2 # for reading MINC2 files

using ImageView
using Gtk.ShortNames
using ImageCore
using Statistics

function mriview(path, zoom::Int=5)
    a = Minc2.open_minc_file(path)
    hdr, mri=Minc2.read_minc_volume(a)

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


mriview("/home/vfonov/mni/icbm152_model_09c/mni_icbm152_t1_tal_nlin_asym_09c.mnc")
