from _widget import Psf_widget 
import sys
print(__name__)


import napari
viewer = napari.Viewer()
widget = Psf_widget(viewer)
viewer.window.add_dock_widget(widget,
                                name = 'PSF Simulator @Polimi',
                                add_vertical_stretch = True)

napari.run()
    