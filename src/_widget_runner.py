from napari_psf_simulator._widget import Psf_widget 
'''
Script that runs the napari plugin from the IDE. 
It is not executed when the plugin runs.
'''


import napari
viewer = napari.Viewer()
widget = Psf_widget(viewer)
viewer.window.add_dock_widget(widget,
                            name = 'PSF Simulator @Polimi',
                            add_vertical_stretch = True)
napari.run()
    