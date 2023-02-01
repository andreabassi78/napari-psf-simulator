from napari_psf_simulator._widget import Psf_widget 
'''
Scripts that runs the napari plugin from the IDE
'''


import napari
viewer = napari.Viewer()
widget = Psf_widget(viewer)
viewer.window.add_dock_widget(widget,
                            name = 'PSF Simulator @Polimi',
                            add_vertical_stretch = True)
napari.run()
    