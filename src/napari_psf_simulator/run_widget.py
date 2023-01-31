from ._widget import Psf_widget

if __name__ == '__main__':
    import napari
    viewer = napari.Viewer()
    widget = Psf_widget(viewer)
    viewer.window.add_dock_widget(widget,
                                  name = 'PSF Simulator @Polimi',
                                  add_vertical_stretch = True)
    napari.run() 