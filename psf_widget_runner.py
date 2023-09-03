from src.napari_psf_simulator._widget import Psf_widget
from src.napari_psf_simulator._calculator_widget import calculate, convolution
from src.napari_psf_simulator._contrast_widget import set_contrast
'''
Script that runs the napari plugin from the IDE. 
It is not executed when the plugin runs.
'''
if __name__ == '__main__':

    import napari
    viewer = napari.Viewer()
    psf_widget = Psf_widget(viewer)
    calculator_widget = calculate()
    convolution_widget = convolution()
    contrast_widget = set_contrast()


    w1 = viewer.window.add_dock_widget(psf_widget,
                                    name = 'PSF Simulator',add_vertical_stretch = True)
    w2 = viewer.window.add_dock_widget(calculator_widget,
                                    name = 'Combine PSFs', add_vertical_stretch = True)
    # w3 = viewer.window.add_dock_widget(convolution_widget,
    #                                 name = 'Propagate', add_vertical_stretch = True)
    w4 = viewer.window.add_dock_widget(contrast_widget,
                                    name = 'Change contrast', add_vertical_stretch = True)
    viewer.window._qt_window.tabifyDockWidget(w1,w2)
    napari.run()
        