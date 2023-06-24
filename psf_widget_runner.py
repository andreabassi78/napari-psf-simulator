from src.napari_psf_simulator._widget import Psf_widget
from src.napari_psf_simulator._calculator import calculate

'''
Script that runs the napari plugin from the IDE. 
It is not executed when the plugin runs.
'''
if __name__ == '__main__':

    import napari
    viewer = napari.Viewer()
    widget = Psf_widget(viewer)
    calculator_widget = calculate()

    w1 = viewer.window.add_dock_widget(widget,
                                  name = 'PSF Simulator @Polimi',
                                  add_vertical_stretch = True)
    w2 = viewer.window.add_dock_widget(calculator_widget, name = 'Combine PSFs', add_vertical_stretch = True)
    viewer.window._qt_window.tabifyDockWidget(w1, w2)
    
    napari.run()
        