'''
Scripts that runs the napari plugin from the IDE.
It adds 2 aberrations to to illustrate how to include multiple aberrations to the PSF simulator. 
This script is just for demonstration. Multiple aberrations are not implemented in the napari widget.
'''
from src.napari_psf_simulator._widget import Psf_widget 
import napari

def my_aberrations_function(n1, thickness, alpha, N, M,weight):
    widget.gen.add_slab_scalar(n1, thickness, alpha) # add first aberration
    widget.gen.add_Zernike_aberration(N, M, weight) # add second aberration

class Psf_widget_with_custom_aberrations(Psf_widget):    
    def setup_aberrations(self):
        super().setup_aberrations()
        self.aberrations.add(name = 'MY_aberrations',
                        phase_aberration_function = my_aberrations_function,
                        n1 = 1.51,
                        thickness = 200.0, thickness_units = '\u03BCm', #um 
                        alpha = 0.0,
                        alpha_units = 'deg',
                        N=3, M=1, 
                        weight=0.6, weight_units = '\u03BB', #lambda
                        )

viewer = napari.Viewer()
widget = Psf_widget_with_custom_aberrations(viewer)
viewer.window.add_dock_widget(widget,
                            name = 'PSF Simulator @Polimi',
                            add_vertical_stretch = True)

napari.run()
    