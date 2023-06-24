# -*- coding: utf-8 -*-
"""
Created on Sat Jun 17 11:16:58 2023

@author: Andrea Bassi @Polimi
"""

import napari
from napari.layers import Image
import numpy as np
from magicgui.widgets import FunctionGui
from magicgui import magic_factory
from enum import Enum
from functools import partial
# from napari.qt.threading import thread_worker

def calculate_widefield(**kwargs):
    psf_det = kwargs['psf_det']
    return 1*psf_det

def calculate_confocal(**kwargs):
    psf_ill = kwargs['psf_ill']
    psf_det = kwargs['psf_det']
    return psf_ill*psf_det

def calculate_2p(**kwargs):
    psf_ill = kwargs['psf_ill']
    return psf_ill**2

def calculate_sted(**kwargs):
    psf_ill = kwargs['psf_ill']
    psf_sted = kwargs['psf_sted']
    sat_ratio = kwargs['sat_ratio']
    I_sted_max = np.amax(psf_sted)
    result = psf_ill * np.exp(-sat_ratio * psf_sted / I_sted_max) # TODO check normalization
    return result

class Microscope(Enum):
    """
    Enum with microscopes types used for selecting microscope in the combobox of the napari widget
    Each microscope type is associated to a function (e.g. partial(calculate_widefield)) .
    The boolean values indicate which ui boxes to show in the napari widget (see _on_microscope_changed function).
    """
    widefield = ( partial(calculate_widefield), False,True,False,False ) 
    confocal = ( partial(calculate_confocal), True,False,True,True )
    two_photons = ( partial(calculate_2p), True,False,False,False )
    sted = ( partial(calculate_sted), True,False,True,True )
    def __new__(cls, *values):
        obj = object.__new__(cls)
        obj.function = values[0]
        obj.illumination_visible =  values[1]
        obj.detection_visible =  values[2]
        obj.sted_visible =  values[3]
        obj.saturation_visible =  values[4]
        return obj


def calc_init(calc_widget: FunctionGui):
    @calc_widget.microscope.changed.connect
    def _on_microscope_changed(microscope: Microscope):
        calc_widget.illumination_psf.visible = microscope.illumination_visible
        calc_widget.detection_psf.visible = microscope.detection_visible
        calc_widget.sted_psf.visible = microscope.sted_visible
        calc_widget.saturation_ratio.visible = microscope.saturation_visible

        
@magic_factory(widget_init=calc_init,
               call_button="Calculate",
               sted_psf={"visible":False},
               saturation_ratio={"visible":False, "min":1e-9, "max":1e9, "step":0.1}) 
def calculate(viewer: napari.Viewer,
            illumination_psf: Image,
            detection_psf: Image,
            sted_psf: Image,
            saturation_ratio: float = 10,
            microscope: Microscope = Microscope.confocal,
            ):
    '''
    Calculates the PSF of  microscopes starting from the illumination and detection PSF.
    Performs operation between the PSF:
    - widefield PSF is equal to the detection PSF 
    - confocal PSF is calculated as product of illumination and detection
    - 2-photons is calculated as squared value of the illumination
    - STED is calculated according to:
        Volker Westphal and Hell, Stefan W. 
        "Nanoscale Resolution in the Focal Plane of an Optical Microscope.” 
        Physical Review Letters 98, 143903, (2005)
    '''
    if illumination_psf is not None:
        psf_function = microscope.function
        result = psf_function(psf_ill = illumination_psf.data,
                                psf_det = detection_psf.data,
                                psf_sted = sted_psf.data,
                                sat_ratio = saturation_ratio)

        viewer.add_image(result, 
                        scale = illumination_psf.scale,
                        colormap = 'twilight',
                        name = microscope.name)
    
    
    