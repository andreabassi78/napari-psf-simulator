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


def set_contrast_init(contrast_widget: FunctionGui):
    @contrast_widget.reference_layer.changed.connect
    def _on_reference_layer_changed(reference_layer:Image):
        if reference_layer is not None:
            layer_data = reference_layer.data
            contrast_widget.min.value = np.amin(layer_data)
            contrast_widget.max.value = np.amax(layer_data)
            

@magic_factory(widget_init = set_contrast_init,
               call_button = "Set contrast",
               min={"min":0.0, "max":np.Inf},
               max={"min":0.0, "max":np.Inf}
               )
def set_contrast(viewer: napari.Viewer,
                reference_layer: Image,
                min: float = 0.0,
                max: float = 1.0,
                propagate_to_all_images: bool = False,
                ):
    
    if reference_layer is not None:
        #reference_layer.contrast_limits_range = [min,max]
        reference_layer.contrast_limits = [min,max]
        viewer.layers.selection.active = reference_layer
        if propagate_to_all_images:
            for layer in viewer.layers:
                if type(layer) is Image:
                    #layer.contrast_limits_range = [min,max]
                    layer.contrast_limits = [min,max]
        
