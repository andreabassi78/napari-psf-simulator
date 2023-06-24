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
import enum
from enum import Enum
from scipy.ndimage import convolve
from napari.qt.threading import thread_worker


class Operation(enum.Enum):
    multiply = 0
    convolve = 1
    square = 2



@magic_factory(call_button="Calculate")
def calculate(viewer: napari.Viewer,
            first_layer: Image,
            second_layer: Image,
            operation: Operation,
            ):
    '''
    Performes operations between stacks
    Parameters
    ----------
    viewer : napari.Viewer
        Current viewer
    '''
    

    def update_layer(image):
        viewer.add_image(image, 
                        scale = first_layer.scale,
                        colormap = 'twilight',
                        name = operation.name)

    @thread_worker(connect={'returned': update_layer})
    def _calculate(first_data,second_data, operation):
            calculate.enabled = False
            if first_data is not None :
                if operation.name == 'square':
                    result = np.power(first_data, 2)
                elif operation.name == 'multiply':
                    if second_data is not None:        
                        result = np.multiply(first_data, second_data)
                elif operation.name == 'convolve':
                    if second_data is not None:
                        result = convolve(first_data, second_data)
            calculate.enabled = True
            return result
    
    _calculate(first_layer.data,second_layer.data,operation)
    
    
    