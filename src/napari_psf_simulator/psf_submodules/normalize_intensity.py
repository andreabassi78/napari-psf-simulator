"""This script contains helper functions to normalize the PSF intensity. 
In order to achieve this, conservation of energy from an incident uniform beam is considered.
Then, the integral of the intensity near the focus must be the same as the one inciding on the focal lens.

To calculate the integral of the intensity near the focus, a 2D numerical integration must be used. 
At the moment, the method employed for this is the trapezoidal approximation.
"""

import numpy as np

def calculate_2D_trapezoidal_method_weight(x_max, y_max, Nx, Ny):
    """Calculates the weights for the trapezoidal 2D integration"""
    h_x=x_max/Nx
    h_y=y_max/Ny
    weight_trapezoid_x=np.zeros(Nx)+h_x
    weight_trapezoid_x[0]=h_x/2
    weight_trapezoid_x[-1]=h_x/2
    weight_trapezoid_y=np.zeros(Ny)+h_y
    weight_trapezoid_y[0]=h_y/2
    weight_trapezoid_y[-1]=h_y/2
    return weight_trapezoid_x*np.vstack(weight_trapezoid_y)#represents the area of each trapezoid for each position in x,y

def calculate_normalizing_factor(intensity_at_the_focus, x_max: float, y_max: float, lens_aperture: float):
    """Calculates the factor needed to normalize the intensity of the field near the focus. So that multiplying by it normalizes the field
     The calculus is based on the conservation of energy since the integral of the intensity near the focus must be the same as the one inciding on the focal lens.

    To calculate the integral of the intensity near the focus, a 2D numerical integration must be used. 
    At the moment, the method employed for this is the trapezoidal approximation.


    Args:
        intensity_at_the_focus (2D array): Array containing the values of the intensity of the field near the focus for each pair of x, y coordinates
        x_max (float): Maximum value of X for wich the field was calculated near the focus (nm)
        y_max (float): Maximum value of Y for wich the field was calculated near the focus (nm)
        lens_aperture (float): Radius of apperture of the focal lens (nm)

    Returns:
        float: Factor by wich to multiply the intensity to normalize it
    """
    E_inc = np.pi*(lens_aperture**2) # Energy inciding at the focal lens (equalt to the integral over a circle of radius "lens_aperture")
    E_inc/=10**12 # Unit transformation from mm^2 to nm^2
    Nz, Nx, Ny = np.shape(intensity_at_the_focus) # Number of division for each coordinate

    weight_trapezoid = calculate_2D_trapezoidal_method_weight(2*x_max, 2*y_max, Nx, Ny) # We multiply by 2 since the field is calculated from -x_max yo +x_max
    E_focus = np.abs(np.sum(intensity_at_the_focus*weight_trapezoid))
    
    return E_inc/E_focus

