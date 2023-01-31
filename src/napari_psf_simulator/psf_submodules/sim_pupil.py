# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 10:20:28 2020

@author: Andrea Bassi
"""
import numpy as np

def gaussian_kernel(X,Y, Wx, Wy, X0=0, Y0=0):
    """ creates a 3D gaussian kernel 
    numpy array with size Ly,Lx
    waist along x, y
    traslated from origin by X0,Y0
    """
    
    kern = 1.0 * np.exp( - ((X-X0)**2)/Wx**2 
                         - ((Y-Y0)**2)/Wy**2
                        )
    return kern


def multiple_gaussians(x, y, waistx, waisty, rhos, amps, thetas, normalized = True):

    #disc = np.zeros((x.shape[0],x.shape[1]), dtype='complex128')
    disc = np.zeros_like(x)
    
    for rho,theta,amp in zip(rhos,thetas,amps):
        disc += amp*gaussian_kernel(x, y, waistx, waisty,
                                rho*np.cos(theta), rho*np.sin(theta))
    
    r = np.sqrt(x**2+y**2) 
    disc = disc * (r <= 1.) 
    if normalized:
        disc = disc / np.sum(disc)
    
    return disc
        

if __name__ == "__main__": 
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    Npixels = 512 # number of pixels
    R = 1 # extent of the xy space
    x = y = np.linspace(-R, +R, Npixels)
    X, Y = np.meshgrid(x,y)
    
    rho = np.sqrt(X**2 + Y**2)
    theta = np.arctan2(Y,X)
    
    n = 5  # Zernike radial order 
    m = 1 # Zernike azimutal frequency
    
    Z = multiple_gaussians(X,Y, 0.1, 0.2, [0.5,1,0.5], [1,1,1], [0, 2*np.pi/3, 4*np.pi/3], normalized = True) 
    
    fig1 = plt.figure(figsize=(9, 9))
    
    plt.title(f'Zernike polynomial of order ({n},{m})')
    
    plt.imshow(np.abs(Z), 
               interpolation='none',
               cmap=cm.gray,
               origin='lower',
               extent = [-R,R,-R,R]
               )
    
    plt.colorbar()
   