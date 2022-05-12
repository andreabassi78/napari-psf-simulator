# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 10:20:28 2020

Defines a function to create Zernike Polynomials.
Based on the code developed by Martin Weigert: 
https://github.com/mpicbg-csbd/phasenet

@author: Andrea Bassi
"""
import numpy as np
from scipy.special import binom

def nm_normalization(n, m):
    """the norm of the zernike mode n,m in born/wolf convetion
    i.e. sqrt( \int | z_nm |^2 )
    """
    return np.sqrt((1.+(m==0))/(2.*n+2))

def nm_polynomial(n, m, rho, theta, normalized=True):
 
    """returns the zernike polyonimal by classical n,m enumeration
    if normed=True, then they form an orthonormal system
        \int z_nm z_n'm' = delta_nn' delta_mm'
        and the first modes are
        z_nm(0,0)  = 1/sqrt(pi)*
        z_nm(1,-1) = 1/sqrt(pi)* 2r cos(phi)
        z_nm(1,1)  = 1/sqrt(pi)* 2r sin(phi)
        z_nm(2,0)  = 1/sqrt(pi)* sqrt(3)(2 r^2 - 1)
        ...
        z_nm(4,0)  = 1/sqrt(pi)* sqrt(5)(6 r^4 - 6 r^2 +1)
        ...
    if normed =False, then they follow the Born/Wolf convention
        (i.e. min/max is always -1/1)
        \int z_nm z_n'm' = (1.+(m==0))/(2*n+2) delta_nn' delta_mm'
        z_nm(0,0)  = 1
        z_nm(1,-1) = r cos(phi)
        z_nm(1,1)  =  r sin(phi)
        z_nm(2,0)  = (2 r^2 - 1)
        ...
        z_nm(4,0)  = (6 r^4 - 6 r^2 +1)
        
    """
    if abs(m) > n:
        raise ValueError(" |m| <= n ! ( %s <= %s)" % (m, n))

    if (n - m) % 2 == 1:
        return 0 * rho + 0 * theta

    radial = 0
    m0 = abs(m)

    for k in range((n - m0) // 2 + 1):
        radial += (-1.) ** k * binom(n - k, k) * binom(n - 2 * k, (n - m0) // 2 - k) * rho ** (n - 2 * k)

    radial = radial * (rho <= 1.) 

    if normalized:
        prefac = 1. / nm_normalization(n, m)
    else:
        prefac = 0.5
    if m >= 0:
        return prefac * radial * np.cos(m0 * theta)
    else:
        return prefac * radial * np.sin(m0 * theta)
        

if __name__ == "__main__": 
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    Npixels = 512 # number of pixels
    R = 1 # extent of the xy space
    x = y = np.linspace(-R, +R, Npixels)
    X, Y = np.meshgrid(x,y)
    
    rho = np.sqrt(X**2 + Y**2)
    theta = np.arctan2(Y,X)
    
    n = 4  # Zernike radial order 
    m = -4 # Zernike azimutal frequency
    
    Z = nm_polynomial(n, m, rho, theta, normalized = False) 
    
    fig1 = plt.figure(figsize=(9, 9))
    
    plt.title(f'Zernike polynomial of order ({n},{m})')
    
    plt.imshow(Z, 
               interpolation='none',
               cmap=cm.gray,
               origin='lower',
               extent = [-R,R,-R,R]
               )
    
    plt.colorbar()
    
    print('Peak to valley distance:', np.amax(Z)-np.amin(Z))
   