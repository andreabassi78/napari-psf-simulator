# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 17:02:23 2021

@author: Andrea Bassi
"""

import numpy as np
from numpy.fft import fft2, ifftshift, fftshift, fftfreq
from warnings import warn
from .zernike_polynomials import nm_normalization, nm_polynomial 

class PSF_simulator():
    '''
    Class to generate 3D Point Spread Functions
    with different pupils and various abberrations. 
    '''
    
    def __init__(self, NA=0.5, n=1, wavelength=0.532, Nxy=127, Nz=3, **kwargs):
        '''
        NA: numerical aperture
        n: refractive index
        wavelength
        Nxy: number of pixels in kx-ky (pixels in x-y is Nxy-1)
        Nz: number of pixels in kz (and in z) 
        over_sampling: ratio between the Abbe resolution and spatial sampling
        aspect_ratio: ratio between z and xy sampling
        
        '''
        if Nxy % 2 == 0:
            Nxy +=1 
            warn(f'Number of pixels Nxy must be odd, changed to {Nxy}')
        if Nz % 2 == 0:
             Nz +=1 #XY number of pixels must be odd 
             warn(f'Number of pixels Nz must be odd, changed to {Nz}')
             
        self.NA = NA # Numerical aperture
        self.n = n # refraction index at the object
        self.wavelength = wavelength
        
        self.Nxy = Nxy
        self.Nz = Nz
        
        DeltaXY = wavelength/2/NA # Diffraction limited transverse resolution
        
        if not 'dr' in kwargs:
            over_sampling = 4 # spatial sampling in xy, chosen to be a fraction (over_sampling) of the resolution
            self.dr = DeltaXY/over_sampling 
        else: 
            self.dr = kwargs['dr']
        
        if not 'dz' in kwargs:
            aspect_ratio = 4
            self.dz = aspect_ratio * self.dr # spatial sampling in z
        else: 
            self.dz = kwargs['dz']

        self.generate_kspace()
        
    def re_init(self,*args,**kwargs):
        """
        reinitializes the attributes of the PSFsimulator object wihout creating a new instance (__new__ is not executed)
        """
        self.__init__(*args,**kwargs)
        
    def generate_kspace(self):
        """
        Generates the k-space used to define the transfer functions
        """
        # generate the real  space (xy is used only for visulaization) 
        # x = y = fftshift(fftfreq(Npixels, dk))
        self.x = self.y = self.dr * (np.arange(self.Nxy) - self.Nxy // 2)
        self.z = self.dz * (np.arange(self.Nz) - self.Nz // 2)
                
        self.k = self.n/self.wavelength # wavenumber
        self.k_cut_off = self.NA/self.wavelength # cut off frequency in the coherent case
        
        
        kx_lin = fftshift(fftfreq(self.Nxy, self.dr))
        ky_lin = fftshift(fftfreq(self.Nxy, self.dr))
        # kx_lin = dk * (np.arange(Nxy) - Nxy // 2)
        # ky_lin = dk * (np.arange(Nxy) - Nxy // 2)
        self.dk = kx_lin[1]-kx_lin[0]
        kx, ky = np.meshgrid(kx_lin,ky_lin) 

        # k-space in radial coordinates
        with np.errstate(invalid='ignore'):    
            self.k_rho = np.sqrt(kx**2 + ky**2)
            self.k_theta = np.arctan2(ky,kx)  
            self.kz = np.sqrt(self.k**2-self.k_rho**2)
        
        self.kx = kx
        self.ky = ky
        
        self.phase = np.zeros_like(self.kz)
        self.amplitude = np.ones_like(self.kz)
        
        self.ATF0 = np.zeros_like(self.kz) 
        self.PSF3D = np.zeros(((self.Nz,self.Nxy,self.Nxy)))

    def add_slab_scalar(self, n1, thickness, alpha):
        """ 
        Calculates the effect of the slab, using scalar theory, without considering s and p polarizations.
        n1: refractive index of the slab
        thickness: thickness
        alpha: angle from the xy plane, conventionally from the y axis 
        """
        
        n0 = self.n
        k= self.k
        k1 = n1/n0 * k
        
        ky=self.ky
        kx=self.kx
        kz=self.kz
            
        # note that, at alpha=0, k_rho remains the same in the slab, in agreement with Snell's law
        
        with np.errstate(invalid='ignore'):
            theta0 = np.arcsin(ky/k)
            theta1 = np.arcsin(n0/n1 * np.sin(theta0 + alpha)) - alpha
            ky1 = k1 * np.sin (theta1)
            k_rho1 = np.sqrt( kx**2 + ky1**2 )
            kz1 = np.sqrt( k1**2 - k_rho1**2 ) 
    
        # additional phase due to propagation in the slab
        phase = 2*np.pi * (kz1-kz) * thickness / np.cos(alpha) 
        
        # Fresnel law of refraction (not used and not important al low NA)
        # Ts01 = 2 * n0 * np.cos(theta0) / (n0* np.cos(theta0) + n1 * np.cos(theta1))
        # Tp01 = 2 * n0 * np.cos(theta0) / (n0* np.cos(theta1) + n1 * np.cos(theta0))
        # Ts10 = 2 * n1 * np.cos(theta1) / (n1* np.cos(theta1) + n0 * np.cos(theta0))
        # Tp10 = 2 * n1 * np.cos(theta1) / (n1* np.cos(theta0) + n0 * np.cos(theta1))
        # transmittance = (Ts01*Ts10 + Tp01*Tp10 ) / 2 # assuming equal s and p polarization components
        # self.amplitude *= transmittance
        
        self.phase += phase
        self.thickness = thickness
        self.alpha = alpha
        self.n1 = n1
        self.correct_slab_defocus()
        
    def correct_slab_defocus(self):    
        """ Calculates the displacement of the focus along z and y 
        as it was calculated by ray-tracing and changes the phase accordingly.
        Correcting for defocus does not change the PSF shape but recenters it around the origin
        """
        n0 = self.n
        NA= self.NA
        
        if hasattr(self, 'thickness'):
            thickness = self.thickness
            n1 = self.n1
            alpha = self.alpha
        else: 
            raise Exception('Slab parameters not specified')
            
        maxtheta0 = np.arcsin(NA/n0)
        maxtheta1 = np.arcsin(NA/n1)
        
        # diplacement along z (it is calculated at alpha==0)
        self.displacementZ = thickness * (1-np.tan(maxtheta1)/np.tan(maxtheta0))
        
        # calculate the displacement of a paraxial ray from the optical axis 
        # (this is correctly calculated as a fucntion of alpha)
        alpha1 = np.arcsin(n0/n1*np.sin(alpha))
        self.displacementY = - thickness/ np.cos(alpha1)*np.sin(alpha-alpha1)
        
        # correct for defocus, to recenter the PSF in z==0 and y==0 (the new ray-tracing focus point)
        self.phase +=  2*np.pi * self.kz * self.displacementZ
        self.phase +=  2*np.pi * self.ky * self.displacementY
        
        # remove piston
        phase = self.phase
        phase = phase[np.isfinite(phase)]
        self.phase = (self.phase - np.min(phase)) 
        
        
    def add_cylindrical_lens(self, f_cyl, f):
        """introduce astigmatism placing a thin lens in the pupil
        which produces a quadratic phase in direction x.
        f_cyl is the focal length of the cylindrical lens
        f is the focal length of the objective lens 
        """
        self.phase += 2*np.pi * (self.kx)**2 /2/f_cyl * f
        
        
    def add_Zernike_aberration(self, N, M, weight):
        # weight of the polynomials in units of lambda (weight 1 means  wavefront abberated of lamba/2)
        self.N = N
        self.M = M
        self.weight = weight
        self.phase += 2* np.pi* weight*nm_polynomial(N, M, 
                                           self.k_rho/self.k_cut_off, 
                                           self.k_theta, normalized = False
                                           ) 
        
    def generate_pupil(self):
        """
        generate the pupil and the transfer functions 
        """
        ATF0 = self.amplitude * np.exp( 1.j*self.phase) 
        #ATF0 = np.exp( 1.j*phase) # Amplitude Transfer Function (pupil)
        cut_idx = (self.k_rho >= self.k_cut_off) # indexes of the evanescent waves (kz is NaN for these indexes)
        ATF0[cut_idx] = 0 # exclude k above the cut off frequency
        self.ATF0 = ATF0
        
    def add_Ndimensional_SIM_pupil(self,     
                     kr = 0.5,
                     waist = 0.01,
                     source_num = 3,
                     add_cw = False
                     ):
        ''' Generates the pupil for mutlidimensional SIM microscopy
            (typically 2,3 or 4 sources are used)
            kr: spatial frequency (radial component) of the sources, 
                relative to the  cutoff frequency self.k_cut_off
            waist: of the gaussian source,
                relative to the  cutoff frequency self.k_cut_off
            
            '''
        from .sim_pupil import multiple_gaussians

        source_theta = list(2*np.pi/source_num * np.arange(source_num))
        source_kr = [kr] * source_num # list with NumSources elements
        source_amp = [1] * source_num # list with NumSources elements
        
        if add_cw:
            source_kr.append(0.)
            source_theta.append(0.)
            source_amp.append(2)#sum(source_amp))
        
        self.amplitude *= multiple_gaussians(self.kx/self.k_cut_off,
                                  self.ky/self.k_cut_off,
                                  waist, waist,
                                  source_kr, source_amp, source_theta)
        
        
    
    def add_lightsheet_pupil(self, 
                     waistx = 0.015,
                     waist_ratio = 20,
                     ):
        ''' Generates the pupil for light sheet illumination
        waistx: of the gaussian source along. 
        waist_ratio: ratio between the waist along y and the one along x
        If waist_x<<1 waist_ratio>>1, a light sheet is formed in the plane xz 
        
        '''
        from .sim_pupil import multiple_gaussians
        waisty = waistx * waist_ratio
        beam= multiple_gaussians(self.kx/self.k_cut_off,
                                  self.ky/self.k_cut_off,
                                  waistx, waisty,
                                  [0.0], [0.0])  
        
        self.amplitude *=beam
        
        
    def add_conical_pupil(self, thickness , n1 = 1.51):
        ''' Simulates the presence of an axicon in the pupil 
        alpha: angle of the cone relative to the xy plane
        '''
        
        delta = 2*np.pi * thickness * (n1 - self.n) / self.wavelength
        ax_phase = delta * (1 - self.k_rho / self.k_cut_off) 
        self.phase += ax_phase 
        
        
    def add_lattice_pupil(self, 
                     cutin = 0.84,
                     cutout = 1,
                     waistx = 0.015,
                     waist_ratio = 20,
                     source_num = 3
                     ):
        ''' Generates the pupil for lattice light sheet microscopy
            All parameters are relative to the cutoff frequency self.k_cut_off
        cutin: minimum radial spatial frequency of the annular ring
        cutout: maximum radial spatial frequency of the annular ring
        waistx: of the gaussian source along x
        waist_ratio: ratio between the waist along y and the one along x
        source_num: order of the lattice
        '''
        from .sim_pupil import multiple_gaussians
        source_rho = [(cutout+cutin)/2] * source_num # repeat list source_num times
        source_theta = 2*np.pi/source_num * np.arange(source_num)
        
        waisty = waistx * waist_ratio
        
        beams= multiple_gaussians(self.kx/self.k_cut_off,
                                  self.ky/self.k_cut_off,
                                  waistx, waisty,
                                  source_rho, source_theta)  
        
        cut_idx = (self.k_rho <= self.k_cut_off*cutin) | (self.k_rho >= self.k_cut_off*cutout)  
        mask = np.ones_like(self.k_rho)
        mask[cut_idx] = 0 # exclude k above the cut off frequency
        
        self.amplitude *=beams*mask
                
    def generate_3D_PSF(self):
        ''' core function of the class.
        Generates the PSF by Fourier Transforming the Amplitude Tranfer Function
        calculated at different z with the Rayleigh-Sommerfield angular spectrum 
        '''
        
        ATF0 = self.ATF0
        
        for idx,zi in enumerate(self.z):
           
            angular_spectrum_propagator = np.exp(1.j*2*np.pi*self.kz*zi)
             
            ATF = ATF0 * angular_spectrum_propagator
        
            evanescent_idx = (self.k_rho > self.k)
            ATF[evanescent_idx] = 0 
            
            ASF = fftshift(fft2(ifftshift(ATF))) #* k**2/f**2 # Amplitude Spread Function
            
            PSF = np.abs(ASF)**2 # Point Spread Function
            
            self.PSF3D[idx,:,:] = PSF
       
    def calculateRMS(self):
        '''calculates the RMS wavefron error
        For Zernike abberrations it is the weight of the 
        rms == 1 indicates 1-wavelength mean wavefront error 
        '''
        cut_idx = self.k_rho >= self.k_cut_off  
        area = self.ATF0.size - np.sum(cut_idx)
        phase = self.phase # /(2*np.pi) #TODO check notation ofr RMS calculation
        #m = np.zeros_like(phase)
        phase[cut_idx] = 0
        phase = phase[np.isfinite(phase)]
        mean_phase = np.mean(phase)
        rms = np.sqrt(np.sum((phase-mean_phase)**2)/area) 
        return(rms)
       
    def print_values(self):
    
        DeltaXY = self.wavelength/2/self.NA # Diffraction limited transverse resolution
        DeltaZ = self.wavelength/self.n/(1-np.sqrt(1-self.NA**2/self.n**2)) # Diffraction limited axial resolution
    
        print(f'The numerical aperture of the system is: {self.NA}') 
        print(f'The transverse resolution is: {DeltaXY:.03f} um') 
        print(f'The axial resolution is: {DeltaZ:.03f} um') 
        print(f'The pixel size is: {self.dr:.03f} um') 
        print(f'The voxel depth is: {self.dz:.03f} um') 
        if hasattr(self, 'displacementZ'):
            print(f'The displacement z from the focus is: {self.displacementZ} um')
        if hasattr(self, 'displacementY'):
            print(f'The displacement y from the optical axis is: {self.displacementY} um')
        if hasattr(self, 'calculateRMS'):
            print(f'The RMS wavefront error is {self.calculateRMS()}')    
            
    def show_pupil(self):
        """ 
        Shows the Amplitude Transfer Function (pupil) in amplitude and phase
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(12, 6), tight_layout=False, dpi=300)
        sup_title =  f'NA = {self.NA}, n = {self.n}'
        if hasattr(self,'thickness'):
             sup_title += f', slab thickness = {self.thickness} $\mu$m, n1 = {self.n1}, alpha = {self.alpha:.02f}'
        fig.suptitle(sup_title)     
        
        for idx in range(2):
            if idx==0:
                data = np.abs(self.ATF0)
                title = 'Pupil (amplitude)'
            elif idx ==1:
                data = np.angle(self.ATF0)
                title = 'Pupil (phase)'
                    
        
            im0=ax[idx].imshow(data, 
                             cmap='pink',
                             extent = [np.amin(self.kx),np.amax(self.kx),np.amin(self.ky),np.amax(self.ky)],
                             origin = 'lower',
                             interpolation = 'none'
                             )
            ax[idx].set_xlabel('kx (1/$\mu$m)')
            ax[idx].set_ylabel('ky (1/$\mu$m)')
            ax[idx].set_title(title)
            fig.colorbar(im0,ax = ax[idx])       
       
    def show_PSF_projections(self, aspect_ratio, mode ='MIP', dpi=150):
        """ Shows the 3D PSF in 3 orthogonal views
        Aspect ratio: between the z axis and the other axes scale
        mode: 'MIP' or 'plane': MAximum intensity projection or plane incercepting the origin
        """
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 3, figsize=(12,6), tight_layout=False, dpi=dpi)
        sup_title =  f'NA = {self.NA}, n = {self.n}'
        if hasattr(self,'thickness'):
             sup_title += f', slab thickness = {self.thickness} $\mu$m, n1 = {self.n1}, alpha = {self.alpha:.02f}'
        fig.suptitle(sup_title)
        
        label_list = ( ('x','y'),('x','z'),('y','z') )
        
        for idx, labels in enumerate(label_list):

            if mode =='mip':
                # create maximum intensity projection
                MIP = np.amax(self.PSF3D, axis=idx)
                im_to_show = MIP
            
            elif mode =='sum':
                # create average intensity projection
                MIP = np.mean(self.PSF3D, axis=idx)
                im_to_show = MIP    
                
            elif mode =='plane':
                PSF = self.PSF3D
                Nz,Ny,Nx = PSF.shape
                Nlist = [Nz,Ny,Nx]
                im_to_show = PSF.take(indices=Nlist[idx]//2 , axis=idx)
            else:
                raise(ValueError, 'Please specify PSF showing mode' )
 
            values0 = getattr(self, labels[0])
            values1 = getattr(self, labels[1])
    
            extent = [np.amin(values0), np.amax(values0),
                      np.amin(values1), np.amax(values1)]
            
            if idx == 0:
                vmin = np.amin(im_to_show)
                vmax = np.amax(im_to_show)
                
            
            axs[idx].imshow(im_to_show,
                             cmap='twilight',
                             extent = extent,
                             origin = 'lower',
                             vmin=vmin, vmax=vmax
                             )
            
            axs[idx].set_xlabel(f'{labels[0]} ($\mu$m)')
            axs[idx].set_ylabel(f'{labels[1]} ($\mu$m)')
            axs[idx].set_title(f'|PSF({labels[0]},{labels[1]})|')  
            
            
            if labels[1] == 'z': 
                axs[idx].set_aspect(1/aspect_ratio)
        
    def plot_phase(self, dpi =150):
        '''
        plot the phase a the pupil along x and y
        '''
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(12, 6), tight_layout=False, dpi=dpi)
        sup_title =  f'NA = {self.NA}, n = {self.n}'
        if hasattr(self,'thickness'):
             sup_title += f', slab thickness = {self.thickness} $\mu$m, n1 = {self.n1}, alpha = {self.alpha:.02f}'
        fig.suptitle(sup_title)
        phase = self.phase
        s = phase.shape
        ax[0].plot( self.kx[s[0]//2,:], phase[s[0]//2,:] )
        ax[0].set_xlabel('kx (1/$\mu$m)')
        ax[0].set_ylabel('phase (rad)')
        
        ax[1].plot( self.ky[:,s[1]//2], phase[:,s[1]//2] )
        ax[1].set_xlabel('ky (1/$\mu$m)')
        ax[1].set_ylabel('phase (rad)')
        
        
    def plot_psf_profile(self, dpi = 150):
        '''
        plot the phase a the pupil along x and y
        '''
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(8, 4), tight_layout=False, dpi=dpi)
        char_size = 12
        sup_title =  f'NA = {self.NA}, n = {self.n}'
        if hasattr(self,'thickness'):
             sup_title += f', slab thickness = {self.thickness} $\mu$m, n1 = {self.n1}, alpha = {self.alpha:.02f}'
        fig.suptitle(sup_title, size=char_size*0.8)
        PSF = self.PSF3D
        Nz, Ny, Nx = PSF.shape
    
        # psf_to_show = PSF.take(indices=Nlist[idx]//2 , axis=idx)
        psf_to_show_x = PSF[Nz//2,Ny//2,:]
        
        psf_to_show_z = PSF[:,Ny//2,Nx//2]
        
        ax[0].plot( self.x, psf_to_show_x,
                   linewidth=1.5)
        ax[0].set_xlabel('x ($\mu$m)',size=char_size)
        ax[0].set_ylabel('PSF', size=char_size)
        ax[0].grid()
        DeltaX = self.wavelength/self.NA/2 # Abbe resolution
        ax[0].plot(np.array([0.,DeltaX,1.22*DeltaX]),
                   np.array([0.,0.,0.]),
                   'o', markersize=2)
                  
        
        ax[1].plot( self.z, psf_to_show_z,
                   linewidth=1.5)
        ax[1].set_xlabel('z ($\mu$m)', size=char_size)
        # ax[1].set_ylabel('PSF')    
        ax[1].grid()
        DeltaZ = self.wavelength/self.n/(1-np.sqrt(1-self.NA**2/self.n**2)) # Diffraction limited axial resolution
        ax[1].plot(DeltaZ, 0., 'o', markersize=2)
        
        for idx in (0,1):
            ax[idx].xaxis.set_tick_params(labelsize=char_size*0.5)
            ax[idx].yaxis.set_tick_params(labelsize=char_size*0.5)
        
        plt.show()
        
    def save_data(self):
        from skimage.external import tifffile as tif
        basename = 'psf'
        filename = '_'.join([basename,
                             f'NA_{self.NA}',
                             f'n_{self.n}'])
        
        if hasattr(self, 'thickness'):
            filename = '_'.join([filename,
                                 f'size_{self.thickness}',
                                 f'alpha_{self.alpha:.2f}',
                                 f'n1_{self.n1}'])
        
        
        
        psf16 = ( self.PSF3D * (2**16-1) / np.amax(self.PSF3D) ).astype('uint16') #normalize and convert to 16 bit
        psf16.shape = 1, self.Nz, 1, self.Nxy, self.Nxy, 1 # dimensions in TZCYXS order
        
        tif.imsave(filename+'.tif', psf16, imagej=True, resolution = (1.0/self.dr, 1.0/self.dr),
                    metadata={'spacing': self.dz, 'unit': 'um'})
        
    def write_name(self, basename =''):
        
        name = '_'.join([basename,
                        f'NA_{self.NA:.1f}',
                        f'n_{self.n:.1f}'])
       
        if all(hasattr(self, attr) for attr in ["thickness","alpha","n1"]): # slab abberation is there
            name = '_'.join([name,
                            f'thk_{self.thickness:.0f}',
                            f'alpha_{self.alpha:.0f}',
                            f'n1_{self.n1:.1f}'])
        
        if all(hasattr(self, attr) for attr in ["N","M","weight"]): # zernike aberration is there
            name = '_'.join([name,
                            f'N{self.N}',
                            f'M_{self.M}',
                            f'w_{self.weight:.1f}'])
        return name


if __name__ == '__main__':

    um = 1.0
    mm = 1000 * um
    deg = np.pi/180
    
    NA = 1.15
    
    wavelength = 0.532 * um 
    n0 = 1.33 # refractive index of the medium
    
    n1 = 1.47 # refractive index of the slab
    thickness = 50 * um # slab thickness
    alpha = 1 * deg # angle of the slab relative to the y axis
    
    SaveData = False
    
    Nxy = 127
    
    Nz = 63
     
    dr = 0.1 * um
    
    dz = 0.2 * um
    
    gen = PSF_simulator(NA, n0, wavelength, Nxy, Nz, 
                        dr = dr, dz = dz
                        )
    
    # gen.add_Ndimensional_SIM_pupil(kr = 0.5,
    #                                 waist = 0.02,
    #                                 source_num = 3,
    #                                 add_cw = True)
    
    # gen.add_lattice_pupil()
    # gen.add_lightsheet_pupil()
    
    # gen.add_slab_scalar(n1, thickness, alpha)
    
    # gen.add_slab_vectorial(n1, thickness, alpha)
    # gen.add_Zernike_aberration(3, 1, weight=1)
    # gen.add_conical_pupil(0.06*mm)
    
    # gen.add_cylindrical_lens(f_cyl=-40, f=20)
    
    # rms = gen.calculateRMS()
    # print(rms)
    
    
    gen.generate_pupil()
    gen.generate_3D_PSF()
    
    # Show results    
    gen.print_values()
    
    gen.show_pupil()
    #gen.plot_phase()
    gen.show_PSF_projections(aspect_ratio=1, mode='plane')
    gen.plot_psf_profile()
       
