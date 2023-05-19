# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 00:16:58 2022

@author: Andrea Bassi @Polimi
"""

from enum import Enum
import numpy as np
from napari.qt.threading import thread_worker
from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QCheckBox,
    QLabel,
    QPushButton,
    QSplitter,
    QVBoxLayout,
    QWidget,
)
from qtpy.QtGui import QColor
from .psf_submodules.gui_utils import Setting, SwitchableSection
from .psf_submodules.psf_generator import PSF_simulator
from .PyFocus.src.napari_adapter.napari_adapter import PyFocusSimulator
# from PyFocus.napari_adapter import PyFocusSimulator  # TODO add reference to pyfocus after debuging

class Psf_widget(QWidget):
    '''
    Napari widget to simulate the point spread function of an imaging system,
    particularly an optical microscope.
    Aberrations currently included are:
        -Zernike phase aberration
        -presence of a layer with different refractive index than the sample 
    The class is constructed in a modular way so that additional aberration
    can be added to the UI (Combo_box) with new Settings containing the parameters
    needed to simulate the phase aberrtion.
    For adding aberrations change the self.setup_aberrations() method, 
    adding a new element to Aberrations as defined in napari_psf_simulator.aberrations  
    '''
    generators = Enum('Gen', 
                      {"Scalar": PSF_simulator,
                       "Vectorial": PyFocusSimulator
                        })
    
    pyfocus_amplitudes = Enum('py_amplitude', 
                      {"Uniform": 0,
                       "Gaussian": 1
                        })
    
    pyfocus_phases = Enum('py_phase', 
                      {"Uniform": 0,
                       "Vortex": 1
                        })
    
    aberrations = Enum('Aberrations',
                       {"NONE": 0,
                        "SLAB":1,
                        "ZERNIKE":2
                        })
    
    def __init__(self, napari_viewer,
                 ):
        self.viewer = napari_viewer
        super().__init__()
        
        self.setup_ui()
        
        
        
    def add_splitter(self, layout, title, color = (50,55,55)):
        splitter = QSplitter(Qt.Vertical)
        layout.addWidget(splitter)      
        lbl= QLabel(title)
        col = QColor(*color)
        lbl.setStyleSheet("QWidget { background-color: %s }" % col.name()) 
        layout.addWidget(lbl)
    
    
    def setup_ui(self): 
        """
        Creates the user interface. 
        Note that the magicgui is not used and the Setting class (gui_utils.Setting)  
        """
            
        # initialize layout
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # add Settings
        self.add_splitter(layout,'Settings', color=(2,10,10))
        settings_layout = QVBoxLayout()
        layout.addLayout(settings_layout)
        self.create_base_Settings(settings_layout)
        self.start_base_simulator()
        self.custom_mask = "1" #TODO: fix this and put in the correct place  
        
        # add plot results checkbox
        self.plot_checkbox = QCheckBox("Plot PSF profile in Console")
        self.plot_checkbox.setChecked(False)
        layout.addWidget(self.plot_checkbox)
        
        # add show Airy disk checkbox
        self.airy_checkbox = QCheckBox("Show Airy disk")
        self.airy_checkbox.setChecked(False)
        layout.addWidget(self.airy_checkbox)

        # add switchable section for Propagation model (scalar, vectorial)
        self.add_splitter(layout,'Propagation model', color=(2,10,10))
        self.generator_section = SwitchableSection(name = 'Model',
                                      baselayout = layout, choices = self.generators,
                                      on_change_function = self.change_generator)
        # add switchable section for aberrations
        self.add_splitter(layout,'Aberrations', color=(2,10,10))
        self.aberration_section = SwitchableSection(name = 'Aberrations',
                                      baselayout = layout, choices = self.aberrations,
                                      on_change_function = self.change_aberration)

        # add PSF calculation section
        self.add_splitter(layout,'Calculate Point Spread Function',color=(2,10,10))
        calculate_layout = QVBoxLayout()
        layout.addLayout(calculate_layout)
        calculate_btn = QPushButton('Calculate PSF')
        calculate_btn.clicked.connect(self.calculate_psf)
        layout.addWidget(calculate_btn)

    def change_generator(self):
        self.generator_section.remove_sub_layout_content()
        selection_name = self.generator_section.combo.text
        if selection_name == "Vectorial":
            self.add_PyFocus_settings(self.generator_section.sub_layout)
        self.reinitialize_simulator()


    def add_PyFocus_settings(self, layout):
        # add show intensity of each component checkbox

        self.add_splitter(layout, 'Pyfocus settings')
        self.pyFocus_component_checkbox = QCheckBox("show x,y,z intensities")
        self.pyFocus_component_checkbox.setChecked(False)
        layout.addWidget(self.pyFocus_component_checkbox)
        self.lens_aperture = Setting(name='lens radius', dtype=float, initial=3, unit='mm', 
                          layout = layout, write_function = self.reinitialize_simulator)
        
        self.amplitude_section = SwitchableSection(name = 'amplitude',
                                      baselayout = layout, choices = self.pyfocus_amplitudes,
                                      on_change_function = self.change_amplitude)
        self.phase_section = SwitchableSection(name = 'phase',
                                      baselayout = layout, choices = self.pyfocus_phases,
                                      on_change_function = self.change_phase)
        self.add_splitter(layout, 'Polarization')
        self.gamma = Setting(name='gamma', dtype=float, initial=45, unit='deg', 
                          layout = layout, write_function = self.reinitialize_simulator)
        self.beta = Setting(name='beta', dtype=float, initial=90, unit='deg',
                          layout = layout, write_function = self.reinitialize_simulator)

    def change_amplitude(self):
        self.amplitude_section.remove_sub_layout_content()
        if self.amplitude_section.combo.text == 'Uniform':
            self.custom_mask = "1" #TODO: fix this and put in the correct place   
        elif self.amplitude_section.combo.text == 'Gaussian':
            self.custom_mask = 'np.exp((rho/waist)**2)' #TODO: fix this (waist is not recongnized) and put in the correct place
            self.waist = Setting(name='waist', dtype=float, initial=2, unit='mm', 
                        layout = self.amplitude_section.sub_layout,
                        write_function = self.reinitialize_simulator)
        self.reinitialize_simulator()

    def change_phase(self):
        self.phase_section.remove_sub_layout_content()
        if self.phase_section.combo.text == 'Uniform':
            pass
            # self.custom_mask = '1' #TODO: fix this and put in the correct place
        elif self.phase_section.combo.text == 'Vortex':
            self.custom_mask = 'np.exp(1j*phi)' #TODO: insert the correct value
            self.order = Setting(name='order', dtype=int, initial=1, 
                        layout = self.phase_section.sub_layout,
                        write_function = self.reinitialize_simulator)
        self.reinitialize_simulator() 


    def change_aberration(self):
        self.aberration_section.remove_sub_layout_content()
        
        if self.aberration_section.combo.text == 'NONE':
            pass 
        
        elif self.aberration_section.combo.text == 'SLAB':           
            self.n1 = Setting(name='n1', dtype=float, initial=1.51, 
                        layout = self.aberration_section.sub_layout,
                        write_function = self.reinitialize_simulator) 
            self.thickness = Setting(name='thickness', dtype=float, initial=100.0, unit = '\u03BCm', 
                        layout = self.aberration_section.sub_layout,
                        write_function = self.reinitialize_simulator)
            self.alpha = Setting(name='alpha', dtype=float, initial=0.0, unit = 'deg', 
                        layout = self.aberration_section.sub_layout,
                        write_function = self.reinitialize_simulator)
        
        elif self.aberration_section.combo.text == 'ZERNIKE':           
            self.N = Setting(name='N', dtype=int, initial=3, 
                        layout = self.aberration_section.sub_layout,
                        write_function = self.reinitialize_simulator)
            self.M = Setting(name='M', dtype=int, initial=1, 
                        layout = self.aberration_section.sub_layout,
                        write_function = self.reinitialize_simulator)
            self.weight = Setting(name='weight', dtype=float, initial=0.6, unit = '\u03BB', 
                        layout = self.aberration_section.sub_layout,
                        write_function = self.reinitialize_simulator)
        self.reinitialize_simulator() 
 

    def create_base_Settings(self,settings_layout):
        """
        Adds the basic settings for psf-simulator with initial values
        
        Parameters
        ----------
        settings_layout : Qlayout
            layout where the Settings will appear. 
            Settings are Qwidgets with numerical or boolean values defined in gui_utils
        """
        self.NA = Setting(name='NA', dtype=float, initial=0.65, 
                          layout = settings_layout, write_function = self.reinitialize_simulator)
        self.n = Setting(name='n', dtype=float, initial=1.33, 
                          layout = settings_layout, write_function = self.reinitialize_simulator)
        self.wavelength = Setting(name='wavelength', dtype=float, initial=0.532,
                          unit = '\u03BCm', spinbox_decimals = 3, 
                          layout = settings_layout, write_function = self.reinitialize_simulator)
        self.Nxy = Setting(name='Nxy', dtype=int, initial=51, vmin=1, vmax = 4095,
                          layout = settings_layout, write_function = self.reinitialize_simulator)
        self.Nz = Setting(name='Nz', dtype=int, initial=3, vmin = 1, vmax = 4095,
                          layout = settings_layout, write_function = self.reinitialize_simulator)
        self.dxy = Setting(name='dxy', dtype=float, initial=0.03, unit = '\u03BCm', 
                          layout = settings_layout, write_function = self.reinitialize_simulator)
        self.dz = Setting(name='dz', dtype=float, initial=0.50, unit = '\u03BCm',
                          layout = settings_layout, write_function = self.reinitialize_simulator)
        
    
    def start_base_simulator(self):
        """Starts the base simulator, using scalar propagation
        """
        self.gen = PSF_simulator(self.NA.val, self.n.val, self.wavelength.val,
                        self.Nxy.val , self.Nz.val, dr = self.dxy.val, dz = self.dz.val)
    

    def reinitialize_simulator(self):
        '''
        Starts the simulator, accorgin to choice made on the combo box,
        Sets the current aberration type.
        '''
        selected_generator = self.generator_section.combo.current_data
        
        if selected_generator is PyFocusSimulator:
            self.gen = selected_generator(self.NA.val, self.n.val, self.wavelength.val,
                                self.Nxy.val , self.Nz.val, dr = self.dxy.val, dz = self.dz.val,
                                gamma = self.gamma.val, beta = self.beta.val,
                                custom_mask = self.custom_mask)
        elif selected_generator is PSF_simulator:
            self.gen = selected_generator(self.NA.val, self.n.val, self.wavelength.val,
                                self.Nxy.val , self.Nz.val, dr = self.dxy.val, dz = self.dz.val)
        self.add_aberration() 
    

    def add_aberration(self):
        '''
        Adds the effect of a phase aberration
        '''  
        if self.aberration_section.combo.text == 'NONE':
            pass
        elif self.aberration_section.combo.text == 'SLAB':
            self.gen.add_slab_scalar(self.n1.val, self.thickness.val, self.alpha.val)
        elif self.aberration_section.combo.text == 'ZERNIKE':
            self.gen.add_Zernike_aberration(self.N.val, self.M.val, self.weight.val)  

    
    def rescaleZ(self, layer):
        """
        Rescales the layer according to the dz/dxy settings values 
        """
        zscaling = self.dz.val / self.dxy.val
        if layer.ndim >=3:
            scale = layer.scale 
            scale[-3] = zscaling
            layer.scale = scale
            
    
    def calculate_psf(self):
        '''
        Calculates the 3D PSF and shows it as a stack
        '''
        @thread_worker
        def generator():
            import warnings
            warnings.filterwarnings('ignore')
            self.gen.generate_pupil()
            self.gen.generate_3D_PSF()
            return self.gen.PSF3D
        
        def update_image(psf):
            #self.viewer.dims.current_step = (0,0,0)
            
            psf_layer = self.viewer.add_image(psf,
                             name=self.gen.write_name(basename = 'Intensity'),
                             colormap='twilight')
            self.rescaleZ(psf_layer)
            
            if self.generator_section.combo.current_data == PyFocusSimulator and self.pyFocus_component_checkbox.checkState():
                x_layer = self.viewer.add_image(np.abs(self.gen.field.Ex)**2,
                                name=self.gen.write_name(basename = 'Ex'),
                                colormap='gray')
                self.rescaleZ(x_layer)
                y_layer = self.viewer.add_image(np.abs(self.gen.field.Ey)**2,
                                name=self.gen.write_name(basename = 'Ey'),
                                colormap='gray')
                self.rescaleZ(y_layer)
                z_layer = self.viewer.add_image(np.abs(self.gen.field.Ez)**2,
                                name=self.gen.write_name(basename = 'Ez'),
                                colormap='gray')
                self.rescaleZ(z_layer)
            
            if self.airy_checkbox.checkState():
                self.show_airy_disk()
            
            if self.plot_checkbox.checkState():
                self.gen.plot_psf_profile()
            
            posxy = self.Nxy.val // 2
            posz = self.Nz.val // 2
            self.viewer.dims.set_point(axis=[0,1,2], value=(0,0,0)) #raises ValueError in napari versions <0.4.13
            self.viewer.dims.current_step = (posz,posxy,posxy) # shows the image center of the stack in 3D
            
            
        worker = generator()  # create "worker" object
        worker.returned.connect(update_image)  # connect callback functions
        worker.start()  # start the thread
        
        
    def show_airy_disk(self):
       '''
       Shows the diffraction limited Airy disk as 3 perpendicular ellipses 
       '''
       deltaR = 1.22*self.wavelength.val/self.NA.val/2 # Rayleigh resolution
       deltaZ = self.wavelength.val/self.n.val/(1-np.sqrt(1-self.NA.val**2/self.n.val**2)) # Diffraction limited axial resolution
       #print(deltaZ)
       #print('paraxial:', 2*self.n.val* self.wavelength.val/ self.NA.val**2)
       
       posxy = self.Nxy.val//2
       posz = self.Nz.val//2
       center = np.array([posz,posxy,posxy])
       deltar = deltaR/self.dxy.val
       deltaz = deltaZ/self.dz.val
       # create three perpendicular ellipses
       bbox_yx = np.array([center+np.array([0, deltar, deltar]),
                           center+np.array([0, deltar,-deltar]),
                           center+np.array([0,-deltar,-deltar]),
                           center+np.array([0,-deltar, deltar])]
                          )
       bbox_zy = np.array([center+np.array([ deltaz, deltar,0]),
                           center+np.array([ deltaz,-deltar,0]),
                           center+np.array([-deltaz,-deltar,0]),
                           center+np.array([-deltaz, deltar,0])]
                          )
       bbox_zx = np.array([center+np.array([ deltaz, 0, deltar]),
                           center+np.array([ deltaz, 0,-deltar]),
                           center+np.array([-deltaz, 0,-deltar]),
                           center+np.array([-deltaz, 0, deltar])]
                          )
       
       if (self.dz.val*self.Nz.val)/2 > deltaZ:
           ellipses = [bbox_yx, bbox_zy, bbox_zx]
       else:
           ellipses = [bbox_yx]
       shapes_layer = self.viewer.add_shapes(name=self.gen.write_name(basename ='AiryDisk'),
                                              edge_width = 0.5,
                                              face_color = [1,1,1,0],
                                              edge_color = 'red')
       shapes_layer.add_ellipses(ellipses)  
       self.rescaleZ(shapes_layer)
       return(shapes_layer)
       
   
    def _show_PSF_projections(self): 
        '''
        Create 2D profiles as new layers.
        Not used
        '''
        PSF = self.gen.PSF3D
        if self.mip_checkbox.checkState():
            # create maximum intensity projection
            im_xy = np.amax(PSF, axis=0)
            im_xz = np.amax(PSF, axis=1)
            text = 'mip'
        else:
            Nz,Ny,Nx = PSF.shape
            im_xy = PSF[Nz//2,:,:]
            im_xz = PSF[:,Ny//2,:]
            text = 'plane'
            
        imageXZ = self.viewer.add_image(im_xz,
                     name= self.gen.write_name(basename =f'xz_{text})'),
                     colormap='twilight')
        imageXZ.scale = (self.dz.val/self.dxy.val, 1)
        
        self.viewer.add_image(im_xy,
                     name=self.gen.write_name(basename = f'xy_{text}'),
                     colormap='twilight') 
        