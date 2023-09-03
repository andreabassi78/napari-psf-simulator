# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 00:16:58 2022

@author: Andrea Bassi @Polimi
"""

from enum import Enum
from traceback import format_exc
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
from .auxiliary_handlers.PyFocus_auxiliary_handler import PyFocusSettingsHandler
from PyFocus.napari_adapter.napari_adapter import PyFocusSimulator

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
    
    aberrations = { # Passing a class as the key to a dictionary might be werid but it works just fine
        PyFocusSimulator: Enum('Aberrations', 
                       {"None": 0,
                        "Interface":1,
                        "Zernike":2
                        }),
        PSF_simulator: Enum('Aberrations',
                       {"None": 0,
                        "Slab":1,
                        "Zernike":2
                        })
    }
    
    # Default internal values for psf_generator for the resolution of the calculation
    Nxy_scalar = 201
    Nz_scalar = 201
    
    def __init__(self, napari_viewer,
                 ):
        self.viewer = napari_viewer
        super().__init__()
        
        self.aberration_changer_mapper = {
            PyFocusSimulator: self.change_vectorial_aberration,
            PSF_simulator: self.change_scalar_aberration
        }
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
        self.add_splitter(self.generator_section.sub_layout,'Aberrations', color=(2,10,10))
        self.aberration_section = SwitchableSection(name = 'Aberrations',
                                      baselayout = self.generator_section.sub_layout, choices = self.aberrations[self.generator_section.combo.current_data],
                                      on_change_function = self.aberration_changer_mapper[self.generator_section.combo.current_data])

        # add PSF calculation section
        self.add_splitter(layout,'Calculate Point Spread Function',color=(2,10,10))
        calculate_layout = QVBoxLayout()
        layout.addLayout(calculate_layout)
        calculate_btn = QPushButton('Calculate PSF')
        calculate_btn.clicked.connect(self.calculate_psf)
        layout.addWidget(calculate_btn)

    def change_generator(self):
        self.generator_section.remove_sub_layout_content()
        selected_generator = self.generator_section.combo.current_data
        if selected_generator is PSF_simulator:
            pass
        elif selected_generator is PyFocusSimulator:
            self.settings_handler = PyFocusSettingsHandler(widget = self)
            self.settings_handler.setup_pyfocus_default_settings_values()
            self.settings_handler.add_PyFocus_settings(self.generator_section.sub_layout)
        self.add_splitter(self.generator_section.sub_layout,'Aberrations', color=(2,10,10))
        self.aberration_section = SwitchableSection(name = 'Aberrations',
                                    baselayout = self.generator_section.sub_layout, choices = self.aberrations[self.generator_section.combo.current_data],
                                    on_change_function = self.aberration_changer_mapper[self.generator_section.combo.current_data])
        self.reinitialize_simulator()

    def change_vectorial_aberration(self):
        self.aberration_section.remove_sub_layout_content()
        
        if self.aberration_section.combo.text == 'None':
            pass 
        
        elif self.aberration_section.combo.text == 'Interface':           
            self.n1 = Setting(name='n1', dtype=float, initial=1.51, 
                        layout = self.aberration_section.sub_layout,
                        write_function = self.reinitialize_simulator) 
            self.axial_position = Setting(name='axial position', dtype=float, initial=0.0, unit = 'nm', 
                        layout = self.aberration_section.sub_layout,
                        write_function = self.reinitialize_simulator)
        
        elif self.aberration_section.combo.text == 'Zernike':           
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
 

    def change_scalar_aberration(self):
        self.aberration_section.remove_sub_layout_content()
        
        if self.aberration_section.combo.text == 'None':
            pass 
        
        elif self.aberration_section.combo.text == 'Slab':           
            self.n1 = Setting(name='n1', dtype=float, initial=1.51, 
                        layout = self.aberration_section.sub_layout,
                        write_function = self.reinitialize_simulator) 
            self.thickness = Setting(name='thickness', dtype=float, initial=100.0, unit = '\u03BCm', 
                        layout = self.aberration_section.sub_layout,
                        write_function = self.reinitialize_simulator)
            self.alpha = Setting(name='alpha', dtype=float, initial=0.0, unit = 'deg', 
                        layout = self.aberration_section.sub_layout,
                        write_function = self.reinitialize_simulator)
        
        elif self.aberration_section.combo.text == 'Zernike':           
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
        self.fov_xy = Setting(name='FOV xy', dtype=float, initial=1.5, unit = '\u03BCm',
                          layout = settings_layout, write_function = self.reinitialize_simulator)
        self.fov_z = Setting(name='FOV z', dtype=float, initial=3.0, unit = '\u03BCm',
                          layout = settings_layout, write_function = self.reinitialize_simulator)
        self.dxy = Setting(name='dxy', dtype=float, initial=0.03, unit = '\u03BCm', 
                          layout = settings_layout, write_function = self.reinitialize_simulator)
        self.dz = Setting(name='dz', dtype=float, initial=0.5, unit = '\u03BCm',
                          layout = settings_layout, write_function = self.reinitialize_simulator)
        self.lens_aperture = Setting(name='lens radius', dtype=float, initial=3, unit='mm', 
                            layout = settings_layout, write_function = self.reinitialize_simulator)

        
    
    def start_base_simulator(self):
        """Starts the base simulator, using scalar propagation
        """
        self._calculate_matrix_dimensions()
        self.gen = PSF_simulator(self.NA.val, self.n.val, self.wavelength.val,
                            self.Nxy_scalar , self.Nz_scalar, dr = self.dxy.val, dz = self.dz.val, crop_Ns = (self.Nxy_show, self.Nz_show))
    

    def reinitialize_simulator(self):
        '''
        Starts the simulator, accorgin to choice made on the combo box,
        Sets the current aberration type.
        '''
        selected_generator = self.generator_section.combo.current_data
        self._calculate_matrix_dimensions()
        
        if selected_generator is PyFocusSimulator:
            self.gen = selected_generator(NA=self.NA.val, n=self.n.val, wavelength=self.wavelength.val, lens_aperture=self.lens_aperture.val,
                                Nxy=self.Nxy_show , Nz=self.Nz_show, dr=self.dxy.val, dz=self.dz.val,
                                gamma = self.gamma, beta = self.beta, incident_amplitude = self.custom_amplitude, 
                                incident_phase = self.custom_phase)
            self.settings_handler.set_incident_energy_ratio()
            self.add_vectorial_aberration()
        elif selected_generator is PSF_simulator:
            self.gen = selected_generator(self.NA.val, self.n.val, self.wavelength.val,
                                self.Nxy_scalar , self.Nz_scalar, dr = self.dxy.val, dz = self.dz.val, crop_Ns = (self.Nxy_show, self.Nz_show))
            self.add_scalar_aberration() 
    
    def _calculate_matrix_dimensions(self):
        """Calculates the dimmensions of the field to show. For scalar approximation, this is usefull since the calculation is done internally with a bigger number of divisions
        Note: This pass should be performed by the psf_generator, but I was unable to perform it there since the field returned was null for certain values of Nxy and Nz and the cause was not found
        """
        self.Nxy_show = int(self.fov_xy.val//self.dxy.val)
        self.Nz_show = int(self.fov_z.val//self.dz.val)
        # The number of points must be odd so that the center pixel corresponds the the (0,0,0) coordinates
        if self.Nxy_show % 2 == 0: self.Nxy_show += 1
        if self.Nz_show % 2 == 0: self.Nz_show += 1
        # We make the fov to show not to be greater than the one we calculate
        if self.Nxy_show >= self.Nxy_scalar: self.Nxy_scalar = self.Nxy_show
        if self.Nz_show >= self.Nz_scalar: self.Nz_scalar = self.Nz_show

    def add_vectorial_aberration(self):
        '''
        Adds the effect of a phase aberration
        '''  
        if self.aberration_section.combo.text == 'None':
            pass
        elif self.aberration_section.combo.text == 'Interface':
            self.gen.add_interface(self.n1.val, self.axial_position.val)
        elif self.aberration_section.combo.text == 'Zernike':
            self.gen.add_Zernike_aberration(self.N.val, self.M.val, self.weight.val)  

    def add_scalar_aberration(self):
        '''
        Adds the effect of a phase aberration
        '''  
        if self.aberration_section.combo.text == 'None':
            pass
        elif self.aberration_section.combo.text == 'Slab':
            self.gen.add_slab_scalar(self.n1.val, self.thickness.val, self.alpha.val)
        elif self.aberration_section.combo.text == 'Zernike':
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
            
    def update_incident_energy_ratio(self):
        if self.generator_section.combo.current_data is PyFocusSimulator:
            try:
                if self.set_incident_energy_ratio_to_one.checkState():
                    self.gen.use_energy_ratio = False # Variable to allow the user to compensate for the inciding energy factor
                else:
                    self.gen.use_energy_ratio = True
            except Exception as e:
                print(format_exc())
                print("No se pudo")
                pass
        else:
            pass
    
    def calculate_psf(self):
        '''
        Calculates the 3D PSF and shows it as a stack
        '''
        @thread_worker
        def generator():
            import warnings
            warnings.filterwarnings('ignore')
            self.gen.generate_pupil()
            self.update_incident_energy_ratio()
            self.gen.generate_3D_PSF()
            return self.gen.PSF3D
        
        def update_image(psf):
            #self.viewer.dims.current_step = (0,0,0)
            
            psf_layer = self.viewer.add_image(psf,
                             name=self.gen.write_name(basename = 'I'),
                             colormap='twilight')
            self.rescaleZ(psf_layer)
            
            if self.generator_section.combo.current_data == PyFocusSimulator and self.pyFocus_component_checkbox.checkState():
                x_layer = self.viewer.add_image(np.abs(self.gen.field.Ex)**2,
                                name=self.gen.write_name(basename = 'Ix'),
                                colormap='gray')
                self.rescaleZ(x_layer)
                y_layer = self.viewer.add_image(np.abs(self.gen.field.Ey)**2,
                                name=self.gen.write_name(basename = 'Iy'),
                                colormap='gray')
                self.rescaleZ(y_layer)
                z_layer = self.viewer.add_image(np.abs(self.gen.field.Ez)**2,
                                name=self.gen.write_name(basename = 'Iz'),
                                colormap='gray')
                self.rescaleZ(z_layer)
            
            if self.airy_checkbox.checkState():
                self.show_airy_disk()
            
            if self.plot_checkbox.checkState():
                self.gen.plot_psf_profile()
            
            posxy = self.Nxy_show // 2
            posz = self.Nz_show // 2
            self.viewer.dims.set_point(axis=[0,1,2], value=(0,0,0)) #raises ValueError in napari versions <0.4.13
            self.viewer.dims.current_step = (posz,posxy,posxy) # shows the image center of the stack in 3D
            self.viewer.dims.axis_labels = ["z", "y","x"]
            
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
       
       posxy = self.Nxy_show.val//2
       posz = self.Nz_show.val//2
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
       
       if (self.dz.val*self.Nz_show.val)/2 > deltaZ:
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
        