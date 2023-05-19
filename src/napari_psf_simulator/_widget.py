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
from .psf_submodules.switchable_layout import SwitchableSection
from .psf_submodules.aberrations import Aberrations 
from .psf_submodules.gui_utils import Combo_box, Setting
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
    
    def __init__(self, napari_viewer,
                 ):
        self.viewer = napari_viewer
        super().__init__()
        
        self.setup_ui()
        self.start_base_simulator()
        # self.setup_aberration_ui()
        
        
    def add_splitter(self, layout, title):
        splitter = QSplitter(Qt.Vertical)
        layout.addWidget(splitter)
        layout.addWidget(QLabel(title))
    
    
    def setup_ui(self): 
        """
        Creates the user interface. 
        Note that the magicgui is not used and the Setting class (gui_utils.Setting)  
        """
            
        # initialize layout
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # add Settings
        self.add_splitter(layout,'Settings')
        settings_layout = QVBoxLayout()
        layout.addLayout(settings_layout)
        self.create_base_Settings(settings_layout)
        
        # add plot results checkbox
        self.plot_checkbox = QCheckBox("Plot PSF profile in Console")
        self.plot_checkbox.setChecked(False)
        layout.addWidget(self.plot_checkbox)
        # add show Airy disk checkbox
        self.airy_checkbox = QCheckBox("Show Airy disk")
        self.airy_checkbox.setChecked(False)
        layout.addWidget(self.airy_checkbox)

        # add use generator combobox with switchable section
        self.add_splitter(layout,'PSF generator')
        self.generator_section = SwitchableSection(name = 'Generator',
                                      baselayout = layout, choices = self.generators,
                                      on_change_function = self.change_generator)
        
        self.add_splitter(layout,'Calculate')
        calculate_layout = QVBoxLayout()
        layout.addLayout(calculate_layout)
        # add calculate psf button
        calculate_btn = QPushButton('Calculate PSF')
        calculate_btn.clicked.connect(self.calculate_psf)
        layout.addWidget(calculate_btn)
        self.layout = layout

    def change_generator(self):
        self.generator_section.remove_sub_layout_content()
        selection_name = self.generator_section.combo.text
        print('Selection name:', selection_name)
        if selection_name == "Vectorial":
            self.add_PyFocus_settings(self.generator_section.sub_layout)
        self.start_simulator()

    def add_PyFocus_settings(self, layout):
        # add show intensity of each component checkbox

        self.add_splitter(layout, 'Pyfocus settings')
        self.pyFocus_component_checkbox = QCheckBox("Show X,Y,Z component Intensity")
        self.pyFocus_component_checkbox.setChecked(False)
        layout.addWidget(self.pyFocus_component_checkbox)
        self.lens_aperture = Setting(name='lens radius', dtype=float, initial=3, unit='mm', 
                          layout = layout, write_function = self.reinitialize_simulator)
        
        self.amplitude_section = SwitchableSection(name = 'Amplitude',
                                      baselayout = layout, choices = self.pyfocus_amplitudes,
                                      on_change_function = self.change_amplitude)
        self.change_amplitude()
        self.phase_section = SwitchableSection(name = 'Phase',
                                      baselayout = layout, choices = self.pyfocus_phases,
                                      on_change_function = self.change_phase)
        self.change_phase()
        self.add_splitter(layout, 'Polarization')
        self.gamma = Setting(name='gamma', dtype=float, initial=45, unit='deg', 
                          layout = layout, write_function = self.reinitialize_simulator)
        self.beta = Setting(name='beta', dtype=float, initial=90, unit='deg',
                          layout = layout, write_function = self.reinitialize_simulator)
        
    def change_amplitude(self):
        self.amplitude_section.remove_sub_layout_content()
        if self.amplitude_section.combo.text == 'Uniform':
            self.custom_mask = "1"    
        elif self.amplitude_section.combo.text == 'Gaussian':
            self.custom_mask = 'np.exp((rho/waist)**2)' #TODO: insert the correct value
            self.waist = Setting(name='waist', dtype=float, initial=2, unit='mm', 
                        layout = self.amplitude_section.sub_layout,
                        write_function = self.reinitialize_simulator)

    def change_phase(self):
        self.phase_section.remove_sub_layout_content()
        if self.phase_section.combo.text == 'Uniform':
            pass
            # self.custom_mask = '1' #TODO: insert the correct value 
        elif self.phase_section.combo.text == 'Vortex':
            self.custom_mask = 'np.exp(1j*phi)' #TODO: insert the correct value
            self.order = Setting(name='order', dtype=int, initial=1, 
                        layout = self.phase_section.sub_layout,
                        write_function = self.reinitialize_simulator) 

    
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
        self.wavelength = Setting(name='wavelength', dtype=float, initial=0.532, unit = '\u03BCm', spinbox_decimals = 3, 
                          layout = settings_layout, write_function = self.reinitialize_simulator)
        self.Nxy = Setting(name='Nxy', dtype=int, initial=51, vmin=1, vmax = 4095,
                          layout = settings_layout, write_function = self.reinitialize_simulator)
        self.Nz = Setting(name='Nz', dtype=int, initial=3, vmin = 1, vmax = 4095,
                          layout = settings_layout, write_function = self.reinitialize_simulator)
        self.dxy = Setting(name='dxy', dtype=float, initial=0.03, unit = '\u03BCm', 
                          layout = settings_layout, write_function = self.reinitialize_simulator)
        self.dz = Setting(name='dz', dtype=float, initial=0.50, unit = '\u03BCm',
                          layout = settings_layout, write_function = self.reinitialize_simulator)
        
        
    def setup_aberration_ui(self):
        '''
        Creates the aberration combobox
        '''
        self.setup_aberrations()
        self.add_splitter(self.layout,'Aberration')
        self.aberration_combo = Combo_box(name = 'Aberration', choices = self.aberrations.enum(),
                                layout = self.layout, write_function= self._on_aberration_change)
        

    def setup_aberrations(self):
        '''
        Defines the available aberrations that will appear in the combobox.
        To define your own aberration use self.aberration.add()
        indicating the:
            -name of the aberration
            -phase_aberration_function: the function to be executed when the aberration is selected
            -the Settings visible when the aberration is selected.
             These will also be the parameters of the phase_aberration_function.
             See the Aberrations and BaseAberration classes in napari_psf_simulator.aberrations
             
        '''
        self.aberrations = Aberrations()

        self.aberrations.add(name = 'SLAB_ABERRATION',
                        phase_aberration_function = self.gen.add_slab_scalar,
                        n1 = 1.51,
                        thickness = 100.0, thickness_units = '\u03BCm', #um 
                        alpha = 0.0,
                        alpha_units = 'deg',
                        )

        self.aberrations.add(name = 'ZERNIKE_ABERRATION',
                        phase_aberration_function = self.gen.add_Zernike_aberration,
                        N=3, M=1, 
                        weight=0.6, weight_units = '\u03BB', #lambda
                        )

        # add additional aberrations
        # self.aberrations.add(name = 'MY aberrations',
        #                phase_aberration_function = my_function,
        #                setting0 = 0.0, setting0_units = 'mm',
        #                setting1 = True)
    
    def start_base_simulator(self):
        """Starts the base simulator, which uses the paraxial approximation
        """
        self.gen = PSF_simulator(self.NA.val, self.n.val, self.wavelength.val,
                        self.Nxy.val , self.Nz.val, dr = self.dxy.val, dz = self.dz.val)
    
    def start_simulator(self):
        """
        Starts the PSF generator selected by the combobox 
        """
        selected_generator = self.generator_section.combo.current_data
        self.gen = selected_generator(self.NA.val, self.n.val, self.wavelength.val,
                                self.Nxy.val , self.Nz.val, dr = self.dxy.val, dz = self.dz.val,
                                gamma = self.gamma.val, beta = self.beta.val,
                                custom_mask = self.custom_mask)
        self.setup_aberrations()
        self.reinitialize_simulator()
    
    def reinitialize_simulator(self):
        '''
        Starts the simulator defining the space in the spatial frequency domain.
        Creates the systems' pupil (self.gen.amplitude and self.gen.phase).
        Sets the current aberration type.
        '''
        
        self.gen.re_init(self.NA.val, self.n.val, self.wavelength.val,
                        self.Nxy.val , self.Nz.val, dr = self.dxy.val, dz = self.dz.val,
                        gamma=self.gamma.val, beta=self.beta.val, custom_mask=self.custom_mask)
        selection = self.aberration_combo.current_data
        self.selected_aberration  = self.aberrations.get_by_idx(selection)
        self.add_aberration(self.selected_aberration) 
    
    def add_aberration(self, selected_aberration):
        '''
        Adds a phase aberration to the pupil (adding a phase term to self.gen.phase)
        '''  
        params = selected_aberration.settings
        function = selected_aberration.phase_aberration_function
        function(**params)
    

    def _on_aberration_change(self, current_data):
        if hasattr(self, 'aberration_layout'):
            self.remove_layout(self.aberration_layout)     
        selection = self.aberration_combo.current_data
        self.selected_aberration  = self.aberrations.get_by_idx(selection)
        aberration_layout = QVBoxLayout()
        self.aberration_layout = aberration_layout
        self.layout.addLayout(aberration_layout)
        self.create_aberration_Settings(self.selected_aberration.settings,
                                        self.selected_aberration.units,                                                                        
                                        aberration_layout,
                                        self._on_aberration_setting_change) 
        self.reinitialize_simulator()
        
        
    def _on_aberration_setting_change(self, *args):
        for key, val in self.selected_aberration.settings.items():
            self.selected_aberration.settings[key] = getattr(self,key).val
        self.reinitialize_simulator()


    def create_aberration_Settings(self,  settings_dict, units_dict, layout, write_function):
        """
        Adds the settings whose name and value is specified in a dictionary,
        to the passed layout

        Parameters
        ----------
        settings_dict : dict
            Dictonary containing the settings name and initial value that will be displayed
        units_dict : dict
            Dictonary containing the units of measurement (if available) for the settings 
        layout : Qlayout
            Layout where the settings will be created
        write_function: function or method
            to be executed when the value of one the aberrations settings is changed

        """
        def find_decimals_num(s):
            return len(str(s).split(".")[1])
            
        for key, val in settings_dict.items():
            name = key 
            unit = ''
            for ukey, uval in units_dict.items():  
                if key in ukey:
                    unit = uval
            decimals = find_decimals_num(val) if type(val)==float else 0
            new_setting = Setting(name=name, dtype=type(val), initial=val, 
                                  spinbox_decimals=decimals, spinbox_step= 0.1,
                                  unit = unit,
                                  layout=layout,
                                  write_function=write_function)
            setattr(self, name, new_setting)    

    
    def remove_layout(self,_layout):
        def delete_items(layout):
            if layout is not None:
                while layout.count():
                    item = layout.takeAt(0)
                    widget = item.widget()
                    if widget is not None:
                        widget.setParent(None)
                    else:
                        delete_items(item.layout())
        delete_items(_layout)
        self.layout.removeItem(_layout)
        delattr(self, 'aberration_layout')

    
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
                                colormap='twilight')
                self.rescaleZ(x_layer)
                y_layer = self.viewer.add_image(np.abs(self.gen.field.Ey)**2,
                                name=self.gen.write_name(basename = 'Ey'),
                                colormap='twilight')
                self.rescaleZ(y_layer)
                z_layer = self.viewer.add_image(np.abs(self.gen.field.Ez)**2,
                                name=self.gen.write_name(basename = 'Ez'),
                                colormap='twilight')
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
        