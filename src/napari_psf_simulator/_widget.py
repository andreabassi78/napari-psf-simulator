# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 00:16:58 2022

@author: Andrea Bassi @Polimi
"""
from .psf_generator import PSF_simulator
from .gui_utils import Setting
from qtpy.QtCore import Qt
from qtpy.QtWidgets import QSplitter, QWidget, QPushButton
from qtpy.QtWidgets import QComboBox,QLabel, QVBoxLayout, QCheckBox
from napari.qt.threading import thread_worker
import numpy as np
 
class Psf_widget(QWidget):
    
    '''
    Napari widget to simulate the point spread function of an imaging system,
    particularly an optical microscope.
    Aberrations currently included are:
        -Zernike phase aberration
        -presence of a layer with different refractive index than the sample 
    '''
    
    ABERRATION_DICT = {0:'No aberrations', 1:'Slab aberrations', 2:'Zernike aberrations'}
    
    def __init__(self, napari_viewer,
                 ):
        self.viewer = napari_viewer
        super().__init__()
      
        self.settings_dict = {'NA': 0.65,
                             'n': 1.33,
                             'wavelength [\u03BCm]': 0.532,
                             'Nxy':127,
                             'Nz':63,
                             'dxy [\u03BCm]':0.05,
                             'dz [\u03BCm]':0.15
                             }
        self.slab_aberrations_dict = {'n1':1.51,
                                      'thickness [\u03BCm]':100.0,
                                      'alpha [deg]':0.0,
                                      }
        self.zernike_aberrations_dict = {'N':3,
                                         'M':1,
                                         'weight [\u03BB]':1.0
                                         }
                
        self.setup_ui() 
        self.initialize_simulator()
    
    def setup_ui(self): 
        """
        Creates the user interface. 
        Note that the magicgui is not used and the Setting class (gui_utils.Setting)  
        """
        
        def add_section(_layout,_title):
            splitter = QSplitter(Qt.Vertical)
            _layout.addWidget(splitter)
            _layout.addWidget(QLabel(_title))
            
        # initialize layout
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # add Settings
        add_section(layout,'Settings')
        settings_layout = QVBoxLayout()
        layout.addLayout(settings_layout)
        self.create_Settings(self.settings_dict, settings_layout)
        
        # add plot results checkbox
        self.plot_checkbox = QCheckBox("Plot PSF profile in Console")
        self.plot_checkbox.setChecked(False)
        layout.addWidget(self.plot_checkbox)
        # add show Airy disk checkbox
        self.airy_checkbox = QCheckBox("Show Airy disk")
        self.airy_checkbox.setChecked(True)
        layout.addWidget(self.airy_checkbox)
        # add calculate psf button
        calculate_btn = QPushButton('Calculate PSF')
        calculate_btn.clicked.connect(self.calculate_psf)
        layout.addWidget(calculate_btn)
        
        # Aberrations
        add_section(layout,'Aberrations')
        self.aberration_combo = QComboBox()
        self.aberration_combo.addItems(list(self.ABERRATION_DICT.values()))
        self.aberration_combo.currentIndexChanged.connect(self.initialize_simulator)
        layout.addWidget(self.aberration_combo)
        
        # Slab aberrations
        add_section(layout,'Slab aberrations')
        slab_layout = QVBoxLayout()
        layout.addLayout(slab_layout)
        self.create_Settings(self.slab_aberrations_dict, slab_layout)
        
        # Zernike aberrations
        add_section(layout,'Zernike aberrations')
        aberration_layout = QVBoxLayout()
        layout.addLayout(aberration_layout)
        self.create_Settings(self.zernike_aberrations_dict, aberration_layout) 
        
        # Other aberrations .... to be implemented
        # add_section(layout,'other aberrtion')
        
    def create_Settings(self,  s_dict,  layout):
        """
        Adds the settings whose name and value is specified in a dictionary to a layout        

        Parameters
        ----------
        s_dict : dict
            Dictonary containing the settings name and initial value that will be displayed
        layout : Qwidget
            Layout (Qwidget) where the settings will be created

        """
        def find_between(s, start, end):
            return (s.split(start))[1].split(end)[0]
            
        def find_decimals_num(s):
            return len(str(s).split(".")[1])
            
        for key, val in s_dict.items():
            name_unit = key.split(" [")
            name = name_unit[0]
            unit = '' if len(name_unit)==1 else find_between(key, '[', ']')
            decimals = find_decimals_num(val) if type(val)==float else 0
            
            new_setting = Setting(name=name, dtype=type(val), initial=val, 
                                  spinbox_decimals=decimals, spinbox_step= 0.1,
                                  unit = unit,
                                  layout=layout,
                                  write_function=self.initialize_simulator)
            setattr(self, name, new_setting)          
    
    def initialize_simulator(self):
        '''
        Starts the simulator defining the space in the spatial frequency domain.
        Creates the systems' pupil (self.gen.amplitude and self.gen.phase).
        Sets the current aberration type.
        '''
        self.gen = PSF_simulator(self.NA.val, self.n.val, self.wavelength.val,
                      self.Nxy.val , self.Nz.val, dr = self.dxy.val, dz = self.dz.val)
        self.gen.generate_kspace()
        
        active_aberration = self.aberration_combo.currentIndex()
        self.add_aberration(active_aberration)
    
    def add_aberration(self, value):
        '''
        Adds a phase aberration to the pupil (adding a phase term to self.gen.phase)
        '''
        
        if value == 1:
            self.gen.add_slab_scalar(self.n1.val, self.thickness.val, self.alpha.val)
        if value == 2:
            self.gen.add_Zernike_aberration(self.N.val, self.M.val, self.weight.val)
    
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
            
            self.viewer.add_image(psf,
                             name=self.gen.write_name(basename = 'stack'),
                             colormap='twilight')
            
            if self.airy_checkbox.checkState():
                self.show_airy_disk()
            
            if self.plot_checkbox.checkState():
                self.gen.plot_psf_profile()  
        
        worker = generator()  # create "worker" object
        worker.returned.connect(update_image)  # connect callback functions
        worker.start()  # start the thread
        
    def show_airy_disk(self):
       '''
       Shows the diffraction limited Airy disk as 3 perpendicular ellipses 
       '''
       deltaR = 1.22*self.wavelength.val/self.NA.val/2 # Rayleigh resolution
       deltaZ = self.wavelength.val/self.n.val/(1-np.sqrt(1-self.NA.val**2/self.n.val**2)) # Diffraction limited axial resolution
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
       self.viewer.dims.current_step = (posz,posxy,posxy) # shows the image center of the stack in 3D
       
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
        
if __name__ == '__main__':
    import napari
    viewer = napari.Viewer()
    widget = Psf_widget(viewer)
    viewer.window.add_dock_widget(widget,
                                  name = 'PSF Simulator @Polimi',
                                  add_vertical_stretch = True)
    napari.run()      