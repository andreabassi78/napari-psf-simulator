# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 00:16:58 2022

@author: Andrea Bassi @Polimi
"""
from .psf_generator import PSF_simulator
from qtpy.QtCore import Qt
from qtpy.QtWidgets import QSplitter, QWidget, QPushButton
from qtpy.QtWidgets import QComboBox,QLabel, QFormLayout, QVBoxLayout, QSpinBox, QDoubleSpinBox, QCheckBox
from napari.qt.threading import thread_worker
import numpy as np


class Setting():
    ''' 
    Auxilliary class to create an numerical attribute with a corresponding Qwidget,
    and update its value as a property (self.val)-
    - name of the QWidget (it contain a label)
    - dtype: Currently supported for int and float 
    - initial_value: stored in the @property self.val
    - vmin, vmax: min and max values of the QWidget
    - layout: parent Qlayout    
    - read function: not implemented
    - write_function is executed on value change of the QWidget
    
    '''
    
    def __init__(self, name ='setting_name',
                 dtype = int,
                 initial_value = 0,
                 vmin = 0,
                 vmax = 2**16-1,
                 layout = None,
                 write_function = None,
                 read_function = None):
        
        self.name= name
        self._val = initial_value
        self.write_function = write_function
        self.read_function = read_function
        self.create_spin_box(layout, dtype, vmin, vmax)
        
    @property    
    def val(self):
        self._val = self.sbox.value()
        return self._val 
    
    @val.setter 
    def val(self, new_val):
        self.sbox.setValue(new_val)
        self._val = new_val
        
    def create_spin_box(self, layout, dtype, vmin, vmax):
        name = self.name
        val = self._val
        if dtype == int:
            sbox = QSpinBox()
            sbox.setMaximum(vmax)
            sbox.setMinimum(vmin)
        elif dtype == float:
            sbox = QDoubleSpinBox()
            sbox.setDecimals(3)
            sbox.setSingleStep(0.1)
            sbox.setMaximum(2**16-1)
        
        else: raise(TypeError, 'Specified setting type not supported')
        sbox.setValue(val)
        if self.write_function is not None:
            sbox.valueChanged.connect(self.write_function)
        settingLayout = QFormLayout()
        settingLayout.addRow(QLabel(name), sbox)
        layout.addLayout(settingLayout)
        self.sbox = sbox
        
class Psf_widget(QWidget):
    
    ABERRATION_DICT = {0:'No aberrations', 1:'Slab aberrations', 2:'Zernike aberrations'}
    
    def __init__(self, napari_viewer,
                 ):
        self.viewer = napari_viewer
        super().__init__()
        
        self.settings_dict = {'NA':float(0.5),
                         'n':float(1.00),
                         'wavelength':float(0.5320),
                         'Nxy':int(127),
                         'Nz':int(63),
                         'dxy':float(0.10),
                         'dz':float(0.2)
                         }
        self.slab_aberrations_dict = {'n1':float(1.51),
                                      'thickness':float(100.0),
                                      'alpha':float(0.0),
                                      }
                
        self.zernike_aberrations_dict = {'N':int(3),
                                         'M':int(1),
                                         'weight':float(1.0)
                                         }
                
        self.setup_ui() 
        self.initialize_simulator()
    
    def setup_ui(self):     
        
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
        self.create_Settings(settings_layout, self.settings_dict)
        
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
        self.create_Settings(slab_layout, self.slab_aberrations_dict)
        
        # Zernike aberrations
        add_section(layout,'Zernike aberrations')
        aberration_layout = QVBoxLayout()
        layout.addLayout(aberration_layout)
        self.create_Settings(aberration_layout, self.zernike_aberrations_dict) 
        
        # Other aberrations ....
        add_section(layout,'')
        
    def create_Settings(self, slayout, s_dict):
        for key, val in s_dict.items():
            new_setting = Setting(name=key, dtype=type(val), initial_value=val,
                                   layout=slayout,
                                   write_function=self.initialize_simulator)
            setattr(self, key, new_setting)          
    
    def initialize_simulator(self):
        self.gen = PSF_simulator(self.NA.val, self.n.val, self.wavelength.val,
                      self.Nxy.val , self.Nz.val, dr = self.dxy.val, dz = self.dz.val)
        self.gen.generate_kspace()
        
        active_aberration = self.aberration_combo.currentIndex()
        self.add_aberration(active_aberration)
    
    def add_aberration(self, value):
        if value == 1:
            self.gen.add_slab_scalar(self.n1.val, self.thickness.val, self.alpha.val)
        if value == 2:
            self.gen.add_Zernike_aberration(self.N.val, self.M.val, self.weight.val)
    
    def calculate_psf(self):
        
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
       
    def _show_PSF_projections(self): #not used 
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