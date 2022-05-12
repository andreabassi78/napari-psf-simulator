# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:24:39 2022

@author: Andrea Bassi
"""
from qtpy.QtWidgets import QLabel, QFormLayout, QSpinBox, QDoubleSpinBox, QCheckBox
from qtpy.QtCore import Qt

class Setting():
    '''
    Auxiliary class to create an numerical or boolean attribute 
    with a corresponding Qwidget (QSpinBox, QDoubleSpinBox or QCheckBox),
    and update its value as a property (self.val). 
    '''
    
    def __init__(self, name ='setting_name',
                 dtype = int,
                 initial = 0,
                 vmin = 0,
                 vmax = 2**16-1,
                 spinbox_decimals=3,
                 spinbox_step=0.05,
                 width = 150,
                 unit = '',
                 layout = None,
                 write_function = None,
                 read_function = None):
        '''
        Parameters
        ----------
        name : str
            Name of the Setting and label shown on the corresponding QWidget.
        dtype : type
            Type of the Setting. Currently supported for int, float and bool.
        initial : int, float or bool
            Initial value of the Setting, stored in the @property self.val. The default is 0.
        vmin : int, float, optional
            Minimum value. The default is 0.
        vmax : int, float, optional
            Maximum value. The default is 2**16-1.
        spinbox_decimals : int
            For objects with dtype==float, it is the number of decimals to show in the DoubleSpinBox. The default is 3.
        spinbox_step : float
            For objects with dtype==float, it is the step of the DoubleSpinBox. The default is 0.05.
        width : int
            Width of the spinbox. The default is 150.
        unit : str
            Unit of measurement of the Setting. The default is ''.
        layout : QWidget
            Parent QWidget layout where the Setting will be shown. The default is None.
        write_function : function or method
            Function/method that is executed on value change of the QWidget
        read_function : function or method
            not implemented

        '''
        self.name= name
        self._val = initial
        self.dtype = dtype
        self.spinbox_decimals = spinbox_decimals
        self.spinbox_step = spinbox_step
        self.unit = unit
        self.write_function = write_function
        # self.read_function = read_function
        self.create_spin_box(layout, dtype, vmin, vmax, unit, width)
        
    def __repr__(self):
        return f'{self.name} : {self._val}'
        
    @property    
    def val(self):
        self._val = self.get_func()
        return self.dtype(self._val) 
    
    @val.setter 
    def val(self, new_val):
        new_val = self.dtype(new_val)
        self.set_func(new_val)
        self._val = new_val
        
    def create_spin_box(self, layout, dtype, vmin, vmax, unit, width):
        name = self.name
        val = self._val
        if dtype == int:
            sbox = QSpinBox()
            sbox.setMaximum(vmax)
            sbox.setMinimum(vmin)
            self.set_func = sbox.setValue
            self.get_func = sbox.value
            sbox.setSuffix(unit)
            change_func = sbox.valueChanged
            sbox.setFixedWidth(width)
        elif dtype == float:
            sbox = QDoubleSpinBox()
            sbox.setDecimals(self.spinbox_decimals)
            sbox.setSingleStep(self.spinbox_step)
            sbox.setMaximum(vmax)
            sbox.setMinimum(vmin)
            sbox.setSuffix(' '+unit)
            self.set_func = sbox.setValue
            self.get_func = sbox.value
            change_func = sbox.valueChanged
            sbox.setFixedWidth(width)
        elif dtype == bool:
            sbox = QCheckBox()
            self.set_func = sbox.setChecked
            self.get_func = sbox.checkState
            change_func = sbox.stateChanged
        
        else: raise(TypeError, 'Specified setting type not supported')
        
        self.set_func(val)
        if self.write_function is not None:
            change_func.connect(self.write_function)
        settingLayout = QFormLayout()
        settingLayout.setFormAlignment(Qt.AlignLeft)
        lab = QLabel(name)
        lab.setWordWrap(False)
        settingLayout.addRow(sbox,lab)
        layout.addLayout(settingLayout)
        self.sbox = sbox  

