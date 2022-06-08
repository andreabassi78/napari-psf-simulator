# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 16:34:41 2022

@author: Andrea Bassi @ Polimi
"""
from qtpy.QtWidgets import QLabel, QFormLayout, QSpinBox, QDoubleSpinBox, QCheckBox, QComboBox
from qtpy.QtCore import Qt
from enum import Enum, EnumMeta

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
                 spinbox_decimals=2,
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
            Parent QWidget layout where the Setting will be shown. The default is None. Needs to be specified
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
        
        else: raise(TypeError('Specified setting type not supported'))
        
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
    
    def set_min_max(self, vmin = 0, vmax = 3):
        vmax = self.dtype(vmax)
        vmin = self.dtype(vmin)
        self.sbox.setMaximum(vmax)
        self.sbox.setMinimum(vmin)


class Combo_box():
    '''
    Auxiliary class to create an combobox. 
    '''
    def __init__(self, name ='combo name',
                 initial = '_',
                 choices = ['_','one','two'],
                 userdata = [],
                 layout = None,
                 width = 150,
                 write_function = None,
                 read_function = None):
        '''
        Parameters
        ----------
        name : str
            Name of the combobox and label shown on the corresponding widget.
        initial : str
            Initial value of the combox. The default is "_".
        width : int
            Width of the combobox. The default is 150.
        choices : list(str) or Enum
            Combobox choices.
        userdata : list
            Additional data stored in the combobox.
        layout : QWidget
            Parent QWidget layout where the combobox will be shown.
            The default is None. Needs to be specified
        write_function : function or method
            Function/method that is executed on value change of the combobox
        read_function : function or method
            not implemented
        '''
        self.name = name
        self.write_function = write_function
        self.create_combo_box(name, choices, userdata, width, layout)
        self.choices = choices
    
    @property    
    def val(self):
        return self.combo.currentIndex()
    
    @property    
    def text(self):
        _text = self.combo.currentText()
        return str(_text) 
    
    @property 
    def current_data(self):
            _data = self.combo.currentData()
            return _data    
       
    def create_combo_box(self, name, choices, userdata, width, layout):
        combo = QComboBox()
        if type(choices) is list:  
            choices_names = choices
        elif type(choices) is EnumMeta:
            choices_names = choices._member_names_
            userdata = list(choices._value2member_map_.keys())
        assert len(userdata) in (0,len(choices_names)), f'Uncorrect userdata in {self.name} Combobox'    
        for idx, choice in enumerate(choices_names):
            shown_text = choice.replace('_',' ')
            if len(userdata) == len(choices):
                data = userdata[idx]
                combo.addItem(shown_text, userData=data)
            else:
                combo.addItem(shown_text)

        #combo.setEditable(True)
        #combo.lineEdit().setAlignment(Qt.AlignCenter)
        comboLayout = QFormLayout()
        comboLayout.setFormAlignment(Qt.AlignLeft)
        lab = QLabel(name)
        lab.setWordWrap(False)
        comboLayout.addRow(combo,lab)
        layout.addLayout(comboLayout)
        if self.write_function is not None:
            combo.currentIndexChanged.connect(self.write_function)
        # combo.setFixedWidth(width)
        self.combo = combo


class BaseEnum(Enum):
    
    @classmethod
    def to_dict(cls):
        """Returns a dictionary representation of the enum."""
        return {e.name: e.value for e in cls}
    
    @classmethod
    def keys(cls):
        """Returns a list of all the enum keys."""
        return cls._member_names_
    
    @classmethod
    def values(cls):
        """Returns a list of all the enum values."""
        return list(cls._value2member_map_.keys())


def add_timer(function):
    import time
    """
    Function decorator to mesaure the execution time of a function or a method.
    """ 
    def inner(*args,**kwargs):
        
        print(f'\nStarting method "{function.__name__}" ...') 
        start_time = time.time() 
        result = function(*args,**kwargs) 
        end_time = time.time() 
        print(f'Execution time for method "{function.__name__}": {end_time-start_time:.6f} s') 
        
        
        return result
    inner.__name__ = function.__name__
    return inner 
