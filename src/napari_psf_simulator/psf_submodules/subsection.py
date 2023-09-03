# -*- coding: utf-8 -*-
"""
Created on Fri May 27 23:31:57 2022

@author: andrea
"""
from enum  import Enum 
    
class BaseSubSection:
    '''
    Base Subsection class
    '''
    def __init__(self, name='MY_SUBSECTION',
                       subsection_function = None,
                       **kwargs):
        '''
        Parameters
        ----------
        name : str
            name of the aberration, it will also appean in the Combo_box
            showing the available aberrations in napari_psf_simulator._widget
        subsection_function : function
            function/method to be used when section is active 
        **kwargs : int or float or bool or str
            these are the Settings that will be created in the subsection layout of napari_psf_simulator._widget
            and the parameters passed to the subsection fucntion

        '''
        def _pass(*args,**kwargs):
            pass
        
        self.name = name
        if subsection_function == None:
            self.subsection_function = _pass
        else: 
            self.subsection_function = subsection_function
        
        self.settings,self.units = self.define_settings(kwargs)
        
    @staticmethod     
    def define_settings(settings_with_units):
        '''
        Defines if the value passed are Settings or their units.
        Considers a unit each passed value whose name is a str type and ends with '_units' 
        Parameters
        ----------
        settings_with_units : dict
            If the key ends with '_units' (and the value is a str type) 
            considers this a unit, otherwise considers it a Setting  
        '''
        units = {key:val for (key,val) in settings_with_units.items() if '_units' in key}
        settings = {key:val for (key,val) in settings_with_units.items() if '_units' not in key}
        return settings, units
   
class SubSection():
    """
    Class that contains multiple Base_Subsection, as attributes
    __num is a class attribute indicating the total number of available option in the section conmbo box 
    and it is used to index the aberrations.
    """
    __num = 0
    
    def __init__(self):
        
        self.add('None')
    
    def add(self, name='MY_SECTION',
                       subsection_function = None,
                       **kwargs):
        '''
        adds a new Baseaberration to the object  with the specified 
        name, phase_aberration_function and settings/units in kwargs
        '''
        self._add(BaseSubSection(name, subsection_function, **kwargs))
        
    def _add(self, subsection:BaseSubSection):
        setattr(self, subsection.name, subsection) 
        subsection._index = self.__num 
        self.__num +=1
    
    def enum(self):
        ''' Returns an Enum with the aberrations names and index'''
        attrs = {}
        for name in vars(self):
            attr = getattr(self, name)
            if type(attr) is BaseSubSection:
                attrs[getattr(attr,'name')] = getattr(attr,'_index')   
        return Enum('Abber', attrs)
     
    def get_by_name(self,name):
        ''' Returns the BaseSubSection by name '''
        subsection = getattr(self, name)
        return subsection 
    
    def get_by_idx(self,idx):
        ''' Returns the BaseSubSection by index '''
        enum = self.enum()
        name = enum(idx).name
        subsection = getattr(self, name)
        return subsection 
        