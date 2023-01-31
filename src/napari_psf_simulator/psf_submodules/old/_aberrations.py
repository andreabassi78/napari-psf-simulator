# -*- coding: utf-8 -*-
"""
Created on Fri May 27 23:31:57 2022

@author: andrea
"""
from enum  import Enum 

    
   
class BaseAberration:
    
    def __init__(self, name='MY_ABERRATION',
                       phase_aberration_function = None,
                       **kwargs):
        self.name = name
        self.phase_aberration_function = phase_aberration_function
        self.settings,self.units = self.define_settings(kwargs)
        self._index = 0
  
    def define_settings(self, settings_with_units):
        units = {key:val for (key,val) in settings_with_units.items() if '_units' in key}
        settings = {key:val for (key,val) in settings_with_units.items() if '_units' not in key}
        return settings, units            
        
 
   
class Aberrations:
    
    __num = 0
    
    def __init__(self):
        pass
    
    @classmethod        
    def add_aberration(cls, aberration:BaseAberration):
        setattr(cls, aberration.name, aberration) 
        aberration._index = cls.__num 
        cls.__num +=1
    
    #@property
    @classmethod
    def enum(cls):
        attrs = {}
        for name in vars(cls):
            attr = getattr(cls, name)
            if type(attr) is BaseAberration:
                attrs[getattr(attr,'name')] = getattr(attr,'_index')   
        return Enum('Abber', attrs)
     
    @classmethod 
    def get_by_name(cls,name):
        aberration = getattr(cls, name)
        return aberration 
    
    @classmethod 
    def get_by_idx(cls,idx):
        enum = cls.enum()
        name = enum(idx).name
        aberration = getattr(cls, name)
        return aberration 
        