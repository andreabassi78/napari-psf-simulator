from enum import Enum
from qtpy.QtWidgets import QCheckBox
from ..psf_submodules.gui_utils import Setting, SwitchableSection

class PyFocusSettingsHandler:
    '''Auxiliary class for generating the combo boxes for amplitude, polariation and aberrations'''
    pyfocus_amplitudes = Enum('py_amplitude', 
                        {"Uniform": 0,
                        "Gaussian": 1
                        })

    pyfocus_phases = Enum('py_phase', 
                        {"Uniform": 0,
                        "Vortex": 1
                        })

    pyfocus_polarizations = Enum('py_polarization', 
                        {"X linear": 0,
                        "Y linear": 1,
                        "Right handed circular": 2,
                        "Left handed circular": 3,
                        "Custom": 4
                        })
    
    def __init__(self, widget):
        self.widget=widget
    
    def setup_pyfocus_default_settings_values(self):
        """Sets the default pyfocus parameters in the widget"""
        self.widget.gamma = 0
        self.widget.beta = 0
        self.widget.custom_amplitude = 1
        self.widget.custom_phase = 0

    def add_PyFocus_settings(self, layout):
        """Adds the settings portion in the widget"""
        # add show intensity of each component checkbox
        self.widget.add_splitter(layout, 'Pyfocus settings')
        self.widget.pyFocus_component_checkbox = QCheckBox("show x,y,z intensities")
        self.widget.pyFocus_component_checkbox.setChecked(False)
        layout.addWidget(self.widget.pyFocus_component_checkbox)
        self.widget.lens_aperture = Setting(name='lens radius', dtype=float, initial=3, unit='mm', 
                            layout = layout, write_function = self.widget.reinitialize_simulator)
        
        self.widget.amplitude_section = SwitchableSection(name = 'amplitude',
                                        baselayout = layout, choices = self.pyfocus_amplitudes,
                                        on_change_function = self.change_amplitude)
        self.widget.phase_section = SwitchableSection(name = 'phase',
                                        baselayout = layout, choices = self.pyfocus_phases,
                                        on_change_function = self.change_phase)
        self.widget.add_splitter(layout, 'Polarization')
        self.widget.polarization_section = SwitchableSection(name = 'polarization',
                                        baselayout = layout, choices = self.pyfocus_polarizations,
                                        on_change_function = self.change_polarization)


    def set_polarization(self):
        self.widget.gamma = self.widget.gamma_setting.val
        self.widget.beta = self.widget.beta_setting.val
        self.widget.reinitialize_simulator()

    def change_polarization(self):
        self.widget.polarization_section.remove_sub_layout_content()
        if self.widget.polarization_section.combo.text == "X linear":
            self.widget.gamma = 0
            self.widget.beta = 0
            self.widget.reinitialize_simulator()
        elif self.widget.polarization_section.combo.text == "Y linear":
            self.widget.gamma = 90
            self.widget.beta = 0
            self.widget.reinitialize_simulator()
        elif self.widget.polarization_section.combo.text == "Right handed circular":
            self.widget.gamma = 45
            self.widget.beta = 90
            self.widget.reinitialize_simulator()
        elif self.widget.polarization_section.combo.text == "Left handed circular":
            self.widget.gamma = 45
            self.widget.beta = -90
            self.widget.reinitialize_simulator()
        elif self.widget.polarization_section.combo.text == "Custom":
            self.widget.gamma_setting = Setting(name='gamma', dtype=float, initial=45, unit='deg', 
                            layout = self.widget.polarization_section.sub_layout, write_function = self.set_polarization)
            self.widget.beta_setting = Setting(name='beta', dtype=float, initial=90, unit='deg',
                            layout = self.widget.polarization_section.sub_layout, write_function = self.set_polarization)
            self.set_polarization()
        else:
            raise Exception(f"Combo text for polarization has a wrong value: {self.widget.polarization_section.combo.text}")

    def set_gaussian_amplitude(self):
        self.widget.custom_amplitude = f'np.exp(-(rho/{self.widget.waist.val*1000000})**2/2)' # The *1000000 factor comes from a passage from mm to nm
        self.widget.reinitialize_simulator()

    def change_amplitude(self):
        self.widget.amplitude_section.remove_sub_layout_content()
        if self.widget.amplitude_section.combo.text == 'Uniform':
            self.widget.custom_amplitude = "1"
            self.widget.reinitialize_simulator()
        elif self.widget.amplitude_section.combo.text == 'Gaussian':
            self.widget.waist = Setting(name='waist', dtype=float, initial=2, unit='mm', 
                        layout = self.widget.amplitude_section.sub_layout,
                        write_function = self.set_gaussian_amplitude)
            self.set_gaussian_amplitude()
        else:
            raise Exception(f"Combo text for amplitude has a wrong value: {self.widget.amplitude_section.combo.text}")

    def set_vortex_phase_mask(self):
        self.widget.custom_phase = self.widget.custom_mask = f'phi*{int(self.widget.order.val)}'
        self.widget.reinitialize_simulator()
        
    def change_phase(self):
        self.widget.phase_section.remove_sub_layout_content()
        if self.widget.phase_section.combo.text == 'Uniform':
            self.widget.custom_phase = '0' # Any constant value would do
            self.widget.reinitialize_simulator()
        elif self.widget.phase_section.combo.text == 'Vortex':
            self.widget.order = Setting(name='order', dtype=int, initial=1, 
                        layout = self.widget.phase_section.sub_layout,
                        write_function = self.set_vortex_phase_mask)
            self.set_vortex_phase_mask()
        else:
            raise Exception(f"Combo text for phase has a wrong value: {self.widget.phase_section.combo.text}")
