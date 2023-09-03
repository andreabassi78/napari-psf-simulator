# -*- coding: utf-8 -*-
"""
Created on Fri May 27 23:31:57 2022

@author: andrea
"""
from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QLabel,
    QSplitter,
    QVBoxLayout,
    )
from .gui_utils import Combo_box

class SwitchableSection:
    '''
    SwitchableSection class
    '''
    def __init__(self,
                 name,
                 choices,
                 baselayout,
                 on_change_function,
                 **kwargs):
                 
        '''
        Parameters
        ----------
        ....

        '''
        
        self.name = name
        self.choices = choices
        self.baselayout = baselayout
        self.on_change_function = on_change_function
        self.combo = Combo_box(name = name, choices = choices,
                                layout = baselayout, write_function= on_change_function)
        self.add_sub_layout()      

    def add_splitter(self, baselayout, title):
        splitter = QSplitter(Qt.Vertical)
        baselayout.addWidget(splitter)
        baselayout.addWidget(QLabel(title))

    def add_sub_layout(self):
        _layout = QVBoxLayout()

        self.baselayout.addLayout(_layout)

        self.sub_layout = _layout

    def remove_sub_layout_content(self):
        if hasattr(self, 'sub_layout'):
            _layout = self.sub_layout
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

        