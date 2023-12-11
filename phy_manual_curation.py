# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 02:09:23 2022

@author: masah
"""
import os
#path = 'D:/Neuropixel_recording/M22069/20221130/kilosort'
#path = 'X:\ibn-vision\DATA\SUBJECTS\M22069\ephys\20221201\kilosort'
path = '/research/DATA/SUBJECTS/M22069/ephys/20221201/kilosort'
os.chdir(path)


from phy.apps.template import template_gui
template_gui("params.py")