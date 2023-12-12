# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

Use waveform extraction function from phy to get waveform template
phy extract-waveforms params.py

"""
import os
from phy.apps.template import template_gui

os.chdir('Z:/ibn-vision/DATA/SUBJECTS/M23028/ephys/20230704/kilosort_probe_2')
os.chdir('Z:/ibn-vision/DATA/SUBJECTS/M23017/ephys/20230630/kilosort')
os.chdir('D:/Neuropixel_recording/DATA/SUBJECTS/M23028/ephys/20230706/kilosort_probe_1')



import os
from phy.apps.template import template_gui
#os.chdir('Z:/ibn-vision/DATA/SUBJECTS/M23028/ephys/20230703/kilosort_probe_1')
os.chdir('Z:/ibn-vision/DATA/SUBJECTS/M23017/ephys/20230630/kilosort')
#os.chdir('Z:/ibn-vision/DATA/SUBJECTS/M22069/ephys/20221201/kilosort')
template_gui("params.py")
