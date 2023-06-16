#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 21:13:03 2023
@author: Hamzeh Mohammadigheymasi, Peidong Shi, Xiao Zhuowei
"""


import subprocess
subprocess.run(['python', './1-Data_Preparation.py'],)
subprocess.run(['python', './2-Run_EqT.py'],)
subprocess.run(['python', './3-Run_S-EqT.py'],)
subprocess.run(['python', './4-Run_REAL.py'],)
subprocess.run(['python', './5-Run_MIL.py'])

