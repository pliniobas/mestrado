# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 15:10:55 2018

@author: Plinio Bueno
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

x = np.linspace(0,1,10,endpoint = False)
y = 10**x

plt.plot(x,y)
plt.plot(x,np.log10(y))
#y = np.log10(temp)


