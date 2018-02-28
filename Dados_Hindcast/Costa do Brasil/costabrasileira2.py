# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 08:04:12 2017

@author: VPNUser
"""

import numpy as np

f = open('costabrasileira2.txt','r')
d = f.readlines()

lat = d[0].strip('\n').strip(' ').split(',')
lon = d[1].strip('\n').strip(' ').split(',')

temp = []
temp2 = []
for aux in range(len(lat)):
    try: 
        temp.append(float(lat[aux]))
        temp2.append(float(lon[aux]))
        pass
    except:
        print aux
        pass
    pass

lat = temp
lon = temp2

