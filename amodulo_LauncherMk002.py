# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 09:40:41 2017

@author: VPNUser
"""

import subprocess as sp

l = [3,3,3,3,3,3,3,3] #dir 
#l = [7,7,7,7,7,7,7,7]
l = [1,2,3,4,5]
l = [1,2]
len(l)

p = []

#1º argv = CtrlDirNum
#2º argv = CtrlDirQuad - pode ser uma lista
#3º argv = thread
#4º argv = cpuname 

#cpuname = 'Lenovo Corei7 3632QM'
cpuname = 'FX8350'
#cpuname = 'Corei7 3770'

for ix,aux in enumerate(l):
    p.append(sp.Popen(("python","C:/Google Drive/0Mestrado Dissertacao/0Trabalhando Geral/amodulo.py",' 8',' %d'%aux,' %d'%ix,' %s'%cpuname),
    stdout=sp.PIPE,
    stdin=sp.PIPE,
    stderr=sp.PIPE))
    print(ix)

#%%
#p[0].communicate()
#dir(p[0])

#a = p[0].poll()
p[0].communicate()
