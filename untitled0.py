# -*- coding: utf-8 -*-
"""
Created on Tue May 22 10:10:24 2018

@author: pliniobas
"""

def teste(**k):
    if k.has_key('teste'):
        print('abacate')
        
teste(teste=None)
