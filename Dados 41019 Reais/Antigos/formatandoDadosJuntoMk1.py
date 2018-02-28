# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 18:42:05 2017

@author: Plinio Bueno Andrade Silva
"""

#==============================================================================
# YY	Year
# MM	Month
# DD	Day
# hh	Hour
# mm	Minute
# WDIR	Wind direction (the direction the wind is coming from in degrees clockwise from true N)
# WSPD	Wind speed, averaged over an eight-minute period (m/s)
# GST	Peak 5 or 8 second gust speed measured during the eight-minute or two-minute period(m/s)
# WVHT	Significant wave height (meters)
# DPD	Dominant wave period (seconds)
# APD	Average wave period (seconds)
# MWD	The direction from which the waves at the dominant period (DPD) are coming. The units are degrees from true North, increasing clockwise, with North as 0 (zero) degrees and East as 90 degrees.
# PRES	Sea level pressure (hPa)
# ATMP	Air temperature (Celsius)
# WTMP	Sea surface temperature (Celsius)
# DEWP	Dewpoint temperature (Celsius)
# VIS	visibility (nautical miles)
# TIDE	Water level in feet above or below Mean Lower Low Water, MLLW (feet)
#==============================================================================

lat = 27.537
lon = -62.945
latint =
lonint =


from datetime import datetime as dt
import numpy as np
import time

listAno = []




for aux in range(11,17,1):
#    listAno.append(dt.strptime('201%d'%aux,'%Y'))
    listAno.append('20%d.txt'%aux)
    pass


#==============================================================================
#%% Le o arquivo inicial
#==============================================================================
d = []
for aux in listAno:
    f = open(aux,'r')
    h = f.readline() #descarta primeira linha heading
    h2 = f.readline() #descarta segunga linha heading
    f.readline() #descarta terceira linha repetida do ano anterior
    temp = f.readlines() 
    for aux2 in temp:
        d.append(aux2)
    f.close()
    pass

#==============================================================================
#%% Formata a tabela de dados com split nos dados. Cria objeto datetime
#==============================================================================

d1 = []
for aux in d:
    temp = aux.split()
    temp2 = temp[0] + temp[1]+ temp[2]+ temp[3]+ temp[4]
#    temp3 = time.mktime(dt.strptime(temp2,'%Y%m%d%H%M').timetuple())
    temp3 = dt.strptime(temp2,'%Y%m%d%H%M')
    temp4 = [float(temp[aux2]) for aux2 in range(5,len(temp),1)]
    temp5 = [temp3] + temp4
    d1.append(temp5)
    pass

temp = []
temp = h.split()
temp[0] = 'data'


#==============================================================================
#%% Separa os maiores periodos sem interrupção nos dados de onda.
#==============================================================================

p = []
aux2 = 0
aux3 = False

for aux in d1:
    
    if aux[4] < 90.0 and aux[5] < 90.0 and aux[7] < 900.0:
#        print aux[0],aux[4],aux[5],aux[7]
        p.append(aux)
        pass
    else:
        if len(temp) > 100:
#            p.append(np.asarray(temp,dtype='float32'))
#            p.append(temp)
            pass
#        temp = []
        pass    
    pass



#for aux in xrange(len(d1)):
#    if float(aux[4])
#p[1][0]
#p[6][:,7].max()
#%%
#==============================================================================
# Separar a lista p com os maiores periodos de dados sem erro.
#==============================================================================

d3 = ( dict (date = [(aux[0].strftime('%Y/%m/%d %H')) for aux in p],
           hs = [aux[4] for aux in p],
           tp = [aux[5] for aux in p],
           dp = [aux[7] for aux in p]
           )
        )

f = open('41049w_%d.Dict.txt'%len(aux2),'w')
f.write(str(d3))
f.close()
        
    
