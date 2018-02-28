# -*- coding: utf-8 -*-
"""
Spyder Editor

"""

#no windows, tem que executar o script em modo admin (abrir CMD como admin)
#f = file(r'\\.\PhysicalDrive0','rb')

import sys
import os
import threading
import time

t0 = time.clock()
print time.clock()

progName = 'DistriProbGeralMk64.py'
CrtlMaxThreadNumbers = 8
CrtlFolderInsideDadosHindcast = 'Dados_Hindcast/Long -25 graus'

if '__file__' in globals():
    FILE = __file__



ncDirName = os.path.join(os.path.dirname(os.path.normpath(FILE)).replace('\\','/'),CrtlFolderInsideDadosHindcast)
progDirName = os.path.dirname(os.path.normpath(FILE)).replace('\\','/')
progName = progDirName + '/' + progName

ncFileList = [str(fil) for fil in os.listdir(ncDirName) if os.path.isfile(os.path.join(ncDirName, fil))]


class multi(threading.Thread):
    def __init__(self,progName,ncFileName,fileDir):
        threading.Thread.__init__(self)
        self.progName = '\"' + progName + '\"'
        self.ncFileName = '\"' + fileDir + ncFileName + '\"' 
        self.resultado = ''
        pass
    def run(self):
        print ('python ' + self.progName +' '+ self.ncFileName)
        os.system('python ' + self.progName +' '+ self.ncFileName)
        self.resultado = 'ok'
        pass
    pass

ThreadList = [multi(progName,aux,ncDirName) for aux in ncFileList]

#objetoThread = multi(progName,ncFileList[10])
#objetoThread.start()

#ThreadActive = [[] for _ in range(CrtlMaxThreadNumbers)]
ThreadActive = []

for aux in xrange(len(ThreadList)):
    ThreadActive.append(ThreadList[aux])
    ThreadActive[-1].start()
    print 'Thread iniciada: ',aux
    while len(ThreadActive) > CrtlMaxThreadNumbers-1:
        for aux2 in xrange(len(ThreadActive)):
            if ThreadActive[aux2].isAlive():
                pass
            else:
                del(ThreadActive[aux2])
                break
            time.sleep(1)
            pass
        pass
    pass
                
print time.clock() - t0
print 'fim'


#print(objetoThread.entered)

#threading.enumerate()


