# -*- coding: latin-1 -*-
"""
Created on Sun Nov 19 01:21:45 2017

@author: Plinio Bueno Andrade Silva
"""

from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
import numpy as np 
import sys
#from scipy.stats import (kurtosis,skew,norm,lognorm,
#gumbel_r,gumbel_l,dweibull,weibull_max,weibull_min,frechet_r,spearmanr,rayleigh,gamma)

from scipy.stats import norm,lognorm,weibull_min
from scipy.special import gamma
from scipy.optimize import root,fsolve
from scipy.optimize import curve_fit
#from scipy.interpolate import interp1d,interp2d
from scipy.integrate import quad,nquad
import scipy
#import pandas as pd
#from sympy import *
#from math import sqrt,log,exp
import time
#from numba import jit
import os
#from fractions import gcd
import gc
import pandas as pd
#import class_ln_wei
#from class_ln_wei import wei

# =============================================================================
# Classe mdc - metodo da distribuicao condicional
# =============================================================================


#%%
#import multiprocessing
#from multiprocessing import Process,Pool

class mdc():
    
    def __init__(s):
        s.name = 'Modelo Distribuicao Conjunta\n'
        pass
    
    
    
    def fit(s,H,T,**kwargs):
        #Criar distribuicao de H Lognormal
        Hscale,Hloc,Hshape = lognorm.fit(H,floc = 0) # fit wei
        s.fln_H = lognorm(Hscale,Hloc,Hshape) #cria a distribuicao
        s.fln_H_fit = 'xi = %.4f, lambda = %.4f'%(Hscale,np.log(Hshape))
        
        
        #Criar distribuicao de H Weibull
        Hscale,Hloc,Hshape = weibull_min.fit(H) #fit wei
        s.fwei_H = weibull_min(Hscale,Hloc,Hshape) #cria a distribuicao
        s.fwei_H_fit = 'shape = %.4f,loc = %.4f, scale = %.4f'%(Hscale,Hloc,Hshape)
        
        
###############################################################################  
#        #Trocando o fit da biblio weibull_min por wei
#        s.fwei_H = wei()
#        s.fwei_H.fit(H)
#        s.fwei_H_fit = 'lambda = %.4f, alpha = %.4f'%(s.fwei_H.lambW,s.fwei_H.alphaW)
###############################################################################

        s.H = H
        s.T = T
        
        #Determinar o hitograma
        s.bins = 20
        s.rangeH = [0,10]
        s.rangeT = [0,20]
        s.polydegree = 3
        
        if kwargs.has_key('bins'):
            s.bins = kwargs['bins']
        if kwargs.has_key('rangeH'):
            s.rangeH = kwargs['rangeH']
            pass
        if kwargs.has_key('rangeT'):
            s.rangeT = kwargs['rangeT']
            pass
        if kwargs.has_key('polydegree'):
            s.polydegree = kwargs['polydegree']
            pass

        #determinando o histograma para fazer o scatter.
        s.H_hist_y,s.H_hist_x = np.histogram(H,s.bins,density = True,range = s.rangeH)
        s.H_hist_xM = s.H_hist_x[:-1] + s.H_hist_x[0:2].mean()
        
        s.T_hist_y,s.T_hist_x = np.histogram(T,s.bins,density = True,range = s.rangeT)
        s.T_hist_xM = s.T_hist_x[:-1] + s.T_hist_x[0:2].mean()
        
                #Separando T condicional a H e calculando os parametros da distribuicao
        dft = pd.DataFrame(dict(H=H,T=T))
        
        ln_param = []
        wei_param = []
        
        # se for imprimir o grafico, mostra o ajuste de Hs em lognormal e weibull
        if kwargs.has_key('print_cdf') and kwargs['print_cdf'] == True:
            temp = np.linspace(0,max(s.H+1),20)
            figHs = plt.figure()
            ax = figHs.add_subplot(111)
            ax.plot(temp,s.fln_H.cdf(temp),'--',label='ln')
            ax.plot(temp,s.fwei_H.cdf(temp),':',label='wei')
            ax.plot(s.H_hist_xM,np.cumsum(s.H_hist_y/np.sum(s.H_hist_y)),'-',label = 'raw')
            ax.legend()
            ax.set_title('Fits de Hs com weibull ou lognormal')
        

       #se for imprimir o grafico (kwargs.has_key('print_cdf')), d e dd sao as colunas e linhas do grafico
       # determina d e dd que Ã© o tamanho do grafico a ser criado ax.subplot(d,dd,ddd)
        if kwargs.has_key('print_cdf') and kwargs['print_cdf'] == True:
            dd = 0
            d = 1
            while 24 / d > dd:
                dd = d
                d += 1
                pass
            
            figlog = plt.figure(figsize=(d*2,dd*2))
            x = np.linspace(s.rangeT[0],s.rangeT[1],100,endpoint = False)              
        ddd = 0 ### ddd Ã© a posicao que o grafico vai assumir na matrix definida por d e dd           
        
        #circula entre as classes de bins de hs
        for ix,aux in enumerate(s.H_hist_x): 
            ddd += 1
            
            if ix == len(s.H_hist_x) - 1:
                break
            temp = dft[np.logical_and(dft.H>s.H_hist_x[ix],dft.H<s.H_hist_x[ix+1])]['T'].values #dataframe com valores de tp condicionados a hs
            if len(temp) > 50:
                #lista contem [xi,loc,lamb*e,posicaoX,hs_condicionador]
                #xi = [0], loc=[1], lamb =[2], xpos = [3]
                ln_param.append(lognorm.fit(temp,floc=0) + tuple([s.H_hist_xM[ix]])) #acha os parametros da distribuiÃ§Ã£o lognormal
                #a lista contem [lamb,loc,alpha,hs_condicionador]
                #lambw=[0],loc=[1], alpha=[2], xpos=[3]
#                wei_param.append(weibull_min.fit(temp,floc=0) + tuple([s.H_hist_xM[ix]])) #acha os parametros da distribuiÃ§Ã£o weibull
                wei_param.append(weibull_min.fit(temp) + tuple([s.H_hist_xM[ix]])) #acha os parametros da distribuiÃ§Ã£o weibull

                #fazendo o print do cdf de cada tp condicionado a hs
                if kwargs.has_key('print_cdf') and kwargs['print_cdf'] == True:
                    ax = figlog.add_subplot(d,dd,ddd)

                    logtemp = lognorm(ln_param[-1][0],ln_param[-1][1],ln_param[-1][2])
                    weitemp = weibull_min(wei_param[-1][0],wei_param[-1][1],wei_param[-1][2])
                    
                    #print das cumulativas dos histogramas
                    temp_T_hist_y,temp_T_hist_x = np.histogram(temp,s.bins,density = True,range = s.rangeT)
                    temp_T_hist_xM = temp_T_hist_x[:-1] + temp_T_hist_x[0:2].mean()
                    
                    ax.plot(temp_T_hist_xM,np.cumsum(temp_T_hist_y/np.sum(temp_T_hist_y)),'-',label='hist')
                                       
                    #print das cumulativas das funcoes ajustadas
                    ax.plot(x,logtemp.cdf(x),'--',label = 'ln')
                    ax.plot(x,weitemp.cdf(x),':',label = 'wei')
                    ax.legend()
              
                    ax.set_title("Classe de Hs %.2f"%s.H_hist_x[ix])
                    figlog.tight_layout()                    
                    pass
                pass
            pass


        
###############################################################################
                #Trocando o fit da biblio weibull_min por wei
#                w = wei()
#                w.fit(temp)
#                wei_param.append((w.lambW,0,w.alphaW)+ tuple([s.H_hist_xM[ix]])) 
###############################################################################

        
        #Criando as funcoes dos parametros de distribuicao
#        s.Tfxi = np.poly1d(np.polyfit([aux[3] for aux in ln_param],[aux[0] for aux in ln_param],s.polydegree))
#        s.Tflamb = np.poly1d(np.polyfit([aux[3] for aux in ln_param],[aux[2] for aux in ln_param],s.polydegree))
#
        s.Tflambw = np.poly1d(np.polyfit([aux[3] for aux in wei_param],[aux[0] for aux in wei_param],s.polydegree))
        s.Tfalphaw = np.poly1d(np.polyfit([aux[3] for aux in wei_param],[aux[2] for aux in wei_param],s.polydegree))

  
        
        #Criando as funcoes dos parametros de distribuicao
        s.dl = pd.DataFrame(ln_param)
        s.dw = pd.DataFrame(wei_param)
        
        s.dl.columns = "shape,loc,scale,HsBinClass".split(',')
        s.dw.columns = "shape,loc,scale,HsBinClass".split(',')
        
        xdata = s.dl.loc[:,'HsBinClass'].values
        ydata = s.dl.loc[:,'shape'].values
#        
#        print(s.dl)
        #Função que define o xi = sigma = shape
        s.Tfxi = lambda h,a1,a2,b1,b2: a1 * np.exp(-b1*h) + a2 * np.exp(-b2*h)
        s.ln_func_param_xi,temp = curve_fit(s.Tfxi,xdata,ydata)
#        
        xdata = s.dl.loc[:,'HsBinClass'].values
        ydata = s.dl.loc[:,'scale'].values
        #Função que defune o lambda = sigma^2 = scale
#        s.Tflamb = lambda h,a1,a2,b1,b2: a1 * np.exp(-b1*h) + a2 * np.exp(-b2*h)
        s.Tflamb = lambda h,a1,a2,b: a1 + (a2 * (h ** b))
        s.ln_func_param_lamb,temp = curve_fit(s.Tflamb,xdata,ydata)


        
        
        #definindo a Função de distribuiÃ§Ã£o de Hs
        if kwargs.has_key('tipofH'):
            if kwargs['tipofH'] == 'weibull':
               s.fH = s.fwei_H
               s.fHtype = 'weibull'
            if kwargs['tipofH'] == 'lognormal':
                s.fH = s.fln_H
                s.fHtype = 'lognormal'
                pass
        else: 
            s.fH = s.fln_H
            s.fHtype = 'lognormal'
            pass
        
        #definindo a Função de distribuiÃ§Ã£o de Tp condicional ao Hs
        if kwargs.has_key('tipofT'):
            if kwargs['tipofT'] == 'weibull':
                s.fT = lambda h: weibull_min(s.Tflambw(h),0,s.Tfalphaw(h))
                s.fTtype = 'weibull'
            if kwargs['tipofT'] == 'lognormal':
                s.fT = lambda h: lognorm(s.Tfxi(h,*s.ln_func_param_xi),0,s.Tflamb(h,*s.ln_func_param_lamb))
                s.fTtype = 'lognormal'
        else:
            s.fT = lambda h: lognorm(s.Tfxi(h,*s.ln_func_param_xi),0,s.Tflamb(h,*s.ln_func_param_lamb))
            s.fTtype = 'lognormal'
            

    def pdf(s,h,t):
#        print(s.fH.pdf(h))
#        print(s.fT(h).pdf(t))
        temp = s.fH.pdf(h) * s.fT(h).pdf(t)
        if np.isnan(temp):
            return 0
        else:
            return temp
        
    def pdf_series(s,H,T):
        
        #parametros do histograma de H
        hs_hist_y,hs_hist_x = np.histogram(H,s.bins,density = True,range = s.rangeH)
        hs_hist_xM = hs_hist_x[:-1] + hs_hist_x[0:2].mean()
        
        #parametro do histograma de T
        tp_hist_y,tp_hist_x = np.histogram(T,s.bins,density = True,range = s.rangeT)
        tp_hist_xM = tp_hist_x[:-1] + tp_hist_x[0:2].mean()

        #parametros para plotagem dos dados da tabela df                            
        xpos,ypos = np.meshgrid(hs_hist_xM,tp_hist_xM)
        dfR = pd.DataFrame(np.ndarray((len(xpos),len(ypos))),index = tp_hist_xM, columns = hs_hist_xM)
        dfR[:] = np.nan
        
        for h in hs_hist_xM:
            for t in tp_hist_xM:
                dfR.loc[t,h] = s.pdf(h,t)
                #observar que as posicoes do indice realmente sao [t,h] porque h esta no eixo x
                pass
            pass
        #dfR = dfR.T
        return dict(xpos=xpos,ypos=ypos, zvalues = dfR)
        pass
        
      
    def integral(s):
        ran = [[s.H.min(),s.H.max()],[s.T.min(),s.T.max()]]
        s.i = nquad(lambda h,t: s.pdf(h,t),ran)
        print('calculo de integral terminado = ',s.i)
        return s.i


    # =============================================================================
    # Calculo de rho
    # =============================================================================

    def rho(s):
        ran = [[s.H.min(),s.H.max()],[s.T.min(),s.T.max()]]
        print('------------------------------------------')
        print('calculando media H')
        s.mh = nquad(lambda h,t: h* s.pdf(h,t),ran)[0]
        if s.H.mean() - s.mh > 1 or s.H.mean() - s.mh < -1:
            return 'erro no calculo de mh = %.6f'%s.mh
        print('calculando media H = %.6f'%s.mh)
        
        print('calculando media**2 H')
        s.m2h = nquad(lambda h,t: h** 2 * s.pdf(h,t),ran)[0]
        print('calculando media**2 H = %.6f'%s.m2h)        
    
        print('calculando media T')
        s.mt = nquad(lambda h,t: t * s.pdf(h,t),ran)[0]
        print('calculando media T = %.6f'%s.mt)

        print('calculando media **T')
        s.m2t = nquad(lambda h,t: t ** 2 * s.pdf(h,t),ran)[0]
        print('calculando media **T = %.6f'%s.m2t)
                
        print('calculando E')
        s.mht = nquad(lambda h,t: h * t * s.pdf(h,t),ran)[0]    
        print('calculando mht = %.6f'%s.mht)
        
        print('calculando desvio padrao H')
        s.dph = np.sqrt(s.m2h - s.mh**2) #desvio padrao
        print('calculando desvio padrao H = %.6f'%s.dph)

        print('calculando desvio padrao de T')
        s.dpt = np.sqrt(s.m2t - s.mt**2)    
        print('calculando desvio padrao de T = %.6f'%s.dpt)
        
        s.rhomdc = (s.mht - s.mh*s.mt) / (s.dph * s.dpt)
        print('calculo rho terminado = ',s.rhomdc )
        print('------------------------------------------')
        return s.rhomdc   
                 
    def printContour(s,**kwargs):
        
        #determinando serie distribuicao conjunta / nataf
        obj_hist = s.pdf_series(s.H,s.T)
        s.matriz = obj_hist

        #determinando serie bruta
        seriebruta, xedges, yedges = np.histogram2d(s.H, s.T, bins = s.bins, range=[s.rangeH,s.rangeT], normed=True)
        seriebruta = seriebruta.T
        xedges = xedges[:-1] + xedges[0:2].mean()
        yedges = yedges[:-1] + yedges[0:2].mean()
        xgrid,ygrid = np.meshgrid(xedges,yedges)
        
        #passango geral para escala log
#        xgrid = np.log10(xgrid)
#        ygrid = np.log10(ygrid)
#        seriebruta = np.log10(seriebruta)
#        obj_hist['zvalues'] = np.log10(obj_hist['zvalues'])
        
        
            
        #parametros de cores do contour
        cmap2 = cm.get_cmap('jet')
        cmap2.set_under('w')
        cmapbias = cm.get_cmap('seismic')
        levels = np.arange(0.001,0.3,0.01)
        levelsBias = np.arange(-0.11,0.11,0.001)
#        levels = np.arange(0.01,0.3,0.03)
#        levelsBias = np.arange(-0.11,0.11,0.02)
            
        
        #iniciando a figura
        fig = plt.figure(figsize=(12,6))
        
        ##############################################
        ax = fig.add_subplot(231)
        
        colbar = ax.contour(xgrid,ygrid, (seriebruta),cmap=cmap2,vmin = 0.001, vmax = 0.3,levels=levels)
        
        ax.set_xlabel('Hs')
        ax.set_ylabel('Tp')
        ax.set_title('Histograma Densidade Probabilidade\n Serie Bruta')
        
        yticks = obj_hist['zvalues'].index.values.astype(float)
        yticklabels = ['%.1f'%aux for aux in yticks]
        ax.set_yticks(yticks[::2])
        ax.set_yticklabels(yticklabels[::2])
                
        xticks = obj_hist['zvalues'].columns.values.astype(float)
        xticklabels = ['%.1f'%aux for aux in xticks]
        ax.set_xticks(xticks[::2])
        ax.set_xticklabels(xticklabels[::2])
#        #ax.set_xscale('log')   
#        #ax.set_yscale('log')
        ax.set_ylim(4,20)
             
        ax.grid()
        fig.colorbar(colbar, shrink=0.5, aspect=5,)
        
        ###############################################
        ax = fig.add_subplot(232)
        colbar = ax.contour(obj_hist['xpos'],obj_hist['ypos'],(obj_hist['zvalues']),cmap=cmap2, vmin = 0.001, vmax = 0.3,levels=levels)
        ax.set_xlabel('Hs')
        ax.set_ylabel('Tp')
        ax.set_title('Modelo Distribuicao Conjunta\nHs ajustado em %s\nTp ajustado em %s'%(s.fHtype,s.fTtype))
        
        ax.set_xticks(xticks[::2])
        ax.set_xticklabels(xticklabels[::2])
        ax.set_yticks(yticks[::2])
        ax.set_yticklabels(yticklabels[::2])
#        #ax.set_xscale('log')   
#        #ax.set_yscale('log')
        ax.set_ylim(4,20)
        
        ax.grid()
        fig.colorbar(colbar, shrink=0.5, aspect=5)

        
        ######################################################
        ax = fig.add_subplot(233)
        bias = (seriebruta) - (obj_hist['zvalues'])
        colbar = ax.contour(obj_hist['xpos'],obj_hist['ypos'], bias,cmap=cmapbias, vmin = -0.15, vmax = 0.15,levels=levelsBias)
        ax.set_xlabel('Hs')
        ax.set_ylabel('Tp')
        ax.set_title('Bias Brutos - MDC')
        ax.set_xticks(xticks[::2])
        ax.set_xticklabels(xticklabels[::2])
        ax.set_yticks(yticks[::2])
        ax.set_yticklabels(yticklabels[::2])
        
#        #ax.set_xscale('log')   
#        #ax.set_yscale('log')
        ax.set_ylim(4,20)
        
        ax.grid()
        fig.colorbar(colbar, shrink=0.5, aspect=5)   
        
       
        #plotando series unidas
        
        ax = fig.add_subplot(234)
        ax.plot(seriebruta.flatten(),label = 'Serie Bruta',ls='--')
        ax.plot(obj_hist['zvalues'].values.flatten(),label = 'Serie Ajustada',ls=':')
        
        xticks = np.linspace(0,s.bins**2,s.bins,endpoint=False)
        temp = obj_hist['zvalues'].index.values.astype(float)
        xticklabels = ['%.1f'%aux for aux in temp]
        
        ax.set_xticks(xticks[::2])
        ax.set_xticklabels(xticklabels[::2])
        ax.grid()
        ax.set_title(u'Secções de Hs em Função de Tp')
        ax.set_xlabel(u'Secções de Tp')
        ax.set_ylabel('Densidade Probabilidade')
        ax.legend()

        # =============================================================================
        #  Escrevendo parametros da funcao na figura.       
        # =============================================================================
        
        ax  = fig.add_subplot(235)
        temp = 0.94
        
        if kwargs.has_key('filename'):
            filename = kwargs['filename']
            ax.text(0.05,temp,u'Arquivo:\n%s'%filename)
            temp -= 0.06
        
        if kwargs.has_key('dirI') and kwargs.has_key('dirM') and kwargs.has_key('dirF'):
            dirI = str(kwargs['dirI'])
            dirM = str(kwargs['dirM'])
            dirF = str(kwargs['dirF'])
            ax.text(0.05,temp,u'NÂº coletas %d || %s < %s < %s'%(len(s.H),dirI,dirM,dirF))
            temp -= 0.1
        
        ax.text(0.05,temp,u'rho dados brutos = %.3f'%(scipy.stats.pearsonr(s.H,s.T)[0]))
        temp -= 0.06
        
        try:
            ax.text(0.05,temp,u'mh = %.5f'%(s.mh))
            temp -= 0.06
            ax.text(0.05,temp,u'm2h = %.5f'%(s.m2h))
            temp -= 0.06
            ax.text(0.05,temp,u'mt = %.5f'%(s.mt))
            temp -= 0.06
            ax.text(0.05,temp,u'm2t = %.5f'%(s.m2t))
            temp -= 0.06
            ax.text(0.05,temp,u'mht = %.5f'%(s.mht))
            temp -= 0.06
            ax.text(0.05,temp,u'dph = %.5f'%(s.dph))
            temp -= 0.06
            ax.text(0.05,temp,u'dpt = %.5f'%(s.dpt))
            temp -= 0.06
            ax.text(0.05,temp,u'rho funcao ajustada = %.5f'%(s.rhomdc))
            temp -= 0.1    
        
        except:
            ax.text(0.05,temp,u'para rho rode obj.start() antes do obj.contourf()')
            temp -= 0.1    
            pass
        
        ax.text(0.05,temp,u'parametros:')
        temp -= 0.06
        
        if s.fHtype == 'lognormal':
            s.fln_H_fit
            ax.text(0.05,temp,u'fln_H_fit: %s'%s.fln_H_fit)
            temp -= 0.06
            pass
        
        if s.fHtype == 'weibull':
            
            ax.text(0.05,temp,u'fwei_H_fit: %s'%s.fwei_H_fit)
            temp -= 0.06   
            pass
        
        if s.fTtype == 'lognormal':
            ax.text(0.05,temp,u'Tfxi.coef: %s'%s.ln_func_param_xi)
            temp -= 0.06    
            ax.text(0.05,temp,u'Tflamb.coef %s:'%s.ln_func_param_lamb)
            temp -= 0.06   
        pass
        
        if s.fTtype == 'weibull':
            ax.text(0.05,temp,u'Tflambw.coef %s:'%s.Tflambw.coef)
            temp -= 0.06    
            ax.text(0.05,temp,u'Tfalphaw.coef %s:'%s.Tfalphaw.coef)
            temp -= 0.06   
            pass
        
        ax = fig.add_subplot(236)
#        s.ln_func_param_lamb = np.log(s.ln_func_param_lamb)
        if s.fTtype == 'lognormal':
            x = np.linspace(0,H.max(),100)
            ax.plot(x,s.Tfxi(x,*s.ln_func_param_xi),'r', label = 'fxi(h)')
            ax.plot(x,s.Tflamb(x,*s.ln_func_param_lamb),'g',label = 'flambda(h)')
            ax.plot(s.dl.loc[:,'HsBinClass'],s.dl.loc[:,'shape'],'o',color = 'r',label = 'xi')
            ax.plot(s.dl.loc[:,'HsBinClass'],s.dl.loc[:,'scale'],'o',color = 'g',label = 'lambda')
            ax.set_title('fxi(h) e flamb(h)')
            ax.legend()
        else :
            x = np.linspace(0,H.max(),100)
            ax.plot(x,s.Tflambw(x),'b',label = 'lambdaw = flamb(h)')
            ax.plot(x,s.Tfalphaw(x),'y',label = 'alphaw')
            ax.plot(s.dw.iloc[:,3].values,s.dw.iloc[:,0].values,'o',color ='b',label = 'lambda')
            ax.plot(s.dw.iloc[:,3].values,s.dw.iloc[:,2].values,'o',color ='y',label = 'alphas')
            ax.legend()
            
        fig.tight_layout()
        
        return fig

#m = mdc()
#m.fit(H,T,tipofH=aux,tipofT=aux2,bins=CtrlBinClasses,rangeH=CtrlRangeHs,rangeT=CtrlRangeTp,polydegree=CtrlPolyFitGrau,print_cdf=True)

#%%Fim

# =============================================================================
# Classe Nataf
# =============================================================================
   
class nataf():
    def __init__(s):
#        scipy.stats.__init__(self)
#        threading.Thread.__init__(s)
        pass
    
    def __doc__():
        '''
        Ã‰ necessario importar as seguintes bibliotecas
        import scipy
        import pandas as pd
        import numpy as np
        '''
        pass
    
    
    def integral(s):
        ran = [[s.H.min(),s.H.max()],[s.T.min(),s.T.max()]]
        s.i = nquad(lambda h,t: s.pdf(h,t),ran)
        print('calculo de integral terminado = ',s.i)
        return s.i
    
    def fit(s,H,T,**kwargs):
        
        
        #faz o fit para os parametros de hs em distribuicao Ln e Wei
        
        shapeH,lH,scaleH = lognorm.fit(H,floc = 0) 
        s.fln_H = lognorm(shapeH,lH,scaleH) #cria a distribuicao
        s.fln_H_fit = '%.4f,%.4f'%(shapeH,np.log(scaleH))
        
        shapeH,lH,scaleH = weibull_min.fit(H) #fit wei
        s.fwei_H = weibull_min(shapeH,lH,scaleH) #cria a distribuicao
        s.fwei_H_fit = '%.4f,%.4f,%.4f'%(shapeH,lH,scaleH)
###############################################################################        
        #Trocando o fit da biblio weibull_min por wei
#        s.fwei_H = wei()
#        s.fwei_H.fit(H)
#        s.fwei_H_fit = '%.4f,%.4f'%(s.fwei_H.lambW,s.fwei_H.alphaW)
###############################################################################
        
        #faz o fit para os parametros de hs em distribuicao Ln e Wei
        
        shapeT,lT,scaleT = lognorm.fit(T,floc = 0)
        s.fln_T = lognorm(shapeT,lT,scaleT) #cria a ditribuicao ln de tp
        s.fln_T_fit = '%.4f,%.4f'%(shapeT,np.log(scaleT))
        
        shapeT,lT,scaleT = weibull_min.fit(T)
        s.fwei_T = weibull_min(shapeT,lT,scaleT)
        s.fwei_T_fit = '%.4f,%.4f,%.4f'%(shapeT,lT,scaleT)
        
###############################################################################        
#        #Trocando o fit da biblio weibull_min por wei
#        s.fwei_T = wei()
#        s.fwei_T.fit(H)
#        s.fwei_T_fit = '%.4f,%.4f'%(s.fwei_T.lambW,s.fwei_T.alphaW)
###############################################################################
        
        
        s.H = H #guarda a serie de hs internamente
        s.T = T #guarda a serie de tp internamente 
               
        #Usa a distribuicao conforme escolha do usuario
        if kwargs.has_key('tipofH'):
            if kwargs['tipofH'] == 'weibull':
                s.fH = s.fwei_H
                s.fHtype = 'weibull'
            if kwargs['tipofH'] == 'lognormal':
                s.fH = s.fln_H
                s.fHtype = 'lognormal'
                pass
        else: 
            s.fH = s.fln_H
            s.fHtype = 'lognormal'
            pass
        
        
        if kwargs.has_key('tipofT'):
            if kwargs['tipofT'] == 'weibull':
                s.fT = s.fwei_T
                s.fTtype = 'weibull'
            if kwargs['tipofT'] == 'lognormal':
                s.fT = s.fln_T
                s.fTtype = 'lognormal'
        else:
            s.fT = s.fln_T
            s.fTtype = 'lognormal'
        
        
        # Define automaticamente os bins e ranges ou usa os kwargs.
        s.bins = 20
        s.rangeH = [0,10]
        s.rangeT = [0,20]
      
        if kwargs.has_key('bins'):
            s.bins = kwargs['bins']
        if kwargs.has_key('rangeH'):
            s.rangeH = kwargs['rangeH']
            pass
        if kwargs.has_key('rangeT'):
            s.rangeT = kwargs['rangeT']
            pass
                
        #acha a rhobru, correlacao entre H e T para ser usada na pdf
        s.rhobru = scipy.stats.pearsonr(H,T)[0]
#        s.rhobru = scipy.stats.spearmanr(H,T)[0]

        #define uma funcao normal padrao para criar u1 e u2, ou uh e ut
        s.N = norm(0,1)
        
        s.phi_1 = lambda u: 1/(np.sqrt(2*np.pi)) * np.exp(-u**2/2)
       
        s.phi_2 = lambda u1,u2,rhobru: (2*np.pi*np.sqrt(1-rhobru**2))**-1 * np.exp((-2*(1-rhobru**2))**-1 * (u1**2 + u2**2 - 2*rhobru*u1*u2))
        pass
#    
    def rmse(s,X,Y):
        return np.sqrt(((X-Y)**2).mean())
        
    
    def pdf(s,h,t):
        
        uh = s.N.ppf(s.fH.cdf(h))
        if np.isinf(uh):
            uh = 0.9999999999999999
        
        ut = s.N.ppf(s.fT.cdf(t))
        if np.isinf(ut):
            ut = 0.9999999999999999
  
        
        t1 = s.fH.pdf(h) * s.fT.pdf(t) * s.phi_2(uh,ut,s.rhobru)
        
#        print(t1)
        t2 = (s.N.pdf(uh) * s.N.pdf(ut))   
        
#        print(t2)
        return t1/t2
    
    
    def pdf_series(s,H,T):
        
        #parametros do histograma de H
        hs_hist_y,hs_hist_x = np.histogram(H,s.bins,density = True,range = s.rangeH)
        hs_hist_xM = hs_hist_x[:-1] + hs_hist_x[0:2].mean()
        
        #parametro do histograma de T
        tp_hist_y,tp_hist_x = np.histogram(s.T,s.bins,density = True,range = s.rangeT)
        tp_hist_xM = tp_hist_x[:-1] + tp_hist_x[0:2].mean()



        #parametros para plotagem dos dados da tabela df                            
        xpos,ypos = np.meshgrid(hs_hist_xM,tp_hist_xM)
        dfR = pd.DataFrame(np.ndarray((len(xpos),len(ypos))),index = tp_hist_xM, columns = hs_hist_xM)
        dfR[:] = np.nan
        
        for h in hs_hist_xM:
            for t in tp_hist_xM:
                uh = s.N.ppf(s.fH.cdf(h))
                if np.isinf(uh):
                    uh = 0.9999999999999999
                ut = s.N.ppf(s.fT.cdf(t))
                if np.isinf(ut):
                    ut = 0.9999999999999999
                dfR.loc[t,h] = s.fH.pdf(h) * s.fT.pdf(t) * s.phi_2(uh,ut,s.rhobru) / (s.phi_1(uh) * s.phi_1(ut))
                #observar que as posicoes do indice realmente sao [t,h] porque h esta no eixo x
                pass
            pass
        #dfR = dfR.T
        return dict(xpos=xpos,ypos=ypos,zvalues = dfR)
        pass
        
    def rho(s):
        ran = [[s.H.min()*0.8,s.H.max()*1.2],[s.T.min()*0.8,s.T.max()*1.2]]
        
       
        print('------------------------------------------')
        print('calculando media H')
        s.mh = nquad(lambda h,t: h* s.pdf(h,t),ran)[0]
        if s.H.mean() - s.mh > 1 or s.H.mean() - s.mh < -1:
            return 'erro no calculo de mh = %.6f'%s.mh
        print('calculando media H = %.6f'%s.mh)
        
        print('calculando media**2 H')
        s.m2h = nquad(lambda h,t: h** 2 * s.pdf(h,t),ran)[0]
        print('calculando media**2 H = %.6f'%s.m2h)        
    
        print('calculando media T')
        s.mt = nquad(lambda h,t: t * s.pdf(h,t),ran)[0]
        print('calculando media T = %.6f'%s.mt)

        print('calculando media **T')
        s.m2t = nquad(lambda h,t: t ** 2 * s.pdf(h,t),ran)[0]
        print('calculando media **T = %.6f'%s.m2t)
                
        print('calculando E')
        s.mht = nquad(lambda h,t: h * t * s.pdf(h,t),ran)[0]    
        print('calculando mht = %.6f'%s.mht)
        
        print('calculando desvio padrao H')
        s.dph = np.sqrt(s.m2h - s.mh**2) #desvio padrao
        print('calculando desvio padrao H = %.6f'%s.dph)

        print('calculando desvio padrao de T')
        s.dpt = np.sqrt(s.m2t - s.mt**2)    
        print('calculando desvio padrao de T = %.6f'%s.dpt)
        
        s.rhomdc = (s.mht - s.mh*s.mt) / (s.dph * s.dpt)
        print('calculo rho terminado = ',s.rhomdc )
        print('------------------------------------------')
        return s.rhomdc


    def printContour(s,**kwargs): #NATAFFFFF
        
        #determinando serie distribuicao conjunta / nataf
        obj_hist = s.pdf_series(s.H,s.T)
        s.matriz = obj_hist
#        obj_hist['xpos'] = obj_hist['xpos'][::-1]
#        obj_hist['ypos'] = obj_hist['ypos'][::-1]
        
        #determinando serie bruta
        seriebruta, xedges, yedges = np.histogram2d(s.H, s.T, bins = s.bins, range=[s.rangeH,s.rangeT], normed=True)
        seriebruta = seriebruta.T
        xedges = xedges[:-1] + xedges[0:2].mean()
        yedges = yedges[:-1] + yedges[0:2].mean()
        xgrid,ygrid = np.meshgrid(xedges,yedges)
        s.bruta = [seriebruta, xedges, yedges]
        
        
        #Fazendo o x ticks labels
        xtl = []
        for aux in xrange(len(xgrid[0])):
            if aux % 4 == 0:
                xtl.append('%.1f'%xgrid[0][aux])
                pass
            else:
                xtl.append('')
                pass
            pass
            
        #parametros de cores do contour
        cmap2 = cm.get_cmap('jet')
        cmap2.set_under('w')
        cmapbias = cm.get_cmap('seismic')
        levels = np.arange(0.001,0.3,0.01)
        levelsBias = np.arange(-0.11,0.11,0.001)
        
        #iniciando a figura
        fig = plt.figure(figsize=(12,6))
        
        ####################################
        ax = fig.add_subplot(231)
        
        colbar = ax.contour(xgrid,ygrid, seriebruta,cmap=cmap2,vmin = 0.001, vmax = 0.3,levels=levels)
        #ax.set_xscale('log')   
        #ax.set_yscale('log')
        ax.set_xlabel('Hs')
        ax.set_ylabel('Tp')

        ax.set_title('Histograma Densidade Probabilidade\n Serie Bruta')
        
        yticks = obj_hist['zvalues'].index.values.astype(float)
        yticklabels = ['%.1f'%aux for aux in yticks]
        
        ax.set_yticks(yticks[::2])
        ax.set_yticklabels(yticklabels[::2])
        
        
        xticks = obj_hist['zvalues'].columns.values.astype(float)
        xticklabels = ['%.1f'%aux for aux in xticks]
        
        ax.set_xticks(xticks[::2])
        ax.set_xticklabels(xticklabels[::2])
#        ax.set_xlim(0.25,11.7)
        ax.set_ylim(4,20)
        
        ax.grid()
        fig.colorbar(colbar, shrink=0.5, aspect=5,)
        
        #################################
        ax = fig.add_subplot(232)
        colbar = ax.contour(obj_hist['xpos'],obj_hist['ypos'],obj_hist['zvalues'],cmap=cmap2, vmin = 0.001, vmax = 0.3,levels=levels)
        #ax.set_xscale('log')   
        #ax.set_yscale('log')
        ax.set_xlabel('Hs')
        ax.set_ylabel('Tp')

        ax.set_title('Nataf\nHs ajustado em %s\nTp ajustado em %s'%(s.fHtype,s.fTtype))
        
        ax.set_xticks(xticks[::2])
        ax.set_xticklabels(xticklabels[::2])
        ax.set_yticks(yticks[::2])
        ax.set_yticklabels(yticklabels[::2])
#        ax.set_xlim(0.25,11.7)
        ax.set_ylim(4,20)
        
        
        ax.grid()
        fig.colorbar(colbar, shrink=0.5, aspect=5)
        
        #################################
        ax = fig.add_subplot(233)
        bias = seriebruta - obj_hist['zvalues']
        colbar = ax.contour(obj_hist['xpos'],obj_hist['ypos'], bias,cmap=cmapbias, vmin = -0.15, vmax = 0.15,levels=levelsBias)
        #ax.set_xscale('log')   
        ##ax.set_yscale('log')

        ax.set_xlabel('Hs')
        ax.set_ylabel('Tp')
        ax.set_title('Bias Brutos - MDC')

        ax.set_xticks(xticks[::2])
        ax.set_xticklabels(xticklabels[::2])
        ax.set_yticks(yticks[::2])
        ax.set_yticklabels(yticklabels[::2])
#        ax.set_xlim(0.25,11.7)
        ax.set_ylim(4,20)
        

        ax.grid()
        fig.colorbar(colbar, shrink=0.5, aspect=5)   
        
               
        #plotando series unidas
        
        
        ax = fig.add_subplot(234)
        ax.plot(seriebruta.flatten(),label = 'Serie Bruta',ls='--')
        ax.plot(obj_hist['zvalues'].values.flatten(),label = 'Serie Ajustada',ls=':')
        
        xticks = np.linspace(0,s.bins**2,s.bins,endpoint=False)
        temp = obj_hist['zvalues'].index.values.astype(float)
        xticklabels = ['%.1f'%aux for aux in temp]
        
        ax.set_xticks(xticks[::2])
        ax.set_xticklabels(xticklabels[::2])        
        ax.grid()
        ax.set_title(u'Secções de Hs em Função de Tp')
        ax.set_xlabel(u'Secções de Tp')
        ax.set_ylabel('Densidade Probabilidade')
        ax.legend()
        
            
        # =============================================================================
        #  Escrevendo parametros da funcao na figura.       
        # =============================================================================
        
        ax  = fig.add_subplot(235)
        temp = 0.94
        
        if kwargs.has_key('filename'):
            filename = kwargs['filename']
            ax.text(0.05,temp,u'Arquivo: %s'%filename)
            temp -= 0.06
        
        if kwargs.has_key('dirI') and kwargs.has_key('dirM') and kwargs.has_key('dirF'):
            dirI = str(kwargs['dirI'])
            dirM = str(kwargs['dirM'])
            dirF = str(kwargs['dirF'])
            ax.text(0.05,temp,u'NÂº coletas %d || %s < %s < %s'%(len(s.H),dirI,dirM,dirF))
            temp -= 0.1
        
        ax.text(0.05,temp,u'rho dados brutos = %.3f'%(scipy.stats.pearsonr(s.H,s.T)[0]))
        temp -= 0.06
        
        try:
            ax.text(0.05,temp,u'mh = %.5f'%(s.mh))
            temp -= 0.06
            ax.text(0.05,temp,u'm2h = %.5f'%(s.m2h))
            temp -= 0.06
            ax.text(0.05,temp,u'mt = %.5f'%(s.mt))
            temp -= 0.06
            ax.text(0.05,temp,u'm2t = %.5f'%(s.m2t))
            temp -= 0.06
            ax.text(0.05,temp,u'mht = %.5f'%(s.mht))
            temp -= 0.06
            ax.text(0.05,temp,u'dph = %.5f'%(s.dph))
            temp -= 0.06
            ax.text(0.05,temp,u'dpt = %.5f'%(s.dpt))
            temp -= 0.06
            ax.text(0.05,temp,u'rho funcao ajustada = %.5f'%(s.rhomdc))
            temp -= 0.1    
        
        except:
            ax.text(0.05,temp,u'para rho rode obj.start() antes do obj.contourf()')
            temp -= 0.1    
            pass
        
        ax.text(0.05,temp,u'parametros:')
        temp -= 0.06
        
        if s.fHtype == 'lognormal':
            s.fln_H_fit
            ax.text(0.05,temp,u'fln_H_fit: %s'%s.fln_H_fit)
            temp -= 0.06
            pass
        
        if s.fHtype == 'weibull':
            
            ax.text(0.05,temp,u'fwei_H_fit: %s'%s.fwei_H_fit)
            temp -= 0.06   
            pass
        
        if s.fTtype == 'lognormal':
            ax.text(0.05,temp,u'fln_T_fit: %s'%s.fln_T_fit)
            temp -= 0.06    
            
               
        if s.fTtype == 'weibull':
            ax.text(0.05,temp,u'fwei_T_fit %s:'%s.fwei_T_fit)
            temp -= 0.06    
            pass
        
        ax = fig.add_subplot(236)
        
        fig.tight_layout()
        
        return fig
    
#    
#n = nataf()
#n.fit(H,T,tipofH = aux,tipofT = aux2,bins=CtrlBinClasses,rangeH=CtrlRangeHs,rangeT=CtrlRangeTp)
#
#fig = n.printContour(filename = DirFileList[CtrlfileNum],
#                     dirI = direcaofinal[dirInt[dirQuad]['dirI']],
#                     dirM = direcaofinal[dirInt[dirQuad]['dirM']],
#                     dirF = direcaofinal[dirInt[dirQuad]['dirF']])


#%% Main
if __name__ == '__main__':
    
    
    #%% Abrindo o arquivo. Escolha manual no except ou chamado pela linha de comando no try.
    if '__file__' in globals():
        FILE = __file__
    
    DirName = os.path.dirname(os.path.normpath(FILE)) + '/Dados_Hindcast/' 
    DirFigures = os.path.dirname(os.path.normpath(FILE)) + '/Dados_Hindcast/0Figuras/'
     
    proName = os.path.basename(os.path.abspath(FILE)) 
    DirFileList = [fil for fil in os.listdir(DirName) if os.path.isfile(os.path.join(DirName, fil))]
    #%%
    
    CtrlfileNum = 6 #o arquivo que serÃ¡ rodado
    CtrlfileNum = 0 #o arquivo que serÃ¡ rodado
    
    print(DirFileList[CtrlfileNum]) # nome do arquivo que sera analisado
    
#    sys.exit(15)
    
    CtrlDirNum = 1 #numero de direções que serÃ£o divididos os scatter. Os intervalos calculados ficam armazenados em dirInt (direcao intervalo)
#    CtrlDirQuad = [1,2,3,4,5,6,7,8] #Quais as direcoes que serao executadas
    CtrlDirQuad = [1]
    
    #controle de entrada que entrarao nas funcoes
    CtrlBinClasses = 20 #definiÃ§Ã£o do numero de bins no histograma inicial
    CtrlPolyFitGrau = 1 
    CtrlMinDataLen = 50 #quantidade de coletas na direcao para ser analisado
    CtrlRangeHs = [0,8] #Define o range de 0.1 a x para Hs
    CtrlRangeTp = [0,20] #Define o range de 0.1 a x para Tp
#    CtrlDistTypeH = ['lognormal','weibull'] # ou so ['lognormal'] ou sÃ³ ['weibull']
#    CtrlDistTypeT = ['lognormal','weibull'] # ou so ['lognormal'] ou sÃ³ ['weibull']
    CtrlDistTypeH = ['weibull'] # ou so ['lognormal'] ou sÃ³ ['weibull']
    CtrlDistTypeT = ['lognormal'] # ou so ['lognormal'] ou sÃ³ ['weibull']
    CtrlPrintCDF = True
    CtrlCalcRho = True    
    
    auxTT = 0 #faz com que o mapa so seja plotado uma vez
    
    font = {#'family' : 'normal',
            #'weight' : 'bold',
            'size'   : 8}
    plt.rc('font', **font)
        
    #%% Abrindo arquivos           
    temp = 0
    try:
        print 'abrindo arquivo: '
        fileName =  sys.argv[1]
        print 'fileName: ',fileName
        with open(fileName,'r') as f:
            temp = f.read()
            temp.replace('masked','nan')
            f.close()
            d = eval(temp)
            del(temp)
            DirFileList[CtrlfileNum] = os.path.basename(fileName) #necessario para salvar os arquivos com o nome certo. 
            DirFigures = os.path.dirname(fileName) + '/0Figuras/'
            print 'DirFigures: ' + DirFigures
            pass
        pass
    except:
        print 'escolhendo o arquivo manualmente: ', (DirName+DirFileList[CtrlfileNum])
        with open(DirName+DirFileList[CtrlfileNum],'r') as f:
            filename = DirName+DirFileList[CtrlfileNum]
            temp = f.read()
            temp.replace('masked','nan')
            f.close()
            d = eval(temp)
            del(temp)
            pass
        pass
    pass
    
    df = pd.DataFrame(d)
    
    #%%Criando lista de direcoes
    
        
    temp = []
    for aux in range(16):
        temp.append(aux * 22.5)
        pass
    temp.append(360.0)
                
    temp2 = ['N','NNE','NE','ENE','E','ESE','SE','SSE',
                     'S','SSO','SO','OSO','O','ONO','NO','NNO','N']
    direcaofinal = dict(zip(temp,temp2))
    
    dirInt = [] #dicionario que guarda o valor dos intervalos de direcoes. 
    dirInt.append(None)
    
    for ix,aux in enumerate(xrange(CtrlDirNum)):
        dirI = (aux + 1) * 360/CtrlDirNum - 360/(CtrlDirNum*2) #determinando a direcao inicial
        if dirI < 0:
            dirI = dirI + 360
            pass
        dirF = (aux + 1) * 360/CtrlDirNum + 360/CtrlDirNum - 360/(CtrlDirNum*2) #determinando a direcao final
        if dirF > 360:
                dirF = dirF - 360
                pass
        dirM = (aux + 1) * 360/CtrlDirNum
            
        if dirI == dirF:
            dirI = 0 
            dirF = 360
            pass
        dirInt.append(dict(dirI=dirI,dirM = dirM, dirF = dirF))
        pass
    

  

    
#%% Exemplos de montagem

    #Jogando as direcoes para as variaveis H e T
#    if dirInt[CtrlDirQuad]['dirI'] > dirInt[CtrlDirQuad]['dirF']:
#        dHT = df.loc[:,['hs','tp']][np.logical_or(dirInt[CtrlDirQuad]['dirI'] < df.dp , df.dp < dirInt[CtrlDirQuad]['dirF'])]
#    else:
#        dHT = df.loc[:,['hs','tp']][np.logical_and(dirInt[CtrlDirQuad]['dirI'] < df.dp , df.dp < dirInt[CtrlDirQuad]['dirF'])]
#
#    
#    H = dHT.hs
#    T = dHT.tp



#    dirname = os.path.dirname(FILE) +'/Dados_Hindcast/0Figuras'
#    filename =  DirFileList[CtrlfileNum] + '_%s_mdc_.jpg'%CtrlDirQuad
#
#    mll = mdc()
#    mll.fit(H,T,tipofH='lognormal',tipofT='lognormal',bins=20,rangeH=[0,H.max()*1.01],rangeT=[0,T.max()*1.01])
#    mll.rho()
#    
#    fig = mll.printContour(filename = filename,
#                     dirI = direcaofinal[dirInt[CtrlDirQuad]['dirI']],
#                     dirM = direcaofinal[dirInt[CtrlDirQuad]['dirM']],
#                     dirF = direcaofinal[dirInt[CtrlDirQuad]['dirF']])
#                     
#    
#    fig.savefig(os.path.join(dirname,filename))
    
    #%%
    
#    dirname = os.path.dirname(FILE) +'/Dados_Hindcast/0Figuras'
#    filename =  DirFileList[CtrlfileNum] + '_%s_nat_.jpg'%CtrlDirQuad
#    
#    nat = nataf()
#    nat.fit(H,T,tipofT = 'lognormal',tipofH = 'lognormal',bins=20,rangeH=[0,H.max()*1.01],rangeT=[0,T.max()*1.01])
#    nat.rho()
#    fig = nat.printContour(filename = DirFileList[CtrlfileNum],
#                     dirI = direcaofinal[dirInt[CtrlDirQuad]['dirI']],
#                     dirM = direcaofinal[dirInt[CtrlDirQuad]['dirM']],
#                     dirF = direcaofinal[dirInt[CtrlDirQuad]['dirF']])
#    
#    fig.savefig(os.path.join(dirname,filename))
#    print(nat.integral())
    
#%%  Exemplo de montagem em serie

    for aux in CtrlDistTypeH:
        for aux2 in CtrlDistTypeT:
#            for CtrlDirQuad in xrange(1,CtrlDirNum,1):
            for dirQuad in CtrlDirQuad:
                
                #Jogando as direcoes para as variaveis H e T
                if dirInt[dirQuad]['dirI'] > dirInt[dirQuad]['dirF']:
                    dHT = df.loc[:,['hs','tp']][np.logical_or(dirInt[dirQuad]['dirI'] < df.dp , df.dp < dirInt[dirQuad]['dirF'])]
                else:
                    dHT = df.loc[:,['hs','tp']][np.logical_and(dirInt[dirQuad]['dirI'] < df.dp , df.dp < dirInt[dirQuad]['dirF'])]
                    pass
                
                H = dHT.hs
                T = dHT.tp
                
                print('O numero de coletas eh de %d'%len(H))
                if len(H) > CtrlMinDataLen:
                    # =============================================================================
                    # Calculando e salvando mdc
                    # =============================================================================
                    t = time.clock()
                    print('Mdc quadrante %d com H ajustado em %s e T ajustado em %s '%(dirQuad,aux,aux2))
                    
                    
                    dirname = os.path.dirname(FILE) +'/Dados_Hindcast/0Figuras'
                    filename =  DirFileList[CtrlfileNum] + '_%s_mdc_%s_%s.jpg'%(dirQuad,aux,aux2)
                    
                    m = mdc()
#                    m.fit(H,T,tipofH=aux,tipofT=aux2,bins=CtrlBinClasses,rangeH=CtrlRangeHs,rangeT=CtrlRangeTp,polydegree=CtrlPolyFitGrau,print_cdf=True)
                    m.fit(H,T,tipofH=aux,tipofT=aux2,bins=CtrlBinClasses,rangeH=CtrlRangeHs,rangeT=CtrlRangeTp,polydegree=CtrlPolyFitGrau,print_cdf=CtrlPrintCDF)
                    if CtrlCalcRho:
                        m.rho()
                    fig = m.printContour(filename = DirFileList[CtrlfileNum],
                                         dirI = direcaofinal[dirInt[dirQuad]['dirI']],
                                         dirM = direcaofinal[dirInt[dirQuad]['dirM']],
                                         dirF = direcaofinal[dirInt[dirQuad]['dirF']])
                    fig.savefig(os.path.join(dirname,filename))
                    print(time.clock() - t)
                    # =============================================================================
                    # Calculando e salvando Nataf                     
                    # =============================================================================
                    t = time.clock()
                    print('Nataf quadrante %d com H ajustado em %s e T ajustado em %s '%(dirQuad,aux,aux2))
                    print('O numero de coletas eh de %d'%len(H))
                    
                    dirname = os.path.dirname(FILE) +'/Dados_Hindcast/0Figuras'
                    filename =  DirFileList[CtrlfileNum] + '_%s_nat_%s_%s.jpg'%(dirQuad,aux,aux2)

                    n = nataf()
                    n.fit(H,T,tipofH = aux,tipofT = aux2,bins=CtrlBinClasses,rangeH=CtrlRangeHs,rangeT=CtrlRangeTp)
                    
                                       
                    if CtrlCalcRho:
                        n.rho()
                        
                    fig = n.printContour(filename = DirFileList[CtrlfileNum],
                                         dirI = direcaofinal[dirInt[dirQuad]['dirI']],
                                         dirM = direcaofinal[dirInt[dirQuad]['dirM']],
                                         dirF = direcaofinal[dirInt[dirQuad]['dirF']])
                    fig.savefig(os.path.join(dirname,filename))
                    print(time.clock() - t)
                    pass
                pass
            pass
        pass

    
    
    
#    CtrlBinClasses = 28 #definiÃ§Ã£o do numero de bins no histograma inicial
#    CtrlPolyFitGrau = 3 
#    CtrlMinDataLen = 50 #quantidade de coletas na direcao para ser analisado
#    CtrlRangeHs = [0,H.max()*1.01] #Define o range de 0.1 a x para Hs
#    CtrlRangeTp = [0,T.max()*1.01]] #Define o range de 0.1 a x para Tp
#    CtrlDistTypeH = ['lognormal','weibull'] # ou so ['lognormal'] ou sÃ³ ['weibull']
#    CtrlDistTypeT = ['lognormal','weibull'] # ou so ['lognormal'] ou sÃ³ ['weibull']
    #%%
        
#    import multiprocessing as mp
#    
#    def teste(q):
#        t = nquad(lambda x,y: x**2 + x + y**2 + y,[[-np.inf,np.inf],[-np.inf,np.inf]])
#        q.put(t)
#        
#    q = mp.Queue()
#    p = mp.Process(target = teste,args = (q,)) 
#    p.start()    
#    print (q.get())
#    p.join()

#%%
    
    
#    from multiprocessing import Process, Pipe
#    
#    def f(conn):
#        conn.send([42, None, 'hello'])
#        conn.close()
#    
#    if __name__ == '__main__':
#        parent_conn, child_conn = Pipe()
#        p = Process(target=f, args=(child_conn,))
#        p.start()
#        print parent_conn.recv()   # prints "[42, None, 'hello']"
#        p.join()
#    
    
    #%%
#    d = m.matriz
#    d1 = n.bruta
    
    sys.exit(10)
    
#%% EXIT
    
m.dl.loc[:,'shape'].values
plt.plot(m.dl.loc[:,'HsBinClass'].values,m.dl.loc[:,'shape'].values,'o')
plt.plot(m.dl.loc[:,'HsBinClass'].values,m.dl.loc[:,'scale'].values,'o')

from scipy.optimize import curve_fit
xdata = m.dl.loc[:,'HsBinClass'].values
ydata = m.dl.loc[:,'shape'].values

f = lambda h,a1,a2,b1,b2: a1 * np.exp(-b1*h) + a2 * np.exp(-b2*h)

cfp,temp = curve_fit(f,xdata,ydata)

x = np.linspace(1,5,100)
plt.plot(x,[f(aux,*cfp) for aux in x])

f = lambda h,a1,a2,b: a1+ a2* h ** b

xdata = m.dl.loc[:,'HsBinClass'].values
ydata = m.dl.loc[:,'scale'].values
cfp,temp = curve_fit(f,xdata,ydata)

x = np.linspace(1,5,100)
plt.plot(x,[f(aux,*cfp) for aux in x])

#%%

cfp

f(2,cfp[0],cfp[1],cfp[2],cfp[3])

f(2,*cfp)

cfp

np.log(np.e)
np.exp(np.e)
e = np.e
e**e
np.log10(10)

#%%

plt.figure()
plt.plot(x,[np.exp(aux) for aux in x])
plt.plot(x,[np.log(aux) for aux in x])
plt.semilogy()
plt.figure()
plt.plot(x,[np.exp(aux) for aux in x])
plt.plot(x,[np.log(aux) for aux in x])
plt.loglog()
#%%
    
from scipy.stats import norm,weibull_min
from scipy.special import gamma
from scipy.optimize import fsolve
from scipy.integrate import quad
import numpy as np
from scipy.optimize import root

class ln():
   
    def __doc__(s):
        print(''' Must import these libraries
        from scipy.integrate import quad
        import numpy as np              
        ''')
    
    def __init__ (s,*args):
        if(len(args)):
            s.fit(args)
            pass
        pass
    
    def fit(s,X):
        s.xi = s.xiLn(X)
        s.lamb = s.lambLn(X)
        return s.xi,s.lamb
    
    def xiLn(s,X):
        return np.sqrt(np.log(1 + (X.std() /  X.mean())**2)) # esse Ã© o calculo do parametro xi
            
    def lambLn(s,X):
        return np.log(X.mean()) - 0.5*s.xi**2 #esse Ã© o calculo do parametro lambda
            
    def pdf(s,x): #funcao lognormal distribuiÃ§Ã£o
        return (1/(np.sqrt(2*np.pi)*x*s.xi))*np.exp(-0.5*(((np.log(x)-s.lamb)**2)/(s.xi**2)))
         
    def cdf(s,x):
        return quad(lambda x: s.pdf(x),0,x)[0]
    pass


class wei():
       
    def __doc__(s):
        print(''' Must import these libraries
        from scipy.special import gamma
        from scipy.optimize import fsolve
        from scipy.integrate import quad,nquad
        import numpy as np              
        ''')
        
    def __init__ (s,*args):
        if(len(args)):
            s.fit(args)
            pass
            
    def fit(s,X):
        s.lambW = s.lambWei(X)
        s.alphaW = s.alphaWei(X)
        return s.lambW,s.alphaW
   
#    def fitex(s,X):
#        s.lambW = s.lambWei(X)
#        s.alphaW = s.alphaWei(X)
#        return s.alphaW,s.lambW
    
    def lambWei(s,X):
        cofVar = np.nanstd(X)/np.nanmean(X) 
        lambW = fsolve(lambda y: (np.sqrt(gamma((2/y) + 1) - gamma((1/y) + 1)**2) / gamma((1/y)+1)) - cofVar ,3)
        return lambW
    
    def alphaWei(s,X):
        cofVar = np.nanstd(X)/np.nanmean(X) 
        lambW = fsolve(lambda y: (np.sqrt(gamma((2/y) + 1) - gamma((1/y) + 1)**2) / gamma((1/y)+1)) - cofVar ,3)
        return np.nanmean(X)/gamma(1/lambW + 1)

    def pdf(s,x):
        return ((x**(s.lambW-1)) / (s.alphaW**s.lambW)) * s.lambW * np.exp(-(x/s.alphaW)**s.lambW)

    def cdf(s,x):
        return quad(lambda x: s.pdf(x),0,x)[0]

#%%
                   
a,b,c = weibull_min.fit(dHT.hs,floc = 0)
weifloc = weibull_min(a,b,c)

a,b,c = weibull_min.fit(dHT.hs)
weinorm = weibull_min(a,b,c)

wei2 = wei()
wei2.fit(dHT.hs)


x = np.linspace(0,8,100,endpoint=False)

yl = weifloc.cdf(x)
yn = weinorm.cdf(x)    
y2 = [wei2.cdf(aux) for aux in x]



#%% Criando valores da distribuiÃ§Ã£o ajustada para comparar com valores do scatter dos valores brutos. 


z,x = np.histogram(H,bins=CtrlBinClasses*4,range=CtrlRangeHs)
xm = x[:-1] + (x[1] - x[0])/2
z,y = np.histogram(T,bins=CtrlBinClasses*4,range=CtrlRangeTp)
ym = y[:-1] + (y[1] - y[0])/2


t0 = time.clock()
grid = [[m.pdf(ahs,atp) for ahs in xm] for atp in ym]
df2 = pd.DataFrame(data=grid,columns=xm,index=ym)
print(time.clock() - t0)

#t0 = time.clock()
#
#df = np.zeros([len(xm),len(ym)])
#df = pd.DataFrame(data=df,columns=xm,index=ym)
#for ahs in xm:
#    for atp in ym:
#        df.loc[atp,ahs] = m.pdf(ahs,atp)
#        pass
#
#print(time.clock() - t0)

#df == df2
np.log(m.dl.loc[:,"scale"])
