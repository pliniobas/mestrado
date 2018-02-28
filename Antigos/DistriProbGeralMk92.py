# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 23:11:24 2016

@author: Plinio Bueno Andrade Silva
"""

from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm, gridspec
from mpl_toolkits.basemap import Basemap
import numpy as np 
import sys
from scipy.stats import (kurtosis,skew,norm,lognorm,
gumbel_r,gumbel_l,dweibull,weibull_max,weibull_min,frechet_r,spearmanr,rayleigh,gamma)

from scipy.special import gamma
from scipy.optimize import root,fsolve
from scipy.interpolate import interp1d,interp2d
from scipy.integrate import quad,dblquad,nquad
import scipy
#import pandas as pd
#from sympy import *
from math import sqrt,log,exp
import time
from numba import jit
import os
from fractions import gcd
import gc
import pandas as pd


t0 = time.clock()
#print (time.clock() - t0)

# As variaveis que começam con Ctrl são variavies de controle do programa. Ajuste os valores conforme vontade ou necessario.
###Ctrl


#%% Abrindo o arquivo. Escolha manual no except ou chamado pela linha de comando no try.
if '__file__' in globals():
    FILE = __file__

DirName = os.path.dirname(os.path.normpath(FILE)) + '/Dados_Hindcast/' 
DirFigures = os.path.dirname(os.path.normpath(FILE)) + '/Dados_Hindcast/0Figuras/'
 
proName = os.path.basename(os.path.abspath(FILE)) 
DirFileList = [fil for fil in os.listdir(DirName) if os.path.isfile(os.path.join(DirName, fil))]
#%%

CtrlfileNum = 2 #o arquivo que será rodado
DirFileList[CtrlfileNum]


CtrlBinClasses = 28 #definição do numero de bins no histograma inicial
CtrlRangeHs = 10 #Define o range de 0.1 a x para Hs
CtrlRangeTp = 20 #Define o range de 0.1 a x para Tp
CtrlBinCmult = 1 #multiplica a definição, somente para exibição
CtrlDirNum = 4 #numero de direções que serão divididos os scatter

#CtrlLimIntegHs = 10 #Limite do valor de Hs na integral para achar o valor E
#CtrlLimIntegTp = 25 #Limite do valor de Tp na integral para achar o valor E
CtrlLimMinColetasTpHs = 10 #Na função TpHs, controla quantos valores precisa ter....
CtrlPolyFitGrau = 3 #
CtrlMinDataLen = 200 #quantidade de coletas na direcao para ser analisado

CtrlLinhas = int(round(np.sqrt(CtrlBinClasses)))
CtrlColunas = int(np.sqrt(CtrlBinClasses))+1 #linhas para gráifco TpHs. 

Print1 = False
PrintHist3D = False
PrintTpHs = False
PrintXiLamb = False
PrintTpHs3D = False
PrintSurface = False
PrintPcolor = False
PrintContourf = True
PrintRho = True
PrintI = False

auxTT = 0 #faz com que o mapa so seja plotado uma vez


direcoes = []
for aux in range(16):
    direcoes.append(aux * 22.5)
    pass
            
direcoesnomes = ['N','NNE','NE','ENE','E','ESE','SE','SSE',
                 'S','SSO','SO','OSO','O','ONO','NO','NNO']
direcaofinal = dict(zip(direcoes,direcoesnomes))



font = {#'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 9}
plt.rc('font', **font)



#%%
#del(DirFileList[-1])

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


#==============================================================================
#%% Criando figuras e tabelas para analise do periodo em grafico 2D.
#==============================================================================
if False:
    
    if '__file__' in globals():
        FILE = __file__ 
    #Apd = [pd.DataFrame.from_dict(aux) for aux in d]
    Apd = pd.DataFrame.from_dict(d)
    Apd.index = pd.to_datetime(Apd['date'])
    
    
    ASerie = Apd.iloc[:].loc[:,['hs','dp','tp','lat','lon']]
    
    
    ASerie = ASerie.loc[pd.datetime.strptime('12-09-2014','%d-%m-%Y'):pd.datetime.strptime('22-11-2017','%d-%m-%Y')]
    ASerie.to_html(DirFigures + str(ASerie.iloc[0].loc['lat']) + '_' + str(ASerie.iloc[0].loc['lon']) + '_tabela.html')
    ASerie.to_csv(DirFigures + str(ASerie.iloc[0].loc['lat']) + '_' + str(ASerie.iloc[0].loc['lon']) + '_tabela.csv')
    
    xticks_labels_str = ASerie.index.strftime('%Y-%m-%d')
    
    
    Fig = plt.figure()
    sub1 = Fig.add_subplot(311)
    sub1.plot(ASerie.index.values,ASerie.loc[:,'hs'].values)
    sub1.legend('H')
    sub1.set_xticks(ASerie.index[::8])
    sub1.set_xticklabels(xticks_labels_str[::8],animated=True,rotation=70,size='x-small')
    sub1.grid()
    
    sub2 = Fig.add_subplot(312)
    sub2.plot(ASerie.index.values,ASerie.loc[:,'tp'].values)
    sub2.legend('T')
    sub2.set_xticks(ASerie.index[::8])
    sub2.set_xticklabels(xticks_labels_str[::8],animated=True,rotation=70,size='x-small')
    sub2.grid()
    
    sub3 = Fig.add_subplot(313)
    sub3.plot(ASerie.index.values,ASerie.loc[:,'dp'].values)
    sub3.legend('D')
    sub3.set_xticks(ASerie.index[::8])
    sub3.set_xticklabels(xticks_labels_str[::8],animated=True,rotation=70,size='x-small')
    sub3.grid()
    
    Fig.tight_layout()
    #[aux for aux in dir(sub1) if aux.find('gend') > 0]
    Fig.get_size_inches()
    Fig.set_size_inches(15,5)
    plt.show(Fig)
    Fig.savefig(os.path.dirname(os.path.realpath(FILE)) + '/Dados_Hindcast/0Figuras/'
                +str(ASerie.iloc[0].loc['lat']) + '_' + str(ASerie.iloc[0].loc['lon']) + '_figura.png')
        

lat = np.asarray(d['lat'][0])
lon = np.asarray(d['lon'][0])
latint = np.asarray(d['latint'][0])
lonint = np.asarray(d['lonint'][0])
hs = np.asarray(d['hs'])
tp = np.asarray(d['tp'])
dp = np.asarray(d['dp'])
date = np.asarray(d['date'])

dff = pd.DataFrame(d)
del(d)

fig = []
ax = []


#pd.period_range(ASerie.index[0],ASerie.index[-1], freq='1D')
#pd.date_range(ASerie.index[0],ASerie.index[-1], freq='30min')
#pd.date_range(ASerie.index[0],periods=30,freq='1D')

#==============================================================================
# 
#==============================================================================


#%%

#==============================================================================
# Definindo funções necessárias.
#==============================================================================

#onde 'a' é o array que contem os indices do polinomio e x é um valor dentro série de hs (ou X)
def pol(a,x): #Esse é a função polinomial escrita no formato de somatório \sigma a*x^n com i = 0 a i = n.
    ii = range(len(a)) #o tamanho da lista da o numero de indices [a1,a2,a3]
    ii.reverse() #inverte a ordem da lista para casar a potencia maior com o indice maior.
    soma = 0 #(armazena a somatoria)
    aux = 0
    for i in ii:
        soma += a[aux]*x**i
        aux += 1
    return soma


#Achando lambda

#Nessa função se entra com a série de dados e ela retorna o parametro xi
#para ajustar a distribuição com o método dos mometos

### Funções Lognormal xi e lambda
#@jit
def xiLn(npserie):
    mean = npserie.mean() #media
    std = npserie.std() # desvio padrao
    cofVar = std / mean
    return sqrt(log(1 + cofVar**2)) # esse é o calculo do parametro xi

#@jit
def lambLn(npserie):
    mean = npserie.mean() #media
    std = npserie.std() # desvio padrao
    cofVar = std / mean
    xi = sqrt(log(1 + cofVar**2))
    return log(mean) - 0.5*xi**2 #esse é o calculo do parametro lambda


### Funções Lognormal PDFs

#@jit
def flnPDFHs(aux,xi,lamb): #funcao lognormal distribuição
    return (1/(sqrt(2*np.pi)*aux*xi))*exp(-0.5*(((log(aux)-lamb)**2)/(xi**2)))
    pass

#funcao lognormal condicional. . Depende da função pol para funcionar.
#fXi e fLamb são os indices (a,b,c) de um funcao polinomial ax**2 + bx + c. Podem ser obtidos pelo retorno de uma polyfit
#Ex obtendo fXi: iXiHs = np.polyfit(HsList,TpHsCondList,CtrlPolyFitGrau) #indices da função que retorna valores de xi (de Hs)
#@jit
def flnPDFTpHs(TpAux,HsAux,fXi,fLamb):
    return (1/(sqrt(2*np.pi)*TpAux*pol(fXi,HsAux)))*exp(-0.5*((((log(TpAux)-pol(fLamb,HsAux))/pol(fXi,HsAux))**2)))

#@jit
def flnPDFHs_flnTpHs(TpAux,HsAux,xi,lamb,fXi,fLamb): #funcao conjunta de probabilidade (distribuição de Hs * distruição condicional de TpHs)
    return (flnPDFHs(HsAux,xi,lamb)) * (flnPDFTpHs(TpAux,HsAux,fXi,fLamb))  

#@jit
#def f_HsTpXHT(TpAux,HsAux,xi,lamb,fXi,Flamb): #Função de distribuição conjunta para achar \rho
#    return (TpAux * HsAux) * flnPDFHs_flnTpHs(TpAux,HsAux,xi,lamb,fXi,Flamb)


### Funções Lognormal CDFs
def flnCDFHs(aux,xi,lamb):
    return quad(lambda aux: flnPDFHs(aux,xi,lamb),0,aux)



#==============================================================================
### Funções Weibull
#==============================================================================

def alphaw_lambw(npserie): #onde npserie é a série de Hs
    mean = np.nanmean(npserie)    
    std = np.nanstd(npserie)
    cofVar = std/mean
    lamb = fsolve(lambda y: (sqrt(gamma((2/y) + 1) - gamma((1/y) + 1)**2) / gamma((1/y)+1)) - cofVar ,3)
    alpha = mean/gamma(1/lamb + 1)
    return (alpha),(lamb)

#@jit
def weiPDFHs(HsAux,alphaw,lambw):
    return ((HsAux**(lambw-1)) / (alphaw**lambw)) * lambw * exp(-(HsAux/alphaw)**lambw)

#@jit
def weiPDFHs_flnTpHs(TpAux,HsAux,alphaw,lambw,fXi,fLamb): #funcao conjunta de probabilidade (distribuição de Hs * distruição condicional de TpHs)
    return ((flnPDFTpHs(TpAux,HsAux,fXi,fLamb)) * (weiPDFHs(HsAux,alphaw,lambw)))[0]


# =============================================================================
# Classe mdc - metodo da distribuicao condicional
# =============================================================================

class mdc:
    def __init__(s):
        pass
    
    def fit(s,H,T,**kwargs):
        #Criar distribuicao de H Lognormal
        Hscale,Hloc,Hshape = lognorm.fit(H,floc = 0) # fit wei
        s.fln_H = lognorm(Hscale,Hloc,Hshape) #cria a distribuicao
        s.fln_H_fit = '%.4f,%.4f'%(Hscale,np.log(Hshape))
        
        #Criar distribuicao de H Weibull
        Hscale,Hloc,Hshape = weibull_min.fit(H,floc = 0) #fit wei
        s.fwei_H = weibull_min(Hscale,Hloc,Hshape) #cria a distribuicao
        s.fwei_H_fit = '%.4f,%.4f'%(Hscale,Hshape)

        
        #Determinar o hitograma
        s.bins = 20
        s.rangeH = (0,10)
        s.rangeT = (0,20)
        
        if kwargs.has_key('bins'):
            s.bins = kwargs['bins']
        if kwargs.has_key('rangehs'):
            s.rangeH = (0,kwargs['rangehs'])
            pass
        if kwargs.has_key('rangetp'):
            s.rangeT = (0,kwargs['rangetp'])
            pass

        s.H_hist_y,s.H_hist_x = np.histogram(H,s.bins,density = True,range = s.rangeH)
        s.H_hist_xM = s.H_hist_x[:-1] + s.H_hist_x[0:2].mean()
        
        s.T_hist_y,s.T_hist_x = np.histogram(T,s.bins,density = True,range = s.rangeT)
        s.T_hist_xM = s.T_hist_x[:-1] + s.T_hist_x[0:2].mean()
        
        #Separando T condicional a H e calculando os parametros da distribuicao
        dft = pd.DataFrame(dict(H=H,T=T))
        
        ln_param = []
        wei_param = []        
        for ix,aux in enumerate(s.H_hist_x):
            if ix == len(s.H_hist_x) - 1:
                break
            temp = dft[np.logical_and(dft.H>s.H_hist_x[ix],dft.H<s.H_hist_x[ix+1])]['T'].values
            if len(temp) > 50:
                #lista contem [xi,loc,lamb*e,posicaoX,hs_condicionador]
                #xi = [0], loc=[1], lamb =[2], xpos = [3]
                ln_param.append(lognorm.fit(temp,floc=0) + tuple([s.H_hist_xM[ix]]))
                #a lista contem [lamb,loc,alpha,hs_condicionador]
                #lambw=[0],loc=[1], alpha=[2], xpos=[3]
                wei_param.append(weibull_min.fit(temp,floc=0) + tuple([s.H_hist_xM[ix]]))
            pass
        
        
        #Criando as funcoes dos parametros de distribuicao
        s.Tfxi = np.poly1d(np.polyfit([aux[3] for aux in ln_param],[aux[0] for aux in ln_param],1))
        s.Tflamb = np.poly1d(np.polyfit([aux[3] for aux in ln_param],[aux[2] for aux in ln_param],1))

        s.Tflambw = np.poly1d(np.polyfit([aux[3] for aux in wei_param],[aux[0] for aux in wei_param],1))
        s.Tfalphaw = np.poly1d(np.polyfit([aux[3] for aux in wei_param],[aux[2] for aux in wei_param],1))
        
#        x = np.linspace(0,5,100)
#        plt.plot(x,fxi(x))
#        plt.plot(x,flamb(x))
#        plt.plot(x,flambw(x))
#        plt.plot(x,alphaw(x))
        
        
    def fln_TH(s,h): #funcao condicional Lognormal
        TLn = lognorm(s.Tfxi(h),0,s.Tflamb(h))
        return TLn
    
    def fwei_TH(s,h): #funcao condicional Weibull
        Twei = weibull_min(s.Tflambw(h),0,s.Tfalphaw(h))
        return Twei
    
    def fconjunta_ln(s,h,t):
        return (s.HLn(h)) * s.fcondicional_ln(h,t)
    #%%
H=hs
T=tp
teste = mdc()
teste.fit(H,T)       

nquad(lambda h,t: teste.fln_H.pdf(h) * teste.fln_TH(h).pdf(t),[[0.1,8],[0.1,16]])
nquad(lambda h,t: teste.fwei_H.pdf(h) * teste.fwei_TH(h).pdf(t),[[0.1,8],[0.1,16]])

#%% 
ran = [[0.1,H.max()*2],[0.1,T.max()*2]]

mh = nquad(lambda h,t: h*teste.fwei_H.pdf(h) * teste.fwei_TH(h).pdf(t),ran)[0]
m2h = nquad(lambda h,t: h**2 * teste.fwei_H.pdf(h) * teste.fwei_TH(h).pdf(t),ran)[0]
dph = np.sqrt(m2h - mh**2) #desvio padrao
#tp
mt = nquad(lambda h,t: t*teste.fwei_H.pdf(h) * teste.fwei_TH(h).pdf(t),ran)[0]
m2t = nquad(lambda h,t: t**2 *teste.fwei_H.pdf(h) * teste.fwei_TH(h).pdf(t),ran)[0]
dpt = np.sqrt(m2t - mt**2)    
#
mht = nquad(lambda h,t: h * t *teste.fwei_H.pdf(h) * teste.fwei_TH(h).pdf(t),ran)[0]       
rhowei = (mht - mh*mt) / (dph * dpt)

#%%

mh = nquad(lambda h,t: h*teste.fln_H.pdf(h) * teste.fln_TH(h).pdf(t),ran)[0]
m2h = nquad(lambda h,t: h**2 * teste.fln_H.pdf(h) * teste.fln_TH(h).pdf(t),ran)[0]
dph = np.sqrt(m2h - mh**2) #desvio padrao
#tp
mt = nquad(lambda h,t: t*teste.fln_H.pdf(h) * teste.fln_TH(h).pdf(t),ran)[0]
m2t = nquad(lambda h,t: t**2 *teste.fln_H.pdf(h) * teste.fln_TH(h).pdf(t),ran)[0]
dpt = np.sqrt(m2t - mt**2)    
#
mht = nquad(lambda h,t: h * t *teste.fln_H.pdf(h) * teste.fln_TH(h).pdf(t),ran)[0]       
rholn = (mht - mh*mt) / (dph * dpt)
        
#%%

# =============================================================================
# Classe Nataf
# =============================================================================
class nataf():
    def __init__(s):
#        scipy.stats.__init__(self)
        pass
    
    def __doc__():
        '''
        É necessario importar as seguintes bibliotecas
        import scipy
        import pandas as pd
        import numpy as np
        '''
        pass
    
    def fit(s,H,T,**kwargs):
        
        #faz o fit para os parametros de hs em distribuicao Ln e Wei
        
        scaleH,lH,shapeH = lognorm.fit(H,floc = 0) # fit wei
        s.HLn = lognorm(scaleH,lH,shapeH) #cria a distribuicao
        s.HLnfit = '%.4f,%.4f'%(scaleH,np.log(shapeH))
        
        scaleH,lH,shapeH = weibull_min.fit(H,floc = 0) #fit wei
        s.HWei = weibull_min(scaleH,lH,shapeH) #cria a distribuicao
        s.HWeifit = '%.4f,%.4f'%(scaleH,shapeH)
        
        #faz o fit para os parametros de hs em distribuicao Ln e Wei
        
        scaleT,lT,shapeT = lognorm.fit(T,floc = 0)
        s.TLn = lognorm(scaleT,lT,shapeT) #cria a ditribuicao ln de tp
        s.TLnfit = '%.4f,%.4f'%(scaleT,np.log(shapeT))
        
        scaleT,lT,shapeT = weibull_min.fit(T,floc = 0)
        s.TWei = weibull_min(scaleT,lT,shapeT)
        s.TWeifit = '%.4f,%.4f'%(scaleT,shapeT)
        
        s.H = H #guarda a serie de hs internamente
        s.T = T #guarda a serie de tp internamente 
               
        # =============================================================================
        #         
        # =============================================================================
               
#        s.HLn = ln()
#        s.HLn.fit(H)
#        
#        s.TLn = ln()
#        s.TLn.fit(T)
#                
#        s.HWei = wei()
#        s.HWei.fit(H)
#        
#        s.TWei = wei()
#        s.TWei.fit(T)
       
        
        # =============================================================================
        #         
        # =============================================================================
        #Comparacao entre series com RMSE para ver se vai usar Weibull ou Lognormal
        
        if kwargs.has_key('tipofH'):
            if kwargs['tipofH'] == 'weibull':
                s.fH = s.HWei
                s.Htype = 'weibull'
            if kwargs['tipofH'] == 'lognormal':
                s.fH = s.HLn
                s.Htype = 'lognormal'
                pass
        else: 
            s.fH = s.HLn
            s.Htype = 'lognormal'
            pass
        
        
        if kwargs.has_key('tipofT'):
            if kwargs['tipofT'] == 'weibull':
                s.fT = s.TWei
                s.Ttype = 'weibull'
            if kwargs['tipofT'] == 'lognormal':
                s.fT = s.TLn
                s.Ttype = 'lognormal'
        else:
            s.fT = s.TLn
            s.Ttype = 'lognormal'
        
#        print('HLn',s.rmse(H,s.HLn.pdf(r)),'HWei',s.rmse(H,s.HWei.pdf(r)))
#        print('TLn',s.rmse(T,s.TLn.pdf(r)),'TWei',s.rmse(T,s.TWei.pdf(r)))
        
        # Define automaticamente os bins e ranges ou usa os kwargs.
        s.bins = (0,14)
        s.rangehs = (0,10)
        s.rangetp = (0,20)
      
        if kwargs.has_key('bins'):
            s.bins = kwargs['bins']
        if kwargs.has_key('rangehs'):
            s.rangehs = (0,kwargs['rangehs'])
            pass
        if kwargs.has_key('rangetp'):
            s.rangetp = (0,kwargs['rangetp'])
            pass
                
        #acha a rho, correlacao entre H e T para ser usada na pdf
        s.rho = scipy.stats.pearsonr(H,T)[0]
#        s.rho = scipy.stats.spearmanr(H,T)[0]
#        print('s.rho: ', s.rho)

        #define uma funcao normal padrao para criar u1 e u2, ou uh e ut
        s.N = norm(0,1)
        
        s.phi_1 = lambda u: 1/(np.sqrt(2*np.pi)) * np.exp(-u**2/2)
       
        s.phi_2 = lambda u1,u2,rho: (2*np.pi*np.sqrt(1-rho**2))**-1 * np.exp((-2*(1-rho**2))**-1 * (u1**2 + u2**2 - 2*rho*u1*u2))
        pass
#    
    def rmse(s,X,Y):
        return np.sqrt(((X-Y)**2).mean())
        
    
    def pdf(s,h,t):
        
        uh = s.N.ppf(s.fH.cdf(h))
        
        ut = s.N.ppf(s.fT.cdf(t))

        
        t1 = s.fH.pdf(h) * s.fT.pdf(t) * s.phi_2(uh,ut,s.rho)
        t2 = (s.N.pdf(uh) * s.N.pdf(ut))   

        return t1/t2
    
    def pdf_series(s,H,T):
        
        #parametros do histograma de H
        hs_hist_y,hs_hist_x = np.histogram(H,s.bins,density = True,range = s.rangehs)
        hs_hist_xM = hs_hist_x[:-1] + hs_hist_x[0:2].mean()
        
        #parametro do histograma de T
        tp_hist_y,tp_hist_x = np.histogram(T,s.bins,density = True,range = s.rangetp)
        tp_hist_xM = tp_hist_x[:-1] + tp_hist_x[0:2].mean()

        #parametros para plotagem dos dados da tabela df                            
        xpos,ypos = np.meshgrid(hs_hist_xM,tp_hist_xM)
        dfR = pd.DataFrame(np.ndarray((len(xpos),len(ypos))),index = tp_hist_xM, columns = hs_hist_xM)
        dfR[:] = np.nan
        
        for h in hs_hist_xM:
            for t in tp_hist_xM:
                uh = s.N.ppf(s.fH.cdf(h))
                ut = s.N.ppf(s.fT.cdf(t))
                dfR.loc[t,h] = s.fH.pdf(h) * s.fT.pdf(t) * s.phi_2(uh,ut,s.rho) / (s.phi_1(uh) * s.phi_1(ut))
                #observar que as posicoes do indice realmente sao [t,h] porque h esta no eixo x
                pass
            pass
        return dict(xpos=xpos,ypos=ypos,df = dfR)
        pass


#%%

def rho_calc(func,Hmin,Hmax,Tmin,Tmax):
    
    ran = [[Hmin,Hmax],[Tmin,Tmax]]
    print('calculando media H')
    mh = nquad(lambda h,t: h* func.pdf(h,t),ran)[0] #media
    print('calculando media**2 H')
    m2h = nquad(lambda h,t: h**2 * func.pdf(h,t),ran)[0]
    print('calculando desvio padrao H')
    dph = np.sqrt(m2h - mh**2) #desvio padrao
    #Tp
    print('calculando media T')
    mt = nquad(lambda h,t: t * func.pdf(h,t),ran)[0]
    print('calculando media **T')
    m2t = nquad(lambda h,t: t ** 2 * func.pdf(h,t),ran)[0]
    print('calculando desvio padrao de T')
    dpt = np.sqrt(m2t - mt**2)    
    #
    print('calculando E')
    mht = nquad(lambda h,t: h * t * func.pdf(h,t),ran)[0]    
    rhonat = (mht - mh*mt) / (dph * dpt)
    return rhonat

#%%

#Parece ter sido calculado corretamente
natlw = nataf()
natlw.fit(H,T,tipofH='lognormal',tipofT='weibull')
print(natlw.Htype)
print(natlw.Ttype)
rholw = rho_calc(natlw,0.1,8,0.1,16)    

#%% Dando problema acho que eh o H Max do range o problema. ok para (natwl,0.1,7,0.1,20)
natwl = nataf()
natwl.fit(H,T,tipofH='weibull',tipofT='lognormal')
print(natwl.Htype)
print(natwl.Ttype)
rhowl = rho_calc(natwl,0.1,7,0.1,20)  
    
#%%
natww = nataf()
natww.fit(H,T,tipofH='weibull',tipofT='weibull')
print(natww.Htype)
print(natww.Ttype)
rhoww = rho_calc(natww,0.1,8,0.1,16)  
    
#%%

natll = nataf()
natll.fit(H,T,tipofH='weibull',tipofT='weibull')
print(natll.Htype)
print(natll.Ttype)
rholl = rho_calc(natll,0.1,8,0.1,16)

f = open('resultados de rho','a')
f.write('''Para dir = %s 
        natlw %0.6f
        natwl %0.6f
        natww %0.6f
        natll %0.6f
        '''
        ,%(auxT,natlw,natwl,natww,natll))
f.close()


#%%
#from scipy.stats import norm,weibull_min
#from scipy.special import gamma
#from scipy.optimize import fsolve
from scipy.special import gamma
from scipy.optimize import fsolve
from scipy.integrate import quad
import numpy as np


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
        return np.sqrt(np.log(1 + (X.std() /  X.mean())**2)) # esse é o calculo do parametro xi
            
    def lambLn(s,X):
        return np.log(X.mean()) - 0.5*s.xi**2 #esse é o calculo do parametro lambda
            
    def pdf(s,x): #funcao lognormal distribuição
        return (1/(np.sqrt(2*np.pi)*x*s.xi))*np.exp(-0.5*(((np.log(x)-s.lamb)**2)/(s.xi**2)))
         
    def cdf(s,x):
        return quad(lambda x: s.pdf(x),0,x)[0]
    pass


class wei:
    
#    from scipy.stats import norm,weibull_min
#    from scipy.special import gamma
#    from scipy.optimize import fsolve
#    import numpy as np
#    from scipy.optimize import root
    
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
    
    
    

    
    
#==============================================================================
#%% Separando os dados por direção 
#==============================================================================
#hsd = [[]] *17 --> Isso está errado. Vai criar uma copia das sublistas em vez de sublistas independentes
#hsd = [[ ] for aux in xrange(CtrlDirNum+1)]
d = []
temp = 0
for aux2 in range(0,CtrlDirNum):
    hsTemp = []
    tpTemp = []
    dpTemp = []
    dateTemp = []
    
    for aux in range(len(dp)):
        
        di = aux2 * 360/CtrlDirNum - 360/(CtrlDirNum*2) #determinando a direcao inicial
        if di < 0:
            di = di + 360
            pass
        df = aux2 * 360/CtrlDirNum + 360/CtrlDirNum - 360/(CtrlDirNum*2) #determinando a direcao final
        if df > 360:
                df = df - 360
                pass
        dm = aux2 * 360/CtrlDirNum
            
        if di == df:
            di = 0 
            df = 359
            pass
        
        if aux == 0: # Quando aux == 0, di == 338 e df == 27, por isso o or
            
            if di <= dp[aux] or dp[aux] < df: 
                temp+=1
                hsTemp.append((hs[aux]))
                tpTemp.append((tp[aux]))
                dpTemp.append((dp[aux]))
                dateTemp.append((date[aux]))
                pass
            pass
        else: # Quando aux == 0, di == 23 e df == 67, por isso o and
            if di <= dp[aux] and dp[aux] < df:
                temp+=1
                hsTemp.append((hs[aux]))
                tpTemp.append((tp[aux]))
                dpTemp.append((dp[aux]))
                dateTemp.append((date[aux]))
                pass

        pass
    d.append(dict(hs = np.asarray(hsTemp),tp = np.asarray(tpTemp),dp = np.asarray(dpTemp),date = np.asarray(dateTemp),di = di, df = df,dm = dm))
    pass

#sys.exit()
#%% Fazendo analise por cada direção
for auxT in range(len(d)):
    if len(d[auxT]['hs']) > CtrlMinDataLen:
        #%%
        print 'dir %d ------------------------'%auxT
        hs = d[auxT]['hs']
        tp = d[auxT]['tp']
        dirQuad = auxT
    
        #Determinando classes (bins) de Hs
        hsHist, hsBinsEdges = np.histogram(hs,bins=CtrlBinClasses,density = True)
        hsBinsMed = (hsBinsEdges[1] - hsBinsEdges[0])/2
        hsBinsDelta = (hsBinsEdges[1] - hsBinsEdges[0])
        hsBins = [hsBinsEdges[aux] + hsBinsMed for aux in range(len(hsBinsEdges)-1)]
        hsBinsRange = np.arange(hsBins[0],hsBins[-1],0.1)
        hsRange = np.arange(0.1,hsBins[-1]+hsBinsMed,0.1)
        
        
        #Determinando classes (bins) de Tp,
        tpHist, tpBinsEdges = np.histogram(tp,bins=CtrlBinClasses,density = True)
        tpBinsMed = (tpBinsEdges[1] - tpBinsEdges[0])/2
        tpDelta = (tpBinsEdges[1] - tpBinsEdges[0])
        tpBins = [tpBinsEdges[aux] + tpBinsMed for aux in range(len(tpBinsEdges)-1)]
        tpBinsRange = np.arange(tpBins[0],tpBins[-1],0.1)
        tpRange = np.arange(0.1,tpBins[-1]+tpBinsMed,0.1)      
        #==============================================================================
        #%%  Plotando histograma - frequencia
        #==============================================================================
        
        if Print1:
            fig.append(plt.figure())
            ax.append(fig[-1].add_subplot(111))
            pass
        
        
        #print (time.clock() - t0, 'hist freq')
        if Print1:
            tempHsHist,tempHsBins,temp = ax[-1].hist(hs,CtrlBinClasses,color='blue',normed = True,alpha=0.5,
                               stacked=True,histtype='step',label = 'Hist. Dens. Prob.')
            ax[-1].set_xticks(tempHsBins)
            temp = []
            for aux in range(len(tempHsBins)):
                if aux % 2 == 0:
                    temp.append('%.1f'%tempHsBins[aux])
                    pass
                else:
                    temp.append('')
                    pass
                pass
                
                pass
            ax[-1].set_xticklabels(temp)
            pass
        
        #==============================================================================
        # Plotando interp1d dos dados brutos
        #==============================================================================
        
        hsHistInterp = interp1d(hsBins,hsHist,kind='cubic')
        if Print1:    
            ax[-1].plot(hsBinsRange,hsHistInterp(hsBinsRange),'b',label=u'Emprica')
            pass
      
        
        #==============================================================================
        # Plotando a distribuição lognormal
        #==============================================================================
        
        hsXiLn = xiLn(hs)
        hsLambLn = lambLn(hs)
        
        #print (time.clock() - t0, 'dist lognormal met.')
        hsLnPDFy = [flnPDFHs(aux,hsXiLn,hsLambLn) for aux in hsRange]
        
        if Print1:
            ax[-1].plot(hsRange,hsLnPDFy,'k--',label='Dist. Lognormal')
            pass
        #==============================================================================
        ### Plotando distribuição Weibull scipy
        #==============================================================================
        #shape, loc, scale = weibull_min.fit(hs)
        #hsWb = weibull_min(shape,loc,scale)
        #ax[-1].plot(np.arange(0,5,0.1),hsWb.pdf(np.arange(0,5,0.1)),'red',ls='-.',label='Dist. Weibull SciPy')
        
        #shape, loc, scale = frechet_r.fit(hs)
        #hsFr = frechet_r(shape,loc,scale)
        #ax[-1].plot(np.arange(0,5,0.1),hsFr.pdf(np.arange(0,5,0.1)),'red',ls='-.',label='Dist. Freschet SciPy')
        #
        #shape, loc, scale = weibull_min.fit(hs)
        #hsWei = weibull_min(shape,loc,scale)
        #ax[-1].plot(np.arange(0,5,0.1),hsWei.pdf(np.arange(0,5,0.1)),'pink',ls='-.',label='Dist. Weibull SciPy')
        
        #def frechetPDFHs(aux,c,loc,scale):
        #    aux = (aux-loc)/scale
        #    return c * (aux**(c - 1)) * exp(-aux**c)
        #
        #c, loc, scale = frechet_r.fit(hs)
        #hsFrechPDFy = [frechetPDFHs(aux,c,loc,scale) for aux in np.arange(0.1,5,0.1)]
        #ax[-1].plot(np.arange(0.1,5,0.1),hsFrechPDFy,ls='-.',label='Dist. frechet')
        
        ### Funções Weibull
        ### lambw e alphaw
        
        hsWeiAlpha,hsWeiLamb = alphaw_lambw(hs)
        hsWeiPDFy = [weiPDFHs(aux,hsWeiAlpha,hsWeiLamb) for aux in np.arange(0.1,6,0.1)]
        if Print1:    
            ax[-1].plot(np.arange(0.1,6,0.1),hsWeiPDFy,'green',ls='-.',label='Dist. Weibull')
            ax[-1].grid(True)
            #legend = ax[-1].legend(loc='upper right', shadow=True, fontsize='medium')
            legend = ax[-1].legend()
            
            DirFigName = DirFigures + DirFileList[CtrlfileNum]
            DirFigName += 'DQ_%d-%d_HS_PDF.png'%(dirQuad+1,CtrlDirNum)
            
            fig[-1].savefig(DirFigName)
            pass
        #==============================================================================
        #%% Plotando histograma frequencia 3D
        #==============================================================================
        if PrintHist3D:
            #print (time.clock() - t0, 'histograma 3D')
            fig.append(plt.figure())
            ax.append(fig[-1].add_subplot(111, projection='3d'))
            
        
            #Fazendo os dados de histograma / xedges e yedges são os limites de cada bin
            dby, sxedges, syedges = np.histogram2d(hs, tp, bins=CtrlBinClasses, normed=True)
            
            #Construindo a grade para plotagem.
            xpos, ypos = np.meshgrid(sxedges[:-1] , syedges[:-1] ) #retorna a matriz de posições x e y
            
            xpos = xpos.flatten('F')
            ypos = ypos.flatten('F')
            zpos = np.zeros_like(xpos)
            
            # Construct arrays with the dimensions for the 16 bars.
            dx = (sxedges[1] - sxedges[0])* 0.8 * np.ones_like(zpos) #largura da barra em x
            dy = (syedges[1] - syedges[0])* 0.8 * np.ones_like(zpos)
            dz = dby.flatten() # altura da barra.
                 
            ax[-1].bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average',color='gray')
            ax[-1].set_xticks(xpos)
            ax[-1].set_yticks(ypos)
            ax[-1].set_xticklabels(['%2.1f'%aux for aux in xpos])
            ax[-1].set_yticklabels(['%2.1f'%aux for aux in ypos])
            ax[-1].set_xlabel('Hs')
            ax[-1].set_ylabel('Tp')
            ax[-1].set_zlabel('Dens. Probabilidade')
            
            DirFigName =  DirFigures + DirFileList[CtrlfileNum]
            DirFigName += 'DQ_%d-%d_HIST3D.png'%(dirQuad+1,CtrlDirNum)
            DirFigName.replace('/','\\')
            fig[-1].savefig(DirFigName)
            pass
        
        #==============================================================================
        #%% # Separando estados de Tp para cada faixa Hs
        #==============================================================================
        #print (time.clock() - t0, 'Separando estados de Tp para cada faixa Hs')
        tpHsData = [] #tpHsData guarda os dados de onda separando por classe de Hs
        for aux in range(CtrlBinClasses+1):
            temp = []
            for aux2 in range(len(tp)):
                if hsBinsEdges[aux] <= hs[aux2] and hs[aux2] < hsBinsEdges[aux] + hsBinsDelta:
                    temp.append(([hsBinsEdges[aux],hsBinsEdges[aux]+hsBinsDelta,hs[aux2],tp[aux2],dp[aux2]]))
                    pass
                pass
            tpHsData.append(np.asarray(temp))
            pass
        
        #Legenda tpHsData
        #tpHsData[0] #primeira classe de Hs
        #tpHsData[0][0] #primeira onda encontrada dentro da classe de Hs
        #tpHsData[0][0][0] #limite inferior da primeira classe de Hs
        #tpHsData[0][0][1] #limite superior da primeira classe de Hs
        #tpHsData[0][0][2] #Hs da primeira onda encontrada que pertence a primeira classe
        #tpHsData[0][0][3] #Tp da primeira onda encontrada que pertence a primeira classe
        #tpHsData[0][0][4] #Dp da primeira onda encontrada que pertence a primeira classe
        
        
        #%%
        tpHsDataLnPDF = [] #essa é a lista que irá armazenar os parametros sFln,sXi,sLamb
        tpHsDataLnPDFSciPy = [] # essa lista irá amarmazenar os parametros de fit pelo método do Scipy (tentando)
        
        
        if PrintTpHs:
            fig.append(plt.figure(figsize=(CtrlLinhas*2,CtrlColunas)))
            pass
        
        
        for aux in range(len(tpHsData)): #vai circular entre cada faixa de Tp definido dentro de cada Hs
            if len(tpHsData[aux]) <= CtrlLimMinColetasTpHs:
        
                pass
            else:
                if PrintTpHs:
                    ax.append(fig[-1].add_subplot(CtrlLinhas,CtrlColunas,aux+1))
                    ax[-1].set_title('Tp|Hs %1.1f<hs<%1.1f '%(hsBinsEdges[aux],hsBinsEdges[aux] + hsBinsDelta))
                    pass
                
                #calculo de xi e lambda do invervalo.
                tempTpHsXi = xiLn(tpHsData[aux][:,3]) # esse é o calculo do parametro xi do intervalo
                tempTpHsLamb = lambLn(tpHsData[aux][:,3]) #esse é o calculo do parametro lambda do intervalo
                temp1= lognorm.fit(tpHsData[aux][:,3],floc = 0)
                        
                #calculo da serie de dados pra plot.
                tempTpHsHist, tempTpHsBins = np.histogram(tpHsData[aux][:,3],bins=CtrlBinClasses,density=True) # histograma de frequencia
                tempTpHsBinsRange = np.arange(tempTpHsBins[0],tempTpHsBins[-1],0.1)
                tempTpHsy = [flnPDFHs(aux2,tempTpHsXi,tempTpHsLamb) for aux2 in tempTpHsBinsRange]
                tempTpHsBinsDelta = tempTpHsBins[1] - tempTpHsBins[0]
                
                tempTpHsBin = hsBinsEdges[aux] #+ hsBinsDelta/2
                tpHsDataLnPDF.append([0,tempTpHsXi,tempTpHsLamb,tempTpHsBin,temp1[0],np.log(temp1[2])])
        
                #plotando histograma acumulado
                tempTpHsHistCumsum = np.cumsum(tempTpHsHist)*tempTpHsBinsDelta
                
                tempIntepFunc = interp1d(tempTpHsBins[:-1]+tempTpHsBinsDelta,tempTpHsHistCumsum,kind='cubic')
                tempRange = np.arange(tempTpHsBins[0]+tempTpHsBinsDelta,tempTpHsBins[-2]+tempTpHsBinsDelta,0.1)
                
                if PrintTpHs:
                    ax[-1].plot(tempRange,tempIntepFunc(tempRange),'b',label = 'Empirica')
        #            ax[-1].bar(tempTpHsBins[:-1],tempTpHsHistCumsum,tempTpHsBinsDelta,linewidth = 0.01,alpha=0.1,label = 'Hist. dist. acumulada')
                    pass
                
                #plotando função de distribuição acumulada 
                tempRange = np.linspace(0.1,20,10)
                tempTpHsCDFy = [flnCDFHs(aux2,tempTpHsXi,tempTpHsLamb)[0] for aux2 in tempRange]
        #        tempRange = tempRange[:]+tempTpHsBinsDelta
                if PrintTpHs:
                    ax[-1].plot(tempRange,tempTpHsCDFy,'g',ls='--',label = 'Lognormal CDF')
                    pass
                pass
            pass
        
        if PrintTpHs:
            ax[-1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,shadow=True)
            fig[-1].tight_layout()
            
            DirFigName =  DirFigures + DirFileList[CtrlfileNum]
            DirFigName += 'DQ_%d-%d_TPHSCDF.png'%(dirQuad+1,CtrlDirNum)
            DirFigName.replace('/','\\')
            fig[-1].savefig(DirFigName)
            pass
        
        #==============================================================================
        #==============================================================================
        #==============================================================================
        #%% Plotando os valores de xi e lambda da distribução condicional
        #==============================================================================
        #==============================================================================
        #==============================================================================
        
        #print (time.clock() - t0, 'Plotando os valores de xi e lambda da distribuicao condicional')
        
        #Legenda tpHsDataLnPDF
        #tpHsDataLnPDF[0] #guarda dados do primeiro intervalo
        #tpHsDataLnPDF[0][0] #guarda array que forma a curva lognorm.pdf do primeiro intervalo
        #tpHsDataLnPDF[0][1] #guarda o xi do primeiro intervalo
        #tpHsDataLnPDF[0][2] #guarda o lambda do primeiro intervalo
        #tpHsDataLnPDF[0][3] #guarda o valor médio do sBin. (é a posição no eixo x que vou plotar xi e lambda)
        
        tpHsDataLnPDF = np.asarray(tpHsDataLnPDF)
        
        
        if PrintXiLamb:
            fig.append(plt.figure())
            ax.append(fig[-1].add_subplot(111))
            ax[-1].plot(tpHsDataLnPDF[:,3],tpHsDataLnPDF[:,1],'r*',label='xi')
            ax[-1].plot(tpHsDataLnPDF[:,3],tpHsDataLnPDF[:,2],'gp',label='lambda')
            pass
        #==============================================================================
        # Definindo as fuções com polyfit cubicas. O polyfit retorna lista com indices "a b c d" de um polinomio (ax³+bx²+cx+d)
        #==============================================================================
        #tpHsDataLnPDF[:,2] = np.mean(tpHsDataLnPDF[:,2])
        
        fTpHsXi = np.polyfit(tpHsDataLnPDF[:,3].tolist(),tpHsDataLnPDF[:,1].tolist(),CtrlPolyFitGrau) #indices da função que retorna valores de xi (de Hs)
        fTpHsLamb = np.polyfit(tpHsDataLnPDF[:,3].tolist(),tpHsDataLnPDF[:,2].tolist(),CtrlPolyFitGrau) #indices da função que retorna valores de lamb (de Hs)
        
                
        if PrintXiLamb:
            ax[-1].plot(tpHsDataLnPDF[:,3],[pol(fTpHsXi,aux) for aux in tpHsDataLnPDF[:,3]],'r',label=u'pol. 1º Xi Tp|Hs')
            ax[-1].plot(tpHsDataLnPDF[:,3],[pol(fTpHsLamb,aux) for aux in tpHsDataLnPDF[:,3]],'g',label=u'pol. 1º Lamb Tp|Hs')
            #dir(ax[-1])
            ax[-1].grid()
            ax[-1].set_xticks(tpHsDataLnPDF[:,3])
            #ax[-1].set_yticks()
            legend = ax[-1].legend(loc='upper right', shadow=True, fontsize='medium')
            
            DirFigName =  DirFigures + DirFileList[CtrlfileNum]
            DirFigName += 'DQ_%d-%d_LN_XI_lAMB.png'%(dirQuad+1,CtrlDirNum)
            DirFigName.replace('/','\\')
            fig[-1].savefig(DirFigName)
            pass
        
        np.polyval(fTpHsLamb,aux)
        #==============================================================================
        #%% Dados para pcolormesh e surface
        #==============================================================================
        
#        tempBins = (np.linspace(0.01,8.01,CtrlBinClasses)) #,np.linspace(0.01,25.01,CtrlBinClasses)
        tempBins = (CtrlBinClasses,CtrlBinClasses)
#        tempBins = (10,20)
        tempRange = ([0.1,CtrlRangeHs],[0.1,CtrlRangeTp]) 
        dby2, hsBinsEdges2, tpBinsEdges2 = np.histogram2d(hs, tp, bins = tempBins,range=tempRange, normed=True)
        dby2 = dby2.T
        hsBinsEdges2 = hsBinsEdges2[:-1] + (hsBinsEdges2[1] - hsBinsEdges2[0])/2
        tpBinsEdges2 = tpBinsEdges2[:-1] + (tpBinsEdges2[1] - tpBinsEdges2[0])/2
       
        
        xpos2, ypos2 = np.meshgrid(hsBinsEdges2,tpBinsEdges2) #retorna a matriz de posições x e y
        
        #==============================================================================
        #%% Dist. Conj. Hs,Tp flnPDFHs_flnTpHs
        #==============================================================================
        
        temp = []
        for Hs in hsBinsEdges2:
            temp2 = []
            for Tp in tpBinsEdges2:
                temp2.append(flnPDFHs_flnTpHs(Tp,Hs,hsXiLn,hsLambLn,fTpHsXi,fTpHsLamb))
                pass
        #    print(time.clock() - t0,'Distribuicao TpHs Hs =' , Hs)
            temp.append(temp2)
            pass
        
        hs_TpHsLnPDFy = np.asarray(temp).T # probability density function da distribuição calculada.
        
#        print('integral flnPDFHs_flnTpHs: ',nquad(flnPDFHs_flnTpHs,([0,10],[0,20]),args=(hsXiLn,hsLambLn,fTpHsXi,fTpHsLamb)))
        
        #==============================================================================
        #%% Dist. Conj. Hs,Tp usando weiPDFHs_flnTpHs
        #==============================================================================
        
        temp = []
        for Hs in hsBinsEdges2:
            temp2 = []
            for Tp in tpBinsEdges2:
                temp2.append(weiPDFHs_flnTpHs(Tp,Hs,hsWeiAlpha,hsWeiLamb,fTpHsXi,fTpHsLamb))
                pass
            temp.append(temp2)
            pass
        
        hs_TpHsWeiPDFy = np.asarray(temp).T # probability density function da distribuição calculada.

#        print('integra weiPDFHs_flnTpHs: ',nquad(weiPDFHs_flnTpHs,([0,10],[0,20]),args=(hsWeiAlpha,hsWeiLamb,fTpHsXi,fTpHsLamb)))
        
        # =============================================================================
        #%% Criando serie de Hs e Tp baseados na função densidade de probabilidade conjunta Lognormal       
        # =============================================================================
        from scipy.stats import pearsonr
        
        #%%
        hsList = [] 
        tpList = []
        
        xtemp2 = xpos2.flatten()
        ytemp2 = ypos2.flatten()
        
        for ix1,ax1 in enumerate(hs_TpHsLnPDFy.flatten()): #para cada probabilidade de acontecer um determinado Hs e Tp em conjunto
            for _ in range(int(round(1e1 * ax1))): 
                hsList.append(xtemp2[ix1])
                tpList.append(ytemp2[ix1])
        
        
             
        rho  = pearsonr(hsList,tpList)[0]
        rho3 = pearsonr(hs,tp)[0]
#        print 'rho1 (flnPDF) = ',rho
#        print 'rho3 (Brutos) = ', rho3        
             
  
        #%% Criando serie de Hs e Tp baseados na função densidade de probabilidade conjunta Weibull
                
        hsListW = []
        tpListW = []
        
        xtemp2W = xpos2.flatten()
        ytemp2W = ypos2.flatten()
        
        
        for ix1,ax1 in enumerate(hs_TpHsWeiPDFy.flatten()): #para cada probabilidade de acontecer um determinado Hs e Tp em conjunto
            for _ in range(int(round(1e1 * ax1))): 
                hsListW.append(xtemp2[ix1])
                tpListW.append(ytemp2[ix1])
                
                
        rhow  = pearsonr(hsListW,tpListW)[0]
        rhow3 =  rho3
#        print 'rho1 (WeiPDF) = ',rhow
#        print 'rho3 (Brutos) = ', rhow3
            
        #%% Dist. Conj. Hs,Tp Ln Scipy
        
        #def flnHs_flnTpHsSci(TpAux,HsAux,s,scale,fs,fscale):
        #    s2 = pol(fs,HsAux)
        #    scale2 = np.exp(pol(fscale,HsAux))
        #    return lognorm.pdf(HsAux,s,0,scale) * lognorm.pdf(TpAux,s2,0,scale2)
        #    
        #s,loc,scale = lognorm.fit(hs,floc=0)
        #temp = []
        #for Tp in tpBinsEdges2:
        #    temp2 = []
        #    for Hs in hsBinsEdges2:
        #        temp2.append(flnHs_flnTpHsSci(Tp,Hs,s,scale,fTpHsXi,fTpHsLamb))
        #        pass
        ##    print(time.clock() - t0,'Distribuicao TpHs Hs =' , Hs)
        #    temp.append(temp2)
        #    pass
        #
        #hs_TpHsLnPDFy = np.asarray(temp) # probability density function da distribuição calculada.
        #hs_TpHsLnPDFy = np.transpose(hs_TpHsLnPDFy)
        

        
        
        #%% Dist. Conj. Hs,Tp Wei_min Scipy
        
        #def weiHs_flnTpHsSci(TpAux,HsAux,c,scale,fs,fscale):
        #    s2 = pol(fs,HsAux)
        #    scale2 = np.exp(pol(fscale,HsAux))
        #    return weibull_min.pdf(HsAux,c,0,scale) * lognorm.pdf(TpAux,s2,0,scale2)
        ##    return lognorm.pdf(TpAux,s2,scale2)
        #    
        #c,loc,scale = weibull_min.fit(hs,floc = 0)
        #temp = []
        #for Tp in tpBinsEdges2:
        #    temp2 = []
        #    for Hs in hsBinsEdges2:
        #        temp2.append(weiHs_flnTpHsSci(Tp,Hs,c,scale,fTpHsXi,fTpHsLamb))
        #        pass
        ##    print(time.clock() - t0,'Distribuicao TpHs Hs =' , Hs)
        #    temp.append(temp2)
        #    pass
        #
        #hs_TpHsWeiPDFy = np.asarray(temp) # probability density function da distribuição calculada.
        #hs_TpHsWeiPDFy = hs_TpHsWeiPDFy.transpose()
        
        #%% Dist. Conj. Hs,Tp Wei_max Scipy
        
        #def weiHs_flnTpHsSci(TpAux,HsAux,c,loc,scale,fs,fscale):
        #    s2 = pol(fs,HsAux)
        #    scale2 = np.exp(pol(fscale,HsAux))
        #    return weibull_max.pdf(HsAux,c,loc,scale) * lognorm.pdf(TpAux,s2,0,scale2)
        ##    return lognorm.pdf(TpAux,s2,scale2)
        #    
        #c,loc,scale = weibull_max.fit(hs)
        #temp = []
        #for Tp in tpBinsEdges2:
        #    temp2 = []
        #    for Hs in hsBinsEdges2:
        #        temp2.append(weiHs_flnTpHsSci(Tp,Hs,c,loc,scale,fTpHsXi,fTpHsLamb))
        #        pass
        ##    print(time.clock() - t0,'Distribuicao TpHs Hs =' , Hs)
        #    temp.append(temp2)
        #    pass
        #
        #hs_TpHsWeiPDFy = np.asarray(temp) # probability density function da distribuição calculada.
        #hs_TpHsWeiPDFy = hs_TpHsWeiPDFy.transpose()
        
        #%% Dist. Conj. Hs,Tp Gumbel Scipy
        
        #def gumHs_flnTpHsSci(TpAux,HsAux,c,scale,fs,fscale):
        #    s2 = pol(fs,HsAux)
        #    scale2 = np.exp(pol(fscale,HsAux))
        #    return gumbel_r.pdf(HsAux,c,scale) * lognorm.pdf(TpAux,s2,0,scale2)
        ##    return lognorm.pdf(TpAux,s2,scale2)
        #    
        #c,scale = gumbel_r.fit(hs)
        #temp = []
        #for Tp in tpBinsEdges2:
        #    temp2 = []
        #    for Hs in hsBinsEdges2:
        #        temp2.append(gumHs_flnTpHsSci(Tp,Hs,c,scale,fTpHsXi,fTpHsLamb))
        #        pass
        ##    print(time.clock() - t0,'Distribuicao TpHs Hs =' , Hs)
        #    temp.append(temp2)
        #    pass
        #
        #hs_TpHsWeiPDFy = np.asarray(temp) # probability density function da distribuição calculada.
        #hs_TpHsWeiPDFy = hs_TpHsWeiPDFy.transpose()
        
        
        
        #%% plotando pcolormesh da flnPDFHs_flnTpHs
        
        
        #%%Usando as funções de covariancia e correlação do Numpy
        #print (time.clock() - t0, 'calculando rmse')
        
        hs_TpHsLnPDFyf = hs_TpHsLnPDFy.flatten()
        shistf = dby2.flatten()
        
        hs_TpHsLnPDFyf.mean()
        hs_TpHsLnPDFyf.std()
        hs_TpHsLnPDFyf.max()
        shistf.mean()
        shistf.std()
        shistf.max()
        
        temp = 0
        for aux in xrange(len(hs_TpHsLnPDFyf)):
        #    temp += (hs_TpHsLnPDFyf[aux] - shistf[aux])**2
            temp += ( shistf[aux] - hs_TpHsLnPDFyf[aux])**2
            pass
        rmseLn = np.sqrt(temp/len(hs_TpHsLnPDFyf))
#        print 'rmse Lognorm = ', rmseLn
        
        #%% Usando a distribuicao de Weibull
        
        #print (time.clock() - t0, 'calculando rmse')
        
        hs_TpHsWeiPDFyf = hs_TpHsWeiPDFy.flatten()
        shistf = dby2.flatten()
        
        hs_TpHsWeiPDFyf.mean()
        hs_TpHsWeiPDFyf.std()
        hs_TpHsWeiPDFyf.max()
        shistf.mean()
        shistf.std()
        shistf.max()
        
        temp = 0
        for aux in xrange(len(hs_TpHsWeiPDFyf)):
        #    temp += (hs_TpHsWeiPDFyf[aux] - shistf[aux])**2
            temp += ( shistf[aux] - hs_TpHsWeiPDFyf[aux])**2
            pass
        rmseWei = np.sqrt(temp/len(hs_TpHsWeiPDFyf))
#        print 'rmse Weibull = ', rmseWei
                    
            #%% print Contour
        if PrintContourf:
            
            
            cmap2 = cm.get_cmap('jet')
            cmap2.set_under('w')
            cmapbias = cm.get_cmap('seismic')
            levels = np.arange(0.01,0.3,0.03)
            levelsBias = np.arange(-0.11,0.11,0.02)

            fig.append(plt.figure(figsize=(12,9)))
            
            
            #Fazendo o x ticks labels
            xtl = []
            for aux in xrange(len(xpos2[0])):
                if aux % 2 == 0:
                    xtl.append('%.1f'%xpos2[0][aux])
                    pass
                else:
                    xtl.append('')
                    pass
                pass
            
                      
            ax.append(fig[-1].add_subplot(331))
            splot2=ax[-1].contour(xpos2,ypos2, dby2,cmap=cmap2,vmin = 0.001, vmax = 0.3,levels=levels)
            ax[-1].set_xlabel('Hs')
            ax[-1].set_ylabel('Tp')
            ax[-1].set_title('Histograma Densidade Probabilidade')
            ax[-1].set_xticks(xpos2[0])
            ax[-1].set_xticklabels(xtl)
            ax[-1].grid()
            fig[-1].colorbar(splot2, shrink=0.5, aspect=5,)
            
            

############ Modelo Distribuicao Conjunta  ####################################
            
            H = hs
            T = tp
            X= hs
            
            
            def best_dist_fit(X,**kwargs):
                
                bins=30
                rangex = (0,30)
                
                if kwargs.has_key('bins'):
                    bins = kwargs['bins']
                if kwargs.has_key('rangex'):
                    rangex = (0,kwargs['rangex'])
                    pass
                
                #fazendo o histograma de X
                X_hist_y,X_hist_x = np.histogram(X,bins,density = True,range = rangex)
                X_hist_xM = X_hist_x[:-1] + X_hist_x[0:2].mean()
                X_hist_yM = np.cumsum(X_hist_y) * X_hist_x[0:2].mean()*2
                
                #aplicando dist. lognorm.
                lnfit = lognorm.fit(X,floc=0)
                X_ln = lognorm(lnfit[0],lnfit[1],lnfit[2])

                #aplicando dist. weibull_min
                weifit = weibull_min.fit(X,floc=0)
                X_wei = weibull_min(weifit[0],weifit[1],weifit[2])    

                #aplicando gamma
                gamfit = scipy.stats.gamma.fit(X,floc=0)
                X_gam = scipy.stats.gamma(gamfit[0],gamfit[1],gamfit[2])
                
                #aplicando rayleigh
                rayfit = scipy.stats.rayleigh.fit(X,floc=0)     
                X_ray = scipy.stats.rayleigh(rayfit[0],rayfit[1])
                
                #aplicando gumbel
                gumfit = scipy.stats.gumbel_l.fit(X,floc=0)
                X_gum = scipy.stats.gumbel_l(gumfit[0],gumfit[1])
                
                def f_rmse(X,Y):
                    return np.sqrt(((X-Y)**2).mean())
                
                #comparando a serie bruta com a lognormal                
                rmseln = f_rmse(X_hist_yM,X_ln.cdf(X_hist_xM))
                rmsewei = f_rmse(X_hist_yM,X_wei.cdf(X_hist_xM))
                rmsegam = f_rmse(X_hist_yM,X_gam.cdf(X_hist_xM))
                rmseray = f_rmse(X_hist_yM,X_ray.cdf(X_hist_xM))
                rmsegum = f_rmse(X_hist_yM,X_gum.cdf(X_hist_xM))

                
#                fig = plt.figure()
#                ax = fig.add_subplot(111)
#                ax.plot(X_hist_xM,X_hist_yM,label='brut')
#                ax.plot(X_hist_xM,X_ln.cdf(X_hist_xM),label='ln')
#                ax.plot(X_hist_xM,X_wei.cdf(X_hist_xM),label='wei')
#                ax.plot(X_hist_xM,X_gam.cdf(X_hist_xM),label='gam')
#                ax.plot(X_hist_xM,X_ray.cdf(X_hist_xM),label='ray')
#                ax.plot(X_hist_xM,X_gum.cdf(X_hist_xM),label='gum')
#                ax.legend()
#                             
                
                t1 = pd.Series([rmseln,rmsewei,rmseray,rmsegam,rmsegum],
                               index=['lognormal','weibull_min','rayleigh','gamma','gumbel_l'])
                
                return t1.loc[t1.min() == t1].index[0],t1
            
            best_dist_fit(hs,bins=14,rangex=10)
            best_dist_fit(tp,bins=14,rangex=20)      

            def mdc(H,T):
                
#                fH = lambda x
#                fT_H = 0 
                              
                pass

#            def

               

                
                
                #parametro do histograma de T

#                tempTpHsHistCumsum = np.cumsum(tempTpHsHist)*tempTpHsBinsDelta
#                tempIntepFunc = interp1d(tempTpHsBins[:-1]+tempTpHsBinsDelta,tempTpHsHistCumsum,kind='cubic')
                
                pass
            
            
            
            
            
            
            ax.append(fig[-1].add_subplot(332))
            splot1=ax[-1].contour(xpos2,ypos2, hs_TpHsLnPDFy ,cmap=cmap2, vmin = 0.001, vmax = 0.3,levels=levels)
            ax[-1].set_xlabel('Hs')
            ax[-1].set_ylabel('Tp')
            ax[-1].set_title('Funcao Distribuicao Conjunta (Lognormal)')
            ax[-1].set_xticks(xpos2[0])
            ax[-1].set_xticklabels(xtl)
            ax[-1].grid()
            fig[-1].colorbar(splot1, shrink=0.5, aspect=5)
            
            
            ax.append(fig[-1].add_subplot(333))
            bias = dby2 - hs_TpHsLnPDFy
            splot2=ax[-1].contour(xpos2,ypos2, bias,cmap=cmapbias, vmin = -0.15, vmax = 0.15,levels=levelsBias)
            ax[-1].set_xlabel('Hs')
            ax[-1].set_ylabel('Tp')
            ax[-1].set_title('Bias Hist. Dens. Prob. - PDF (Lognormal)')
            ax[-1].set_xticks(xpos2[0])
            ax[-1].set_xticklabels(xtl)
            ax[-1].grid()
            fig[-1].colorbar(splot2, shrink=0.5, aspect=5)
            
            
############ Nataf ############################################################
            nat = nataf()
            
            nat.fit(hs,tp,rangehs = CtrlRangeHs,rangetp = CtrlRangeTp,bins= CtrlBinClasses)
            temp = nat.pdf_series(hs,tp)
            
            dfNataf = pd.DataFrame(temp['df'])
            xpos = temp['xpos']
            ypos = temp['ypos']
            
            
            ax.append(fig[-1].add_subplot(338))
            ax[-1].contour(xpos,ypos,dfNataf)
            ax[-1].set_xlabel('Hs')
            ax[-1].set_ylabel('Tp')
            ax[-1].set_title('Metodo de Nataf')
            ax[-1].set_xticks(xpos2[0])
            ax[-1].set_xticklabels(xtl)
            ax[-1].grid()
            fig[-1].colorbar(splot1, shrink=0.5, aspect=5)

            
            rmseNataf = ((dby2.flatten() - dfNataf.values.flatten())**2).mean() **.5
            
                       
            ax.append(fig[-1].add_subplot(339))
            bias = dby2 - dfNataf
            splot2=ax[-1].contour(xpos2,ypos2, bias,cmap=cmapbias, vmin = -0.15, vmax = 0.15, levels=levelsBias)
            ax[-1].set_xlabel('Hs')
            ax[-1].set_ylabel('Tp')
            ax[-1].set_title('Bias Hist. Den. Prob. - Nataf')
            ax[-1].set_xticks(xpos2[0])
            ax[-1].set_xticklabels(xtl)
            ax[-1].grid()
            fig[-1].colorbar(splot2, shrink=0.5, aspect=5)
                       
            
            ax.append(fig[-1].add_subplot(334))
            ax[-1].plot(dby2.flatten(),label = 'Serie Bruta',ls='solid')
            ax[-1].plot(hs_TpHsLnPDFy.flatten(),label = 'Serie LogPDF.',ls='dashed')
            ax[-1].plot(hs_TpHsWeiPDFy.flatten(),label = 'Serie WeiPDF',ls='dashdot')
            ax[-1].plot(dfNataf.values.flatten(),label = 'Serie Nataf',ls='dotted')
            ax[-1].set_title('PDF todas as series "Flat"')
            ax[-1].set_xlabel('Pontos f(Hs,Tp)')
            ax[-1].set_ylabel('Densidade Probabilidade')
            ax[-1].legend()
            
           
            dby3= dby2.flatten()
            
            ax.append(fig[-1].add_subplot(335))
            splot1=ax[-1].contour(xpos2,ypos2, hs_TpHsWeiPDFy ,cmap=cmap2, vmin = 0.001, vmax = 0.3,levels=levels)
            ax[-1].set_xlabel('Hs')
            ax[-1].set_ylabel('Tp')
            ax[-1].set_title('Funcao Distribuicao Conjunta (Weibull)')
            ax[-1].set_xticks(xpos2[0])
            ax[-1].set_xticklabels(xtl)
            ax[-1].grid()
            fig[-1].colorbar(splot1, shrink=0.5, aspect=5)
            
            
            ax.append(fig[-1].add_subplot(336))
            bias = dby2 - hs_TpHsWeiPDFy
            splot2=ax[-1].contour(xpos2,ypos2, bias,cmap=cmapbias, vmin = -0.15, vmax = 0.15, levels=levelsBias)
            ax[-1].set_xlabel('Hs')
            ax[-1].set_ylabel('Tp')
            ax[-1].set_title('Bias Hist. Den. Prob. - PDF(Weibull)')
            ax[-1].set_xticks(xpos2[0])
            ax[-1].set_xticklabels(xtl)
            ax[-1].grid()
            fig[-1].colorbar(splot2, shrink=0.5, aspect=5)
            
            
                        
            ax.append(fig[-1].add_subplot(337))
            temp = 0.94
            ax[-1].text(0.05,temp,u'Arquivo: %s'%DirFileList[CtrlfileNum])
            temp -= 0.06
            ax[-1].text(0.05,temp,u'Nº coletas %d || %d° < %s < %d°'%(len(hs),d[auxT]['di'],direcaofinal[d[auxT]['dm']],d[auxT]['df']))
            temp -= 0.1
            ax[-1].text(0.05,temp,u'rho (Brutos) = %.3f'%(rho3))
#            temp -= 0.05
#            ax[-1].text(0.05,temp,u'rho (flnPDF) = %.3f'%(rho))
#            temp -= 0.05
#            ax[-1].text(0.05,temp,u'rho (WeiPDF) = %.3f'%(rhow))
#            temp -= 0.05
#            ax[-1].text(0.05,temp,u'rho (Nataf) = %.3f'%(rhoNataf))
            
                      
            temp -= 0.1
            
            ax[-1].text(0.05,temp,u'rmse Lognorm = %.3f'%(rmseLn))
            temp -= 0.05
            ax[-1].text(0.05,temp,u'rmse Weibull = %.3f'%(rmseWei))        
            temp -= 0.05
            ax[-1].text(0.05,temp,u'rmse Nataf = %.3f'%(rmseNataf))        
            
            
            ran = [[1,8],[1,16]]
            if PrintRho:
                
                print('Achando correlacao rholn')
                #hs
                mh = nquad(lambda h,t: h*flnPDFHs_flnTpHs(t,h,hsXiLn,hsLambLn,fTpHsXi,fTpHsLamb),ran)[0]
                m2h = nquad(lambda h,t: h**2 * flnPDFHs_flnTpHs(t,h,hsXiLn,hsLambLn,fTpHsXi,fTpHsLamb),ran)[0]
                dph = np.sqrt(m2h - mh**2) #desvio padrao
                #tp
                mt = nquad(lambda h,t: t*flnPDFHs_flnTpHs(t,h,hsXiLn,hsLambLn,fTpHsXi,fTpHsLamb),ran)[0]
                m2t = nquad(lambda h,t: t**2 *flnPDFHs_flnTpHs(t,h,hsXiLn,hsLambLn,fTpHsXi,fTpHsLamb),ran)[0]
                dpt = np.sqrt(m2t - mt**2)    
                #
                mht = nquad(lambda h,t: h * t *flnPDFHs_flnTpHs(t,h,hsXiLn,hsLambLn,fTpHsXi,fTpHsLamb),ran)[0]       
                rholn = (mht - mh*mt) / (dph * dpt)
    
    
                print('Achando correlacao rhowei')
                mh = nquad(lambda h,t: h*weiPDFHs_flnTpHs(t,h,hsWeiAlpha,hsWeiLamb,fTpHsXi,fTpHsLamb),ran)[0]
                m2h = nquad(lambda h,t: h**2 *weiPDFHs_flnTpHs(t,h,hsWeiAlpha,hsWeiLamb,fTpHsXi,fTpHsLamb),ran)[0]
                dph = np.sqrt(m2h - mh**2) #desvio padrao
                #tp
                mt = nquad(lambda h,t: t*weiPDFHs_flnTpHs(t,h,hsWeiAlpha,hsWeiLamb,fTpHsXi,fTpHsLamb),ran)[0]
                m2t = nquad(lambda h,t: t**2 *weiPDFHs_flnTpHs(t,h,hsWeiAlpha,hsWeiLamb,fTpHsXi,fTpHsLamb),ran)[0]
                dpt = np.sqrt(m2t - mt**2)    
                #
                mht = m2t = nquad(lambda h,t: h * t * weiPDFHs_flnTpHs(t,h,hsWeiAlpha,hsWeiLamb,fTpHsXi,fTpHsLamb),ran)[0]
                rhowei = (mht - mh*mt) / (dph * dpt)
                
                
                print('Achando correlacao rhonat')
                #Hs
                mh = nquad(lambda h,t: h* nat.pdf(h,t),ran)[0] #media
                m2h = nquad(lambda h,t: h**2 * nat.pdf(h,t),ran)[0]
                dph = np.sqrt(m2h - mh**2) #desvio padrao
                #Tp
                mt = nquad(lambda h,t: t * nat.pdf(h,t),ran)[0]
                m2t = nquad(lambda h,t: t ** 2 * nat.pdf(h,t),ran)[0]
                dpt = np.sqrt(m2t - mt**2)    
                #
                mht = nquad(lambda h,t: h * t * nat.pdf(h,t),ran)[0]    
                rhonat = (mht - mh*mt) / (dph * dpt)
                        
            
###              rholn = (nquad(lambda x,y: x*y*flnPDFHs_flnTpHs(y,x,hsXiLn,hsLambLn,fTpHsXi,fTpHsLamb),[[1,8],[1.0,16]])[0] - hs.mean()*tp.mean())/(hs.std()*tp.std())
###              rhowei = (nquad(lambda x,y,t0,t1,t2,t3: x*y*weiPDFHs_flnTpHs(y,x,t0,t1,t2,t3),[[1,8],[1.0,16]],args=(hsWeiAlpha,hsWeiLamb,fTpHsXi,fTpHsLamb))[0] - hs.mean()*tp.mean())/(hs.std()*tp.std())
###                rhonat = (nquad(lambda x,y: x*y*nat.pdf(x,y),[[1,8],[1.0,16]])[0] - hs.mean()*tp.mean())/(hs.std()*tp.std())
                
                temp -= 0.1
                ax[-1].text(0.05,temp,u'Corr. log: %.3f'%rholn)
                temp -= 0.05
                ax[-1].text(0.05,temp,u'Corr. Wei: %.3f'%rhowei)
                temp -= 0.05
                ax[-1].text(0.05,temp,u'Corr Nataf: %.3f'%rhonat)
                
                temp -= 0.05
                ax[-1].text(0.05,temp,u'Htype Nataf: %s'%nat.Htype)
                temp -= 0.05
                ax[-1].text(0.05,temp,u'Ttype Nataf: %s'%nat.Ttype)
            
            if PrintI:
                print('Achando Integracao iLn')
                iln = nquad(lambda x,y: flnPDFHs_flnTpHs(y,x,hsXiLn,hsLambLn,fTpHsXi,fTpHsLamb),[[1,8],[0.5,16]])[0]
                print('Achando Integracao iweo')
                iwei = nquad(lambda x,y,t0,t1,t2,t3: weiPDFHs_flnTpHs(y,x,t0,t1,t2,t3),[[1,8],[1.0,16]],args=(hsWeiAlpha,hsWeiLamb,fTpHsXi,fTpHsLamb))[0]
                print('Achando Integracao inat')
                inat = nquad(lambda x,y: nat.pdf(x,y),[[1,8],[1.0,16]])[0]
    
                
                temp -= 0.1
                ax[-1].text(0.05,temp,u'Integral func Ln: %.3f'%iln)
                temp -= 0.05
                ax[-1].text(0.05,temp,u'Integral func Wei: %.3f'%iwei)
                temp -= 0.05
                ax[-1].text(0.05,temp,u'Integral func Nataf: %.3f'%inat)
                
                
                
            fig[-1].tight_layout()
            
#            plt.plot(dfNataf.values.flatten('F'))
#            plt.plot(dby2.flatten('F'))
#            plt.plot(hs_TpHsLnPDFy.flatten('F'))
            
            
            tempFigName = 'DQ_%d-%d_'%(dirQuad+1,CtrlDirNum)
#            DirFigName =  DirFigures + (DirFileList[CtrlfileNum]).replace('_',tempFigName)
            DirFigName =   DirFigures + (DirFileList[CtrlfileNum]) + tempFigName
            DirFigName += '_CONTOUR.png'
            print DirFigName
            fig[-1].savefig(DirFigName)

#
#            print('Corr. PDF log: ',
#            print('Corr. PDF Wei: ',
#            print('Corr. Nataf: ', 
#        
            
        # =============================================================================
        #%% Posicao do dado - Mapa  
        # =============================================================================
            
            if auxTT == 0:
                auxTT += 1
                fig.append(plt.figure(figsize=(10,4)))  
    #            ax.append(fig[-1].add_subplot(111))
                
                
    #            gs = gridspec.GridSpec(3,3)            
    #            ax.append(fig[-1].add_subplot(gs[2,1:]))
                gs = gridspec.GridSpec(1,1) 
                ax.append(fig[-1].add_subplot(gs[:]))            
                
                m = Basemap(projection='merc',llcrnrlat=-60,urcrnrlat=60,\
                llcrnrlon=-180.,urcrnrlon=180,lat_ts=-23.,resolution='c')
                m.drawcoastlines()
                m.drawstates()
                m.drawcountries()
                parallels = np.arange(-60.,60,15.)
                m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
                meridians = np.arange(-120.,40.,15.)
                m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=6)
                temp1,temp2 = m(lon,lat)
                m.drawmapboundary(fill_color='#99ffff')
                m.fillcontinents(color='#cc9966',lake_color='#99ffff')
                m.scatter(temp1,temp2,120,marker='x',color='#ff0000')
                ax[-1].set_title(u'Posição do dado')
                fig[-1].tight_layout()
                DirFigName =   DirFigures + (DirFileList[CtrlfileNum])
                DirFigName += '_MAPA.png'
                print DirFigName
                fig[-1].savefig(DirFigName)
            
           
        
        #==============================================================================
        #%% Plotando distribuição condicional TpHs. Apenas demosntrativo
        #==============================================================================
        if PrintTpHs3D:
            # Esse é um passo apenas demonstrativo. Não tem nenhum calculo essencial para o programa.
            #print (time.clock() - t0, 'histograma 3D')
            fig.append(plt.figure())
            ax.append(fig[-1].add_subplot(111, projection='3d'))
            
            #Fazendo os dados de histograma / xedges e yedges são os limites de cada bin
            temp1, temp2, temp3 = np.histogram2d(hs, tp, bins=CtrlBinClasses, normed=True)
            
            #Construindo a grade para plotagem.
            xpos3 = temp2[:-1]
            xpos3 = np.linspace(temp2[0],temp2[-1],20)
            ypos3 = temp3[:-1]
            ypos3 = np.linspace(temp3[0],temp3[-1],100)
            
            for aux in xrange(len(xpos3)):
                x = np.asarray([xpos3[aux] for temp4 in range(len(ypos3))])
                y = np.linspace(0.01,ypos3[-1],len(x))
                z = [flnPDFTpHs(aux2,xpos3[aux],fTpHsXi,fTpHsLamb) for aux2 in y]
                ax[-1].plot(x,y,z)
                pass
            
            ax[-1].set_xlabel('Hs')
            ax[-1].set_ylabel('Tp')
            ax[-1].set_zlabel('PDF')
              
#            plt.close(fig[-1])
            pass
        
        print '-------------------------------'
        pass
    pass
    
#%% Fim

print (time.clock() - t0, 'fim')
gc.collect()


sys.exit(10)

#%%
print('mh') #media
mh = nquad(lambda h,t: h* nat.pdf(h,t),[[1,10],[1,20]])[0]

print('m2h') 
m2h = nquad(lambda h,t: h**2 * nat.pdf(h,t),[[1,10],[1,20]])[0]

print('dph') #desvio padrao
dph = np.sqrt(m2h - mh**2)


####################################################
print('mt') #media 
mt = nquad(lambda h,t: t * nat.pdf(h,t),[[1,10],[1,20]])[0]

print('m2t') #
m2t = nquad(lambda h,t: t ** 2 * nat.pdf(h,t),[[1,10],[1,20]])[0]

print('dpt') #desvio padrao
dpt = np.sqrt(m2t - mt**2)

##################################################
print('mht') #
mht = nquad(lambda h,t: h * t * nat.pdf(h,t),[[1,10],[1,20]])[0]    

print('ro')
ro = (mht - mh*mt) / (dph * dpt)





# ro = E - mediaHs * mediaTp / 

#rhonat = nquad(lambda h,t: (h*t*nat.pdf(h,t) - h*nat.pdf(h,t) * t*nat.pdf(h,t)) /( (h ** 2 * nat.pdf(h,t) - (h * nat.pdf(h,t))**2) * (t ** 2 * nat.pdf(h,t) - (t * nat.pdf(h,t))**2) ),[[1,8],[1.0,16]])


#rhonat = nquad(lambda h,t: (h * t * nat.pdf(h,t) - h* nat.pdf(h,t) * t * nat.pdf(h,t)) / ( ( (h**2 * nat.pdf(h,t)) - (h* nat.pdf(h,t))**2 ) * ( (t ** 2 * nat.pdf(h,t)) -  (t * nat.pdf(h,t))**2) ) ,[[1,8],[1.0,16]])

#%%

hs.std()
hs.mean()


nquad(lambda h,t: t * nat.pdf(h,t),[[1,10],[1,20]])


tp.std()
tp.var()
dir(tp)


