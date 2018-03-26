# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 13:45:03 2017

@author: VPNUser
"""

from scipy.stats import norm,weibull_min
from scipy.special import gamma
from scipy.optimize import fsolve
from scipy.integrate import quad
import numpy as np
from scipy.optimize import root




#%% Class ln()
        
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

#%% Class wei

class wei:
    

    
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