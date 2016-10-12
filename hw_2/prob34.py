#!/usr/bin/python
# First bulit by Samantha Chang-Chun Chen 10/08/2016
# Class PHYS 239

import numpy as np
import astropy
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
#########################################
#########################################
    
freq0=320 
variance=1 
resol=1000 
nstep=100

x = np.linspace(freq0 - 5*variance,freq0 + 5*variance,resol)
y=mlab.normpdf(x,freq0,variance)

crsec1 = 3.33*10**-24
crsec2 = 3.33*10**-21
crsec3 = 3.33*10**-18

plt.semilogy(x,y*crsec1,'b-',label=r'$\sigma_0$='+ str(crsec1))
plt.semilogy(x,y*crsec2,'r-',label=r'$\sigma_0$='+ str(crsec2))
plt.semilogy(x,y*crsec3,'g-',label=r'$\sigma_0$='+ str(crsec3))

plt.xlabel(r'Frequency ($GHz$)', size = 10)
plt.ylabel(r'Cross section ($\log_{10}(cm^2)$)', size = 10)
plt.legend()
plt.title(r'Cross section distribution ($\nu_0$ = ' + str(freq0) + 
              r' $GHz$, variance ' + str(variance) + r' $GHz$)', size=14)

plt.show()




#####################################
################Prob.4###########


###############################################
I0 = [1,0,1,3,1,3]
S = [2,2,2,2,2,2]
n = 1
variance = [50,1,1,1,1,1]


ds = 1*3*10**18 ##[cm] 100 steps.

crsec00 = [3*10**-18,3*10**-21,3*10**-21,3*10**-21,3*10**-20,3*10**-20]  
## for tau<1, corss section<10**-21 ###
panel = list ('abcdef')

plt.figure()
for i in range(6):
    plt.subplot(2,3,i+1)
    x = np.linspace(freq0 - 5*variance[i], freq0 + 5*variance[i], resol)
    y=mlab.normpdf(x,freq0,variance[i])
    crsec = crsec00[i]*y ## Gaussian distb
    alpha = n * crsec00[i]*y
    for j in range(nstep):
        if j == 0:
            I = I0[i]
        I += (S[i] - I) * alpha * ds
      
    #######################
    plt.subplot(2,3,i+1)  
 
    plt.plot(x, I, 'r', label=r'$I_\nu(D)$')      
    plt.plot(x, [S[i]]*resol , 'g--', label=r'$S_\nu$')
    plt.plot(x, [I0[i]]*resol , 'b--', label=r'$I_\nu(0)$')
    plt.ylim(min(S[i], I0[i]) - 0.2, max(S[i], I0[i]) + 0.2)
    plt.xlim(314,326)

    plt.title (panel[i])

plt.show()

## (a) tau >>1##

## (b) I0=0, tau <1 ##
## (c) I0 < S, tau<1 ##

