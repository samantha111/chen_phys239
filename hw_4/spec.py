#!/usr/bin/python
#first built by Samantha Chang-Chun Chen 12/08/2016
# Class PHYS 239
### HW4
#from __future__ import absolute_import, division, print_function, unicode_literals
from numpy import *
#from scipy.optimize import curve_fit
#from lmfit import minimize, Parameters, Parameter, report_fit
from pylab import *
from scipy import *
from scipy.special import gamma

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
####== CGS ==========######

from astropy.constants import  e, m_e, h, k_B,L_sun, c
L_sun = L_sun.cgs.value
c = c.cgs.value
q = e.gauss.value
m_e = m_e.cgs.value
h = h.cgs.value
k_B = k_B.cgs.value



##data from M82###

M82 = np.loadtxt("m82spec.dat", skiprows=1)

plt.figure()
plt.loglog(M82[:, 0], M82[:, 1])
plt.show()


## star light ####
#############
fig5 = np.loadtxt("fig5", skiprows=200)


f5x = fig5[:, 0].reshape(-1, 1)* 1E-4 ## convert to \nu m

f5y = fig5[:, -3] - np.log10(3.83* 1E29) ## 700Myr
f5y = (fig5[:, 0]**2) * (10**(f5y))  / (3E18) ### devided by speed of light!
f5y = f5y.reshape(-1, 1)


##  dust ##
###
#a =   1E-3  ###SiC21: 0.001 - 10 micron
#
#d = np.loadtxt("dust", skiprows=10)
#Q = d[:, 1].reshape(-1, 1)
#dx = d[:,0].reshape(-1, 1)

#vol = 1E48  ## arbitrary constant used to fit the data


#B = 2 * h * c* 1E4 / wl**3 / (np.exp(h*1E4/dust[:,0]/T/k_B) - 1) ### Plank function


#:dy = 3* np.pi * B * Q / a  *  vol



### free-free ###

G = 3E-6## constant
temp = 45 ### hot, 10^4 is warm
Pff = G * np.exp ( -h * c *1E4 /k_B / temp /M82[:, 0])



## synchrotron ##
p = 2.28      ## gause by eyes
C = 10**12
Bsin = 3

lamb = np.linspace (10**(-1), 10**5,100)
Pt = C * Bsin * gamma(p/4 + 19/12) * gamma(p/4 - 1/12) * (m_e ** (-p/2-1/2)) * (q**(5/2+p/2)) * (c**(-p/2-3/2))* \
              (lamb * Bsin) ** ((p-1)/2)



### plot figure #############

plt.figure()

plt.loglog(lamb, Pt, label =  "synchrotron radiation")
plt.loglog(M82[:, 0], M82[:, 1], label =  "data from observation")
plt.loglog(f5x,f5y, label =  "star light")
plt.loglog (M82[:, 0], Pff, label =  "free-free emission")

plt.xlabel(r"Wavelength ($\mu m$)", size=20)
plt.ylabel(r"Luminosity ($ L_\odot$)", size=20)

plt.xlim([np.min(M82[:, 0]), np.max(M82[:, 0])])
plt.ylim(ymin=1E-8)
plt.legend()
plt.show()