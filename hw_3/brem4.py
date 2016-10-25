#!/usr/bin/python
#first built by Samantha Chang-Chun Chen 10/23/2016
# Class PHYS 239
### HW3

import numpy as np
from astropy.constants import a0, e, m_e, c
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

#### ===================  ####

e = e.gauss.value
a0 = a0.cgs.value
m_e = m_e.cgs.value
####============######
Z = 79
x0 = 8000
Zee = Z * e**2

dt = 1E-13
###########
n = 100

m_f = np.full([n,n], np.nan)

for i in range (n):
    y0 = np.logspace(np.log(1E2),np.log(1E4),100)
    y0 = y0[i]   
    for j in range(n):
        v0 = np.logspace(np.log(-1E6),np.log(-1E9),100)
        v0 = v0[j]
        #### ==== at t=0  ======= ########
        p = [np.array([x0,y0])* a0 ]            ##position(x,y)
        r = [np.sqrt(np.sum(p[0]**2))]
        v = [np.array([v0, 0])]
        a = [-Zee/ m_e / r[0]**3 * p[0]]
        t = [0]


        for k in range (1, 50):
            p.append(p[k-1] + v[k-1] * dt + a[k-1] * dt**2 / 2)
            r.append(np.sqrt(np.sum(p[k]**2)))
            a.append(-Zee/m_e/r[k]**3 * p[k])
            v.append(v[k-1] + a[k-1] * dt)
            t.append( t[k-1] + dt)
        
        a = np.array(a)         
        t = np.array(t)

        ay_omaga = np.fft.fftshift(np.abs(np.fft.fft(a[:,1])))
        ax_omaga = np.abs(np.fft.fftshift(np.fft.fft(a[:,0])))
        freq = np.fft.fftshift(np.fft.fftfreq(len(t), dt))
        a_omaga = np.sqrt(ax_omaga**2 + ay_omaga**2)
    
        m_f[i,j] = np.abs(freq[np.argmax(a_omaga)])

          

plt.figure()
plt.imshow(np.log(m_f), origin='lower',extent=[np.log(-1E6),np.log(-1E9),4,5], aspect="auto")
plt.xlabel(r'Initial x-velocity $log_{10}(v_0 (cm/s))$', size=16)
plt.ylabel(r'Initial y-position $log_{10}(b / a_0)$', size=16)
plt.title(r'Peak frequencies ($log_{10}(f/Hz)$)', size=20)
plt.colorbar()

plt.show()
