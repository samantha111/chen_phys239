#!/usr/bin/python
#first built by Samantha Chang-Chun Chen 10/23/2016
# Class PHYS 239
### HW3

import numpy as np
from astropy.constants import a0, e, m_e, c
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


#### ===CGS=======  ####
e = e.gauss.value
a0 = a0.cgs.value
m_e = m_e.cgs.value
####==Inputs==========######
Z = 100
x0 = 1E6
y0 = 1E4
v0 = -5E7
Zee = Z * e**2
niter = 200
dt = 1E-12


###### ==== at t=0  ======= ########

p = [np.array([x0,y0])* a0 ]            ##position(x,y)
r = [np.sqrt(np.sum(p[0]**2))]      ##distance (cm)
v = [np.array([v0, 0])]
a = [np.array(-Zee/ m_e / r[0]**3 * p[0])]        ###acceleration (ax,ay)
t = [0]

###########


for i in range (1, niter):
    p.append(p[i-1] + v[i-1] * dt + a[i-1] * dt**2 / 2)
    r.append(np.sqrt(np.sum(p[i]**2)))
    v.append(v[i-1] + a[i-1] * dt)
    a.append(-Zee/m_e/r[i]**3 * p[i])
    t.append( t[i-1] + dt)
##########################
p = np.array(p) / a0     # Converting back to a0
r = np.array(r) / a0     # Converting back to a0
v = np.array(v)          # Converting back to cm/s
a = np.array(a)          # Converting back to cm/s**2
t = np.array(t)

xy=['x','y']



######  x-y plot #####

plt.figure()

plt.scatter(p[:, 0], p[:, 1],color="orange")
plt.xlim([np.min(p[:, 0]), np.max(p[:, 0])])
plt.ylim([np.min(p[:, 1]), np.max(p[:, 1])])
#plt.axis([np.min(p[:, 0]), np.max(p[:, 0]),np.min(p[:, 1]), np.max(p[:, 1])])
plt.axis('equal')
plt.xlabel(r'$x$ ($a_0$)', size=16)
plt.ylabel(r'$y$ ($a_0$)', size=16)
plt.annotate('starting point', xy = (x0,y0), xytext=(x0, y0+2E5)
              ,xycoords ='data', arrowprops=dict(color='red',
               arrowstyle="-|>" ))
plt.title('x-y plot', size=20)

plt.show()


######velocity-time plot################
plt.figure()
for i in range(2):

    plt.subplot(1, 2, i + 1)
    plt.plot(t, v[:, i])
    plt.xlim([0, np.max(t)])
    plt.ylim([np.min(v[:, i]), np.max(v[:, i])])
    plt.xlabel(r'$t$ ($s$)', size=16)
    plt.ylabel(r'$v$ ($cm/s$)', size=16)
    plt.title(xy[i]+'-velocity', size=20)
plt.show()

#### t - acceleration plot ################
plt.figure()
for i in range(2):

    plt.subplot(1, 2, i + 1)
    plt.plot(t, a[:, i])
    plt.xlim([0, np.max(t)])
    plt.ylim([np.min(a[:, i]), np.max(a[:, i])])
    plt.xlabel(r'$t$ ($s$)', size=16)
    plt.ylabel(r'$a$ ($cm/s^2$)', size=16)
    plt.title(xy[i]+'-Acceleration', size=20)
plt.show()
########################
plt.figure()
ay_omaga = np.fft.fftshift(np.abs(np.fft.fft(a[:,1])))
ax_omaga = np.abs(np.fft.fftshift(np.fft.fft(a[:,0])))
freq = np.fft.fftshift(np.fft.fftfreq(len(t), dt))
a_omaga = np.sqrt(ax_omaga**2 + ay_omaga**2)


l = int(len(freq) // 2)

freq = freq[l:]
a_omaga = a_omaga[l:]
ax_omaga = ax_omaga[l:]
ay_omaga = ay_omaga[l:]

plt.plot(freq, ax_omaga,color='red',label=r'$a_x(\omega)$')
plt.plot(freq, ay_omaga,color='orange',label=r'$a_y(\omega)$')
plt.plot(freq, a_omaga,color='black',label=r'$a_{total}(\omega)$')
plt.axvline(x=5E10,color='black',linestyle='--')
plt.annotate('f=5E10', xy = (5E10,0), xytext=(6E10,3E18)
              ,xycoords ='data', arrowprops=dict(color='black',
               arrowstyle="-|>" ))

plt.ylabel(r'$a(\omega)$',size=14)
plt.xlabel(r'frequency(Hz)')
       
plt.legend()
plt.show()

m_f = np.abs(freq[np.argmax(a_omaga)])
print (np.argmax(a_omaga),m_f)

       
