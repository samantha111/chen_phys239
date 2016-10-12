#!/usr/bin/python
# First bulit by Samantha Chang-Chun Chen 10/08/2016
# Class PHYS 239

import numpy as np
import astropy

##==========###
D = 100*3.09*10**18
n = 1


####Problem 1#######################
### what's is the column density?####
#################### 
print ('Problem 1:')
N = n*D
print ('Column density is N=%e 1/cm^2' %(N))


d = {10**-3:3.33*10**-24,1:3.33*10**-21,10**3:3.33*10**-18}
for optdep in d:
    crsec = d[optdep]
    print 'Cross section = %s when optical depth is %s' % (crsec,optdep)


###Problem 2##########
###### 
### seperate D into 100 steps##
######absorption######### 
I0= 1
n = 1
crsec = 3*10**-22  ##cross section [cm**2]
ds = 1*3*10**18 ##[cm]
x = crsec*n*ds


count = 0
A = 1
while count < 100:
    A = A*np.exp(-x)
    count +=1
    if count >=100:
        break
I = I0*A
print ('Problem 2:')
print "I at s=D is %.2f" %(I)

####emission##########
S0 = 1
n = 1
crsec = 3*10**-22  ## [cm**2]
ds = 1*3*10**18 ##[cm]
x = crsec*n*ds

T = np.exp(-x) + 1
count = 1
while count < 100: 
    T = T*np.exp(-x) + 1
    count +=1 
    if count >=100: 
        break
S = S0 * T
print "S is %.2f" %(S) 
intensity = I + S 

print "With cross section %s (optically thin), initial intensity I0=%s, and source funtion S0=%.2f, the specific intensity at s=D is I+S=%.2f" %(crsec,I0,S0,intensity) 

###====================###







