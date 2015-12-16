import numpy as np

import pandas as pd

import scipy, statsmodels

# something like a mid-latitude seasonal cycle
y = -5+30*np.sin(np.linspace(0,np.pi,365))

# add some noise
Ta = y+scipy.stats.norm.rvs(0,2.5,len(y))

# get a 10 day moving average
Tg = pd.rolling_mean(Ta,10)

from hunt_resp import mresp

# coefficients for broad leaf trees
r0 = 1.2855
r1 = 0.2061
r2 = 0.0402
nl = 1.5

# initialize array to hold respiration output
rdarkC = np.zeros(len(Ta))

for i in range(0,len(Ta)):
  rdarkC[i] = mresp(r0,r1,r2,nl,Ta[i],Tg[i])

from acme_resp import amresp

# acme constants, could make a function of seasonal cycle
lnc = 1 # mean value from TRY (Kattge, et al. 2011)
elai = -np.log(0.7)/0.3 # to set nscaler = 1 
# could use 4, a very coarse estimate from Harvard Bigfoot project

rdarkA = np.zeros(len(Ta))

for i in range(0,len(Ta)):
  rdarkA[i] = amresp(Ta[i]+273.15,lnc,elai)

#----------Plotting----------

import matplotlib.pyplot as plt

x = np.linspace(1,365,365)

# units? m^-2, but does LAI address this?
plt.plot(x,rdarkC,'k-',x,rdarkA,'b-',x,rdarkN0,'g-')
plt.xlabel('Day of Year')
plt.ylabel('Maintentance Respiration [umolCO2/(m^2 s)]')
plt.legend(['Huntingford,et al.','ALM,N=1.59,elai=4','ALM,N=1,nscal=1'])
plt.axis([1,365,0,1.25])
plt.show()
#plt.savefig("resp_noNnoScal.pdf")
