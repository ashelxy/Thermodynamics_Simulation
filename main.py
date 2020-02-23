# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 15:06:34 2019

@author: xl7817
"""

import numpy as np
import scipy as sp
import matplotlib as plt

import ballclass as bc
import simulation as sim

#%%

C = sim.Simulation(729)
C.run(0, True)
#C.fakerun(1500)
C.hist()
C.boltzmann()
C.frameke()
C.temp()
C.pressure()

#%%

"""

Attempt to obtain results to plot points of pressure against temperature 
with fixed N and V.

For the graphs, a loop of 10 cycles is averaged to form a point on the graph. 

"""

pressure = []
#temperature = []

for x in range(10):
    C = sim.Simulation(10)
#    C.run(0, True)
    C.fakerun(500)
    
#C.hist()
    
#C.boltzmann()
#C.frameke()
#    temperature.append(C.temp())
    pressure.append(C.pressure())
    
#print(np.mean(temperature))
print(np.mean(pressure))



#%%

C = sim.Simulation(729)
C.run(0, True)

#%%

# PT graph:

t = [306.77991,207.3144013,281.9972559,294.1334831,273.4220566,11.48852354,
     13.13291499,12.03043392,120.6665597,106.3788494,103.1505164,93.74892372,
     45.48278038,183.5918577,48.37331904,174.0677235,206.2801366,301.9714243,
     297.0743021,431.4332079,757.8735523,464.3895315,428.9088039,521.4033625,
     798.2360558,500.6976271,562.4568268,556.396558,738.6950039,731.7697466,
     836.9548505,702.4035573]

p = [7.229418e-20,4.707599e-20,4.405013e-20,5.017662e-20,6.977922e-20,
     3.431249e-21,2.813135e-21,1.388555e-21,2.645525e-20,2.089728e-20,
     5.103174e-20,1.259067e-20,1.735941e-20,6.735791e-20,2.984781e-20,
     8.048491e-20,9.614006e-20,1.025727e-19,4.871387e-20,6.546951e-20,
     2.021311e-19,1.086804e-19,1.150370e-19,1.047300e-19,1.940927e-19,
     1.124290e-19,1.668152e-19,1.080624e-19,2.635930e-19,2.916165e-19,
     2.542162e-19,2.409108e-19]

fit,cov = sp.polyfit(t,p,1,cov=True)
fit_values = sp.poly1d(fit)

print('Slope = %.3e +/- %.3e' %(fit[0],sp.sqrt(cov[0,0])))
print('Intercept = %.3e +/- %.3e' %(fit[1],sp.sqrt(cov[1,1])))

plt.plot(t, p, '.')
plt.plot(t,fit_values(t), color = "black")
plt.xlabel("Temperature (K)")
plt.ylabel("Pressure (F/m)")
plt.savefig("pt.png")

#%%

# PV graph:

v = [7.068583470577E+00,1.256637061436E+01,1.963495408494E+01,
     2.827433388231E+01,3.848451000648E+01,5.026548245744E+01,
     6.361725123519E+01,7.853981633974E+01,9.503317777109E+01,
     1.130973355292E+02,1.327322896142E+02,1.539380400259E+02,
     1.767145867644E+02,2.010619298297E+02,2.269800692219E+02,
     2.544690049408E+02,2.835287369865E+02,3.141592653590E+02]


p1 = [6.561234708336E-19,1.470230493397E-19,4.639783491328E-20,
      2.837127475169E-20,1.717892570044E-20,1.254453741785E-20,
      1.113583945332E-20,8.833621373798E-21,7.299228517247E-21,
      6.064996196062E-21,5.150220400653E-21,4.761982233141E-21,
      3.763464621949E-21,2.673451297785E-21,2.375039683338E-21,
      1.726831308932E-21,1.661881574574E-21,1.616552008675E-21]

plt.plot(v,p1, '.')
plt.xlabel("Volume (m^2)")
plt.ylabel("Pressure (F/m)")
plt.savefig("pv.png")
























































































