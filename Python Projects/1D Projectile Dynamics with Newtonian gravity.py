#!/usr/bin/env python
# coding: utf-8

# In[1]:


#1D Projectile Dynamics with newtonian gravity
#Author: Sean Zimmari
#Created on 20th September 2023

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sci


# In[15]:


#Equations of motions
#F = m*a = m*zddot
#z is the altitude of the surface ( in meters )
#zdot is the velocity
#zddot is the acceleration
#Differential equation

#Define some constant parameters
##rocket
mass = 640.0/1000.0 ##rocket mass in kilograms

#Planet Earth
G = 6.6742*10**-11 #Gravitational constant
Mplanet = 5.972e24 #in kilograms
Rplanet = 6357000 #in meters

#Planet Kerbin
#Rplanet = 600000 #in meters
#Mplanet = 5.2915158*10**22 #in kilograms

#Definition of a function which represents gravitional acceleration
def gravity(z):
    global Rplanet,Mplanet
    
    r = np.sqrt(z**2) #norm of the vector of the center of the planet: it represents the distance from the centre of the earth
    
    if r < Rplanet:
        accel = 0.0
    else:
        accel = G*Mplanet/(r**3)*r
        
    print(accel)
        
    return accel
    
#Test Surface gravity
print('Surface gravity (m/s^2) = ', gravity(Rplanet))

def Derivatives(state,t):
    #Global variables
    global mass
    #state vector
    z = state[0]
    velz = state[1]
    
    #Compute zdot
    zdot = velz
    
    #Compute the total forces
    ##Gravity
    gravityF = -gravity(z)*mass
    
    ##Aerodynamics
    aeroF = 0.0
    
    ##Thrust
    thrustF = 0.0
    
    Forces = gravityF + aeroF + thrustF
    
    #Compute acceleration
    zddot = Forces/mass
    
    #Compute the statedot
    statedot = np.asarray([zdot,zddot])
    
    return statedot


# In[16]:


#Main script

#Initial conditions
z0 = Rplanet
velz0 = 33*331.0 #m/s
stateinitial = np.asarray([z0,velz0])

#Time window
tout = np.linspace(0,33.5,1000) #I have assumed that the rocket stay in flight for 30 seconds
stateout = sci.odeint(Derivatives,stateinitial,tout)

#Rename variables
zout = stateout[:,0]
velzout = stateout[:,1]

#Plot altitude
plt.plot(tout,zout)
plt.xlabel('Time(sec)')
plt.ylabel('Altitude(m)')
plt.grid()

#Plot velocity
plt.figure()
plt.plot(tout,velzout)
plt.xlabel('Time(sec)')
plt.ylabel('Normal speed (m/s)')
plt.grid()


# In[ ]:




