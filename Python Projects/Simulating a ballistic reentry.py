#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Simulating a ballistic Reentry 
#Author: Sean Zimmari
#Created on 20th September 2023

import numpy as np #For numerical calcuations
import matplotlib.pyplot as plt #For Plotting
import scipy.integrate as sci #For integration


# In[16]:


#Assumed data:
he = 122000.0 #Entry altitude
ha = np.linspace(0,122000,100000) #Height above the ground
R = 6378000 #Radius of the Earth in meters
Ve = 7500.0 #Initial velocity in m/s
beta = 0.1354/1000.0 #Density constant in 1/m
rhos = 1.225 #Air density at sealevel in kg/m^3
m = 1350.0 #Mass of the spacecraft in kg
CD = 1.5 #Drag coefficient
S = 2.8 #Planform area of entry vehicle in m^2
gammae = -2.0*np.pi/180.0 #Entry flight path angle
gs = 9.81

##Formula that relates the velocity of the entry vehicle with altitude, assuming a lift equals to zero
Va = Ve*np.exp((1.0/(2*beta))*(rhos/(np.sin(gammae)))*(S*CD/m)*np.exp(-beta*ha))

############ANALYTICAL PROCESS################

#Now we have to compute acceleration
dhdt = Va*np.sin(gammae) #This equations assumes that gammae is constant
dVdh = (Va[0:-2]-Va[1:-1])/(ha[0:-2]-ha[1:-1]) #I am considering a little variation of the velocity through altitude
#Now acceleration = dV/dt = (dV/dh)*(dh/dt):
accel = dVdh*dhdt[0:-2]
amax = -(beta*Ve**2)/(2*np.exp(1))*np.sin(gammae)
print(f'Maximum acceleration in Gs = {amax/gs}')

###########NUMERICAL PROCESS##################

def Derivatives(state,t):
    global gs, R
    ##State vector
    V = state[0]
    h = state[1]
    gamma = state[2]
    
    #Gravity model
    g = gs*(R/(R+h))**2
    r = R + h
    
    #Density model
    rho = rhos*np.exp(-beta*h) 
    
    #Aerodynamic model : assuming L = 0
    L = 0.0
    D = 0.5*rho*V**2*S*CD
    
    
    
    #Equations for planar flight
    dVdt = -D/m - (g*np.sin(gamma))
    dgammadt = (L/m - (g-(V**2)/r)*np.cos(gamma))/V
    dhdt = V*np.sin(gamma)
    
    #Return statedot
    statedot = np.asarray([dVdt,dhdt,dgammadt])
    return statedot


#So now we have to integrate in fucntion of time:
stateinitial = np.asarray([Ve,he,gammae])
tout = np.linspace(0,450,10000)
stateout = sci.odeint(Derivatives,stateinitial,tout)
Vout = stateout[:,0]
hout = stateout[:,1]
gammaout = stateout[:,2]

accelout = (Vout[0:-2]-Vout[1:-1])/(tout[0:-2]-tout[1:-1])


# In[17]:


#Plotting:

#Plot of Va through ha:
plt.figure()
plt.plot(ha,Va,label = 'Analytical')
plt.plot(hout,Vout,label = 'Numerical')
plt.xlabel('Altitude (m)')
plt.ylabel('Velocity (m/s)')
plt.grid()
plt.legend()

#Plot of acceleration through time:
plt.figure()
plt.plot(ha[0:-2],accel/gs,label = 'Analytical')
plt.plot(hout[0:-2],accelout/gs,label='Numerical')
plt.xlabel('Altitude (m)')
plt.ylabel("G's")
plt.grid()
plt.legend()

#Plot gamma:
plt.figure()
plt.plot(tout,gammaout,'b-')
plt.xlabel('Time (sec)')
plt.ylabel('Flight Path Angle (rad)')
plt.grid()

#Plot altitude with time:
plt.figure()
plt.plot(tout,hout,'b-')
plt.xlabel('Time (sec)')
plt.ylabel('Altitude (m)')
plt.grid()


# In[ ]:




