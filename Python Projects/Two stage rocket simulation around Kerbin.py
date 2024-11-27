#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Numerical integration of a 2d orbit
#Author: Sean Zimmari
#Created on 20th September 2023

import numpy as np #For numerical calcuations
import matplotlib.pyplot as plt #For Plotting
import scipy.integrate as sci #For integration


# In[14]:


#Equations of motions
#F = m*a = m*zddot
#z is the altitude from the center of the planet along the north pole
#x is the altitude from the center of the earth along the equator
#zdot is the velocity along z
#zddot is the acceleration along z
#Differential equation

#Define some constant parameters
##Rocket : Conditions for sub-orbital flight around Kerbin 
weighttons = 5.3
mass0 = weighttons*2000/2.2 ##rocket mass in kilograms
max_thrust = 167970.0 #Newton
Isp1 = 250.0
Isp2 = 400.0
tMECO = 38.0 #Main Engine Cutoff Time : this means that I am going to burn fuel for 20 seconds

#Planet Earth
G = 6.6742*10**-11 #Gravitational constant
#Mplanet = 5.972e24 #in kilograms
#Rplanet = 6357000 #in meters

#Planet Kerbin
Rplanet = 600000 #in meters
Mplanet = 5.2915158*10**22 #in kilograms

#Initial conditions for single stage rocket
x0 = Rplanet
z0 = 0.0
velz0 = 0.0
velx0 = 0.0
r0 = 200000 + Rplanet
period = 6000 #2*np.pi/np.sqrt(G*Mplanet)*r0**(3.0/2.0)*1.5 #I multiplied 1.5 because the equation of period assume the complete simmetry between axis: if I have an eliptical orbit, then the period will be bigger
tsep1 = 2.0 #Time to remove the first stage
mass1tons = 0.2
mass1 = mass1tons*2000/2.2
ve = Isp1 * 9.81
t2start = 261.0
t2end = t2start + 20.0

#Definition of a function which represents gravitional acceleration
def gravity(x,z):
    global Rplanet,Mplanet
    
    r = np.sqrt(x**2 + z**2) #norm of the vector of the center of the planet: it represents the distance from the centre of the earth
    
    if r < Rplanet:
        accelx = 0.0
        accelz = 0.0
    else:
        accelx = G*Mplanet/(r**3)*x
        accelz = G*Mplanet/(r**3)*z
        
        
    return np.asarray([accelx,accelz])
    
#Test Surface gravity
print('Surface gravity (m/s^2) = ', gravity(0,Rplanet))

def propulsion(t):
    global max_thrust,Isp,tMECO, ve
    ##Timing for thrusters
    theta = 10.0*np.pi/180.0
    if t < tMECO: #Mean Engine Cut-Off
        thrustF = max_thrust #during the lift-off
        mdot = -thrustF/ve
    if t > tMECO and t < (tMECO + tsep1):
        thrustF = 0.0
        #masslost = mass1
        mdot = -mass1/tsep1
    if t > (tMECO + tsep1):
        thrustF = 0.0
        mdot = 0.0
    if t > t2start and t < t2end:
        theta = 90.0*np.pi/180.0  #Because I want to fire sideways
        ve = Isp2*9.81
        thrustF = max_thrust
        mdot = -thrustF/ve
    if t > t2end:
        theta = 0.0
        thrustF = 0.0
        mdot = 0.0
        
    #Angle of my thrust
    thrustx = thrustF*np.cos(theta)
    thrustz = thrustF*np.sin(theta)
    
    
    return np.asarray([thrustx,thrustz]),mdot

def Derivatives(state,t):
    
    #state vector
    x = state[0]
    z = state[1]
    velx = state[2]
    velz = state[3]
    mass = state[4]
    
    #Compute zdot
    zdot = velz
    xdot = velx
    
    #Compute the total forces
    ##Gravity
    gravityF = -gravity(x,z)*mass
    
    ##Aerodynamics
    aeroF = np.asarray([0.0,0.0])
    
    ##Thrust
    thrustF,mdot = propulsion(t)
    
    Forces = gravityF + aeroF + thrustF
    
    if mass > 0:
        #Compute acceleration
        ddot = Forces/mass
    else:
        ddot = 0.0
        mass = 0.0
    
    
    
    #Compute the statedot
    statedot = np.asarray([xdot,zdot,ddot[0],ddot[1],mdot])
    
    return statedot


# In[15]:


#Main script

#Initial conditions - For orbit
#x0 = Rplanet #Because I need to stay outside the earth's surface
#z0 = 0.0
#r0 = np.sqrt(x0**2 + z0**2)
#velz0 = np.sqrt(G*Mplanet/r0)*1.1 #m/s
#velx0 = 100.0 
#stateinitial = np.asarray([x0,z0,velx0,velz0])

#Initial condition for single-stage rocket
stateinitial = np.asarray([x0,z0,velx0,velz0,mass0])


#Time window
tout = np.linspace(0,period,1000) #I have assumed that the rocket stay in flight for 30 seconds
stateout = sci.odeint(Derivatives,stateinitial,tout)

#Rename variables
xout = stateout[:,0]
zout = stateout[:,1]
altitude = np.sqrt(xout**2+zout**2) - Rplanet #Altitude from the earth's surface
velxout = stateout[:,2]
velzout = stateout[:,3]
velout = np.sqrt(velxout**2 + velzout**2)
massout = stateout[:,4]

#Plot altitude
plt.plot(tout,altitude)
plt.xlabel('Time(sec)')
plt.ylabel('Altitude(m)')
plt.grid()

#Plot velocity
plt.figure()
plt.plot(tout,velout)
plt.xlabel('Time(sec)')
plt.ylabel('Total speed (m/s)')
plt.grid()

#Plot mass
plt.figure()
plt.plot(tout,massout)
plt.xlabel('Time (sec)')
plt.ylabel('Mass (kg)')
plt.grid()

#Plot 2D Orbit
plt.figure()
plt.plot(xout,zout,'r-', label='Orbit')
plt.plot(xout[0],zout[0],'g*')
theta = np.linspace(0,2*np.pi,1000)
xplanet = Rplanet*np.sin(theta)
yplanet = Rplanet*np.cos(theta)
plt.plot(xplanet,yplanet,'b-',label='Planet')
plt.grid()
plt.legend()


# In[ ]:




