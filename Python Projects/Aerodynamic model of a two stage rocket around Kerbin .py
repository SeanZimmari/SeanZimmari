#!/usr/bin/env python
# coding: utf-8

# In[2]:


#Numerical integration of a 2d orbit
#Author: Sean Zimmari
#Created on 20th September 2023

import numpy as np #For numerical calcuations
import matplotlib.pyplot as plt #For Plotting
import scipy.integrate as sci #For integration


# In[7]:


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
#name = 'Earth'

#Planet Kerbin
Rplanet = 600000 #in meters
Mplanet = 5.2915158*10**22 #in kilograms
name = 'Kerbin'

#Initial conditions for two stage rocket
x0 = Rplanet
z0 = 0.0
velz0 = 0.0
velx0 = 0.0
r0 = 200000 + Rplanet
period = 6000 #2*np.pi/np.sqrt(G*Mplanet)*r0**(3.0/2.0)*1.5 #I multiplied 1.5 because the equation of period assume the complete simmetry between axis: if I have an eliptical orbit, then the period will be bigger
tsep1 = 2.0 #Time to remove the first stage
mass1tons = 0.2
mass1 = mass1tons*2000/2.2 #Mass of the first stage in kilograms
ve = Isp1 * 9.81
t2start = 261.0
t2end = t2start + 20.0
D = 1.141 #Diameter of the rocket's nozzle
CD = 0.4 #Drag coefficient

#Create and Aerodynamic class
class Aerodynamics():
    
    def __init__(self,name):
        self.name = name
        if name == 'Kerbin':
            #Import aero model for interpolation
            data = np.loadtxt('Kerbin_data_for_aerodynamics')
            self.altitude = data[:,0]
            print(self.altitude)
            self.density = data[:,3]
            print(self.density)
            self.rhos = self.density[0]
            self.beta = 0.0
        elif name == 'Earth':
            #going to use the Earth aero model
            self.beta = 0.1354/1000.0 #Density constant [m-1]
            self.rhos = 1.225 #kg/m3
    
    def getDensity(self,altitude):
        if self.name == 'Kerbin':
            ##Interpolation
            rho = np.interp(altitude,self.altitude,self.density)
        elif self.name == 'Earth':
            ##Equation
            rho = self.rhos*np.exp(-self.beta*altitude)
        return rho

#Create the aerodynamics variable which is an instance of the aerodynamics class
aeroModel = Aerodynamics(name)

#Plot air density as a function of altitude
test_altitude = np.linspace(0,100000,100)
test_rho = aeroModel.getDensity(test_altitude)
plt.figure()
plt.plot(test_altitude,test_rho,'b-')
plt.xlabel('Altitude(m)')
plt.ylabel('Air Density(kg/m^3)')
plt.grid()

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
        
    return np.asarray([accelx,accelz]),r
    
#Test Surface gravity
#print('Surface gravity (m/s^2) = ', gravity(0,Rplanet))

def propulsion(t):
    global max_thrust,Isp,tMECO, ve
    ##Timing for thrusters
    theta = 10.0*np.pi/180.0
    if t < tMECO: #Mean Engine Cut-Off
        thrustF = max_thrust #during the lift-off
        mdot = -thrustF/ve #mass lost
    if t > tMECO and t < (tMECO + tsep1):
        thrustF = 0.0
        #masslost = mass1
        mdot = -mass1/tsep1 #Because it is a small area created by the rectangle y=mass and x=time
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
    global aeroModel, Rplanet
    
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
    accel,r = gravity(x,z)
    gravityF = -accel*mass
    
    ##Aerodynamics
    altitude = r - Rplanet #Altitude above the surface
    rho = aeroModel.getDensity(altitude) #Air Density
    V = np.sqrt(velx**2 + velz**2)
    qinf = (np.pi/8.0)*rho*D**2*abs(V) #Dynamic Pressure
    aeroF = -qinf*CD*np.asarray([velx,velz]) #Draf force
    
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
        


# In[8]:


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
tout = np.linspace(0,period,1000) #I have assumed that the rocket stay in flight for the period
stateout = sci.odeint(Derivatives,stateinitial,tout) #Solve ordinary differential equations

#The result of this method will give an array with 5 elements that have one less differential order
#So for example if the first two elements of 'Derivatives' represent the velocity of the spacecraft
#Then the resulting first two elements will represent the space travelled by the spacecraft

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
#Plotting the orbit
plt.plot(xout,zout,'r-', label='Orbit')
plt.plot(xout[0],zout[0],'g*')
theta = np.linspace(0,2*np.pi,1000)
xplanet = Rplanet*np.sin(theta)
yplanet = Rplanet*np.cos(theta)
plt.plot(xplanet,yplanet,'b-',label='Planet')
plt.grid()
plt.legend()


# In[ ]:




