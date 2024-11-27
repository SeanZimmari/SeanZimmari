#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Simulating a PulseJet thruster
#Author: Sean Zimmari
#Created on 27/09/2023
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sci


# In[59]:


#We have to solve a system of differential equations
pulse_duration = 40.0/1000.0
pulse_wait_duration = 60.0/1000.0

thrusters = 0.0

pulse_off_time = -pulse_wait_duration #Because I want to fire right away
pulse_on_time = 0.0

Np = 0.0 #Number of pulses

def Thrusters(t):
    global thrusters, pulse_off_time,pulse_on_time,pulse_wait_duration,pulse_duration,Np
    
    #Torque on the spacecraft
    ##Model Pulse Jet Thruster
    ###Are thrusters currently on?
    if thrusters == 1:
        #Thrusters are currently on
        ##Do we need to turn them off?
        if t >= pulse_on_time + pulse_duration:
            thrusters = 0
            pulse_off_time = t
    else:
            #Thrusters are currently off
            #Can I turn them back on?
        if t >= pulse_off_time + pulse_wait_duration:
            thrusters = 1
            pulse_on_time = t
            Np += 1
    
    r = 1.0
    F = 50.0*thrusters
    Torque = np.asarray([0,r*F,0])
    return Torque

def Derivatives(state,t,Torque):
    global thrusters, pulse_off_time,pulse_on_time,pulse_wait_duration,pulse_duration,Np
    #Unpack my state vector
    rpy = state[3:6]
    w = state[0:3]
    ##Inertia Matrix I
    Ixx = 800.0
    Iyy = Ixx
    Izz = 500.0
    Inertia = np.asarray([[Ixx,0,0],[0,Iyy,0],[0,0,Izz]])
    Inertia_inverse = np.asarray([[1/Ixx,0,0],[0,1/Iyy,0],[0,0,1/Izz]])#Assuming a symmetric spacecraft
    
    #Torque on the spacecraft
    ##Model Pulse Jet Thruster
    ###Are thrusters currently on?
    if thrusters == 1:
        #Thrusters are currently on
        ##Do we need to turn them off?
        if t >= pulse_on_time + pulse_duration:
            thrusters = 0
            pulse_off_time = t
    else:
            #Thrusters are currently off
            #Can I turn them back on?
        if t >= pulse_off_time + pulse_wait_duration:
            thrusters = 1
            pulse_on_time = t
            Np += 1
    
    r = 1.0
    F = 50.0*thrusters
    Torque = np.asarray([0,r*F,0])
    
    dMdt = Torque - np.cross(w,np.matmul(Inertia,w)) #.cross is for applying vectorial product in there, instead .matmul 
    #is used to represent the scalar product between I and w
    dwdt = np.matmul(Inertia_inverse,dMdt) #Remember that alfa, which is the angular acceleration,
    #equals to Torque divided by Inertia
    
    #Kinematics
    phi = rpy[0]
    theta = rpy[1]
    psi = rpy[2]
    ctheta = np.cos(theta)
    cpsi = np.cos(psi)
    sphi = np.sin(phi)
    stheta = np.sin(theta)
    cphi = np.cos(phi)
    spsi = np.sin(psi)
    ttheta = np.tan(theta)
    J = np.asarray([[1.0,sphi*ttheta,cphi*ttheta],[0.0,cphi,-sphi],[0.0,sphi/ctheta,cphi/ctheta]])
    drpydt = np.matmul(J,w)
    return np.hstack([dwdt,drpydt])

##So now we have to solve Derivatives
tout = np.linspace(0,1000,10000)
dt = tout[2] - tout[1]
stateout = np.zeros((6,len(tout)))

stateout[:,0] = np.asarray([0,0,0*100.0/500.0,0,0,0])


for ctr in range(0,len(tout)-1):
    ti = tout[ctr]
    statei = stateout[:,ctr]
    Torque = Thrusters(ti)
    dstate1 = Derivatives(statei,ti,Torque)
    dstate2 = Derivatives(statei+dstate1*dt/2.0,ti+dt/2.0,Torque)
    dstate3 = Derivatives(statei+dstate2*dt/2.0,ti+dt/2.0,Torque)
    dstate4 = Derivatives(statei+dstate3*dt,ti+dt,Torque)
    chi = (1.0/6.0)*(dstate1 + 2*dstate2 + 2*dstate3 + dstate4)
    stateout[:,ctr+1] = statei + dt*chi
    dwdt = Derivatives(statei,ti,Torque)

    wout = stateout[0:2,:]
    rp = stateout[3:5,:]

plt.figure()
plt.plot(tout,np.transpose(rp))
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Roll,Pitch (rad)')

plt.figure()
plt.plot(tout,np.transpose(wout))
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Angular Velocity (rad/s)')


# In[ ]:




