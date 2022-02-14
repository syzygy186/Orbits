#!/usr/bin/env python3

import math
import matplotlib.pyplot as plt
import numpy as np 
import functools

G = 4*np.pi**2

# Time is measured in years
# Distance is Measured in AU
# Mass is measured in MSun



class Vector:

    '''
    3-Dimensional Vector
    
    Holds 3 numbers and all the operators associated 
    with the algebraic manipulations of Vectors in 3D
    
    '''
    
    # Declare a Vector
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
    
    # Addition
    def __add__(self,other):
        return Vector(self.x+other.x,self.y+other.y,self.z+other.z)
    
    # Subtraction
    def __sub__(self,other):
        return Vector(self.x-other.x,self.y-other.y,self.z-other.z)
    
    # Inner Product
    @functools.singledispatchmethod
    def __mul__(self, other):
        return (self.x*other.x + self.y*other.y+self.z*other.z)
    
    # Multiplication by Scalar
    @__mul__.register(float)
    def _(self,other):
        return Vector(other*self.x,other*self.y,other*self.z)
    
    # Reversed Inner Product
    @functools.singledispatchmethod
    def __rmul__(self, other):
        return (self.x*other.x + self.y*other.y+self.z*other.z)
    
    # Reversed Scalar Multiplication
    @__rmul__.register(float)
    def _(self,other):
        return Vector(self.x*other,self.y*other,self.z*other)
    
    # Cross Product
    def __matmul__(self,other):
        return Vector(self.y*other.z-self.z*other.y,self.z*other.x-self.x*other.z,self.x*other.y-self.y*other.x)
    
    # Print Statement
    def __str__(self):
        return "("+str(self.x)+","+str(self.y)+","+str(self.z)+")"

class Particle:
    
    '''
    Particles for the Simulator
    They hold kinematical variables: Position, Velocity and Acceleration
    And the trajectory of the particle, the time should be managed externally.
    
    '''
    
    # Define a Partile: Name, Mass, Position and Velocity are Required, Acceleration is optional
    def __init__(self,Name,Mass,Position,Velocity,Acceleration=Vector(0,0,0)):
        self.Nam = Name
        self.Mas = Mass
        self.Pos = Position
        self.Vel = Velocity
        self.Acc = Acceleration
        self.Tra = [[],[],[]]
    
    def Evolve(self,dt,Omega=Vector(0,0,0)):
        
        # Euler's Method for Time Evolution
        self.Vel = self.Vel + self.Acc * dt
        self.Pos = self.Pos + self.Vel * dt
        
        # Record the Trajectory
        self.Tra[0].append(self.Pos.x)
        self.Tra[1].append(self.Pos.y)
        self.Tra[2].append(self.Pos.z)

class QParticle: 
    
    '''
    Clone of Particle but faster because Trajectories are not stored
    '''
    
    def __init__(self,Name,Mass,Position,Velocity,Acceleration=Vector(0,0,0)):
        self.Nam = Name
        self.Mas = Mass
        self.Pos = Position
        self.Vel = Velocity
        self.Acc = Acceleration
    
    def Evolve(self,dt,Omega=Vector(0,0,0)):
        
        # Euler's Method for Time Evolution
        self.Vel = self.Vel + self.Acc * dt
        self.Pos = self.Pos + self.Vel * dt
        
    def __str__(self):
        return self.Nam + " Position " + str(self.Pos)
            
class QSimulator:           
    def __init__(self,Particles,Omega=Vector(0,0,0)):
        self.Omega     = Omega
        self.Velocity  = Vector(0,0,0)
        self.Particles = Particles
    
    def GoToCM(self):
        R = Vector(0,0,0)
        P = Vector(0,0,0)
        M = 0
        for Particle in self.Particles:
            R = R + Particle.Mas*Particle.Pos
            P = P + Particle.Mas*Particle.Vel
            M = M + Particle.Mas
        
        self.Velocity = (1.0/M)*P
        self.Position = (1.0/M)*R
        
        for Particle in self.Particles:
            Particle.Vel = Particle.Vel - self.Velocity
            Particle.Pos = Particle.Pos - self.Position
    
    def GoToLab(self):
        for Particle in self.Particles:
            Particle.Vel = Particle.Vel + self.Velocity
    
    def GoToSynodic(self):
        
        self.GoToCM()
        ICM = 0
        L = Vector(0,0,0)
        for Particle in self.Particles:
            ICM += (Particle.Pos*Particle.Pos)*Particle.Mas
            L    = L + Particle.Mas * (Particle.Pos @ Particle.Vel)
        self.Omega = (1.0/ICM)*L
        
        for Particle in self.Particles:
            Particle.Vel = Particle.Vel - self.Omega @ Particle.Pos
        
        print("ICM : ",ICM)
        print("The Anagular Velocity of the Synodic system is:",self.Omega)
    
    def SetOmega(self,Omega):
        self.Omega = Omega
        for Particle in self.Particles:
            Particle.Vel = Particle.Vel - self.Omega @ Particle.Pos
    
    def Simulate(self,dt):
          
        # Interactions
        
        for i, Particle1 in enumerate(self.Particles):
            
            Particle1.Acc = Vector(0,0,0)
                
            # Coriolis Acceleration
            Particle1.Acc = Particle1.Acc - 2.0*(self.Omega @ Particle1.Vel)
                
            # Centrifugal Acceleration
            Particle1.Acc = Particle1.Acc - self.Omega @ (self.Omega @ Particle1.Pos)
                
                            
            for j, Particle2 in enumerate(self.Particles):
                if i == j: 
                    continue
                    
                g = G*Particle2.Mas
                r = Particle2.Pos - Particle1.Pos
                r3 = pow(r*r,1.5)
            
                Particle1.Acc.x += g*r.x/r3
                Particle1.Acc.y += g*r.y/r3
                Particle1.Acc.z += g*r.z/r3
                           
        # Time Evolution
        for Particle in self.Particles:
            Particle.Evolve(dt,self.Omega)
              
class Simulator:
    
    '''
    
    Deals with Time management and Interactions
    
    '''
    
    # Define a Simulator only needs particles, optionally a rotational frquency
    def __init__(self,Particles,Omega=Vector(0,0,0)):
        self.Omega     = Omega
        self.Velocity  = Vector(0,0,0)
        self.Particles = Particles
        self.Time      = []
        self.Energy    = []
        self.Angular   = []
        self.Potential = []
        self.Kinetic   = []
        self.Jacobi    = []
        self.Distance  = []
    
    # Go to the Center of mass
    def GoToCM(self):
        R = Vector(0,0,0)
        P = Vector(0,0,0)
        M = 0
        for Particle in self.Particles:
            R = R + Particle.Mas*Particle.Pos
            P = P + Particle.Mas*Particle.Vel
            M = M + Particle.Mas
        
        self.Velocity = (1.0/M)*P
        self.Position = (1.0/M)*R
        
        for Particle in self.Particles:
            Particle.Vel = Particle.Vel - self.Velocity
            Particle.Pos = Particle.Pos - self.Position
    
    # Go back to the system where things were defined
    def GoToLab(self):
        for Particle in self.Particles:
            Particle.Vel = Particle.Vel + self.Velocity
            Particle.Pos = Particle.Pos + self.Position
    
    
    def ComputeOmega(self):
        
        # L = ICM * Omega
        self.GoToCM()
        ICM = 0
        L = Vector(0,0,0)
        for Particle in self.Particles:
            ICM += (Particle.Pos*Particle.Pos)*Particle.Mas
            L    = L + Particle.Mas * (Particle.Pos @ Particle.Vel)
        self.Omega = (1.0/ICM)*L
        self.GoToLab()
    
    # Goes to Synodic, but Dynamics in the Synodic are unstable
    def GoToSynodic(self):
        
        self.GoToCM()
        ICM = 0
        L = Vector(0,0,0)
        for Particle in self.Particles:
            ICM += (Particle.Pos*Particle.Pos)*Particle.Mas
            L    = L + Particle.Mas * (Particle.Pos @ Particle.Vel)
        self.Omega = (1.0/ICM)*L
        
        for Particle in self.Particles:
            Particle.Vel = Particle.Vel - self.Omega @ Particle.Pos
        
        print("ICM : ",ICM)
        print("The Anagular Velocity of the Synodic system is:",self.Omega)
    
    # Do not use, only for debug
    def SetOmega(self,Omega):
        self.Omega = Omega
        for Particle in self.Particles:
            Particle.Vel = Particle.Vel - self.Omega @ Particle.Pos
    
    def Simulate(self,Time,Steps):
        
        t = 0
        dt = Time/Steps
        
        self.ComputeOmega()
        
        while t < Time: 
            t += dt
            k  = 0
            u  = 0
            lz = 0
            J  = 0
        
            # Interactions
            
            if len(self.Particles)>=2:
                D  = self.Particles[0].Pos - self.Particles[1].Pos
                D  = np.sqrt(D*D)
        
            for i, Particle1 in enumerate(self.Particles):
            
                
            
                if i == 2: 
                    J += (self.Omega*self.Omega)*(Particle1.Pos*Particle1.Pos)
                    r1 = self.Particles[0].Pos - Particle1.Pos
                    r2 = self.Particles[1].Pos - Particle1.Pos
                    J += 2*(self.Particles[0].Mas*G)/(np.sqrt(r1*r1))
                    J += 2*(self.Particles[1].Mas*G)/(np.sqrt(r2*r2))
                    Vel = Particle1.Vel - self.Omega @ Particle1.Pos
                    J += -(Vel*Vel)
            
                Particle1.Acc = Vector(0,0,0)
                
                # Coriolis Acceleration
                #Particle1.Acc = Particle1.Acc - 2.0*(self.Omega @ Particle1.Vel)
                
                # Centrifugal Acceleration
                #Particle1.Acc = Particle1.Acc - self.Omega @ (self.Omega @ Particle1.Pos)
                
                k += (1./2) * Particle1.Mas * (Particle1.Vel * Particle1.Vel)
                lz += Particle1.Mas * ((Particle1.Pos @ Particle1.Vel)).z
            
                for j, Particle2 in enumerate(self.Particles):
                    if i == j: 
                        continue
                    
                    g = G*Particle2.Mas
                    r = Particle2.Pos - Particle1.Pos
                    r3 = pow(r*r,1.5)
            
                    Particle1.Acc.x += g*r.x/r3
                    Particle1.Acc.y += g*r.y/r3
                    Particle1.Acc.z += g*r.z/r3
                
                    u += -g*Particle1.Mas/(2*np.sqrt(r*r))
           
            # Time Evolution
        
            for Particle in self.Particles:
                Particle.Evolve(dt,self.Omega)
        
            self.Energy.append(u+k)
            self.Kinetic.append(k)
            self.Potential.append(u)
            self.Time.append(t)
            self.Angular.append(lz)
            self.Jacobi.append(J)
            self.Distance.append(D)
    
    def GeneratePlots(self):
        plt.clf()
        plt.title("Orbits")
        plt.xlabel("AU")
        plt.ylabel("AU")
        for Particle in self.Particles:
            plt.plot(Particle.Tra[0],Particle.Tra[1],label=Particle.Nam)
        plt.legend(loc="best")
        plt.show()
        
        plt.clf()
        plt.title("Energies")
        plt.xlabel("yr")
        plt.ylabel("MSun AU^2 yr^-2")
        plt.plot(self.Time,self.Energy,label="Total Energy")
        #plt.plot(self.Time,self.Potential,label="Potential Energy")
        #plt.plot(self.Time,self.Kinetic,label="Kinetic Energy")
        plt.legend(loc="best")
        plt.show()
        
        plt.clf()
        plt.title("Jacobi")
        plt.xlabel("yr")
        plt.ylabel("AU^2 yr^-2")
        plt.plot(self.Time,self.Jacobi,label="Jacobi Integral")
        plt.legend(loc="best")
        plt.show()
        
        plt.clf()
        plt.title("Distance Between Stars")
        plt.xlabel("yr")
        plt.ylabel("AU")
        plt.plot(self.Time,self.Distance,label="r")
        plt.legend(loc="best")
        plt.show()
        
        plt.clf()
        plt.title("Angular Momentum")
        plt.xlabel("yr")
        plt.ylabel("MSun AU^2 yr^-1")
        plt.plot(self.Time,self.Angular,label="Lz")
        plt.legend(loc="best")
        plt.show()

class StaticParticle:
    def __init__(self,Name,Mass,Position,Velocity):
        self.Nam = Name
        self.Mas = Mass
        self.Pos = Position
        self.Vel = Velocity
        self.Tra = [[],[],[]]
    
    def Evolve(self,dt):
        
        # Euler's Method for Time Evolution
        self.Vel = self.Vel + self.Acc * dt
        self.Pos = self.Pos + self.Vel * dt
        
        # Record the Trajectory
        self.Tra[0].append(self.Pos.x)
        self.Tra[1].append(self.Pos.y)
        self.Tra[2].append(self.Pos.z)

class StaticSimulator:

    ''' 
        Every particle has an internal clock and Static Simulator synchronizes them.
        Calculate collective properties

    '''

    def __init__(self, particles, forceField = [lambda t,position : -(G*position.x/math.pow(position.x**2+position.y**2+position.z**2,1.5)),
                                                lambda t,position : -(G*position.y/math.pow(position.x**2+position.y**2+position.z**2,1.5)),
                                                lambda t,position : -(G*position.z/math.pow(position.x**2+position.y**2+position.z**2,1.5))]):
        self.forceField = forceField
        self.Particles = particles

    def SimulateEuler(self,Time,Steps):
        
        t = 0
        dt = Time/Steps
        
        while t < Time: 
           
            # Time Evolution
            # d/dt v = fx(t,x,y,z)
            # d/dt x = v(t)
            
            t += dt

            for Particle in self.Particles:
                Force = Vector(self.forceField[0](t,Particle.Pos),self.forceField[1](t,Particle.Pos),self.forceField[2](t,Particle.Pos))
                Particle.Vel += dt*Force
                Particle.Pos += dt*Particle.Vel
                Particle.Tra[0].append(Particle.Pos.x)
                Particle.Tra[1].append(Particle.Pos.y)
                Particle.Tra[2].append(Particle.Pos.z)
        
    def SimulateRungeKutta(self,Time,Steps):
        t = 0
        dt = Time/Steps
        
        while t < Time: 
           
            # Time Evolution
            # d/dt v = fx(t,x,y,z) RK4 
            # d/dt x = v(t)        EULER
            
            t += dt

            for Particle in self.Particles:
                
                # k1: t,x
                # k2: t+dt/2, x + dt*k1/2
                # k3: t+dt/2, x + dt*k2/2
                # k4: t+dt  , x + dt*k3

                k1  = Vector(self.forceField[0](t,Particle.Pos),self.forceField[1](t,Particle.Pos),self.forceField[2](t,Particle.Pos))
                vk1 = (dt/2)*k1
                xk1 = Particle.Pos + (dt/2)*vk1

                k2  = Vector(self.forceField[0](t+dt/2,xk1),self.forceField[1](t+dt/2,xk1),self.forceField[2](t+dt/2,xk1))
                vk2 = (dt/2)*k2
                xk2 = Particle.Pos + (dt/2)*vk2

                k3  = Vector(self.forceField[0](t+dt/2,xk2),self.forceField[1](t+dt/2,xk2),self.forceField[2](t+dt/2,xk2))
                vk3 = (dt/2)*k3
                xk3 = Particle.Pos + (dt/2)*vk3
                
                k4  = Vector(self.forceField[0](t+dt,xk3),self.forceField[1](t+dt,xk3),self.forceField[2](t+dt,xk3))
                vk4 = (dt)*k4
                xk4 = Particle.Pos + (dt)*vk4

                
                Acceleration  = (k1+k2+k2+k3+k3+k4)
                Particle.Vel += (dt/6.)*Acceleration

                # Fixed to Euler
                Particle.Pos += dt*Particle.Vel
                Particle.Tra[0].append(Particle.Pos.x)
                Particle.Tra[1].append(Particle.Pos.y)
                Particle.Tra[2].append(Particle.Pos.z)

    def GeneratePlots(self):
        plt.clf()
        plt.title("Orbits")
        plt.xlabel("AU")
        plt.ylabel("AU")
        for Particle in self.Particles:
            plt.plot(Particle.Tra[0],Particle.Tra[1],label=Particle.Nam)
        plt.legend(loc="best")
        plt.show()
