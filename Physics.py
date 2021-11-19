import math
import matplotlib.pyplot as plt
import numpy as np 
import functools

G = 4*np.pi**2

# Time is measured in years
# Distance is Measured in AU
# Mass is measured in MSun

class Vector:
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
    
    def __add__(self,other):
        return Vector(self.x+other.x,self.y+other.y,self.z+other.z)
    
    def __sub__(self,other):
        return Vector(self.x-other.x,self.y-other.y,self.z-other.z)
    
    @functools.singledispatchmethod
    def __mul__(self, other):
        return (self.x*other.x + self.y*other.y+self.z*other.z)
    
    @__mul__.register(float)
    def _(self,other):
        return Vector(self.x*other,self.y*other,self.z*other)
    
    @functools.singledispatchmethod
    def __rmul__(self, other):
        return (self.x*other.x + self.y*other.y+self.z*other.z)
    
    @__rmul__.register(float)
    def _(self,other):
        return Vector(self.x*other,self.y*other,self.z*other)
    
    def __matmul__(self,other):
        return Vector(self.y*other.z-self.z*other.y,self.z*other.x-self.x*other.z,self.x*other.y-self.y*other.x)
    
    def __str__(self):
        return "("+str(self.x)+","+str(self.y)+","+str(self.z)+")"

class Particle:
    def __init__(self,Name,Mass,Position,Velocity,Acceleration):
        self.Nam = Name
        self.Mas = Mass
        self.Pos = Position
        self.Vel = Velocity
        self.Acc = Acceleration
        self.Tra = [[],[],[]]
    
    def Evolve(self,dt,Omega=Vector(0,0,0)):
        
        self.Vel.x += self.Acc.x * dt
        self.Pos.x += self.Vel.x * dt
        
        self.Vel.y += self.Acc.y * dt
        self.Pos.y += self.Vel.y * dt
        
        self.Vel.z += self.Acc.z * dt
        self.Pos.z += self.Vel.z * dt
        
        self.Tra[0].append(self.Pos.x)
        self.Tra[1].append(self.Pos.y)
        self.Tra[2].append(self.Pos.z)

class Simulator:
    
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
    
    def Simulate(self,Time,Steps):
        
        t = 0
        dt = Time/Steps
        
        while t < Time: 
            t += dt
            k  = 0
            u  = 0
            lz = 0
        
            # Interactions
        
            for i, Particle1 in enumerate(self.Particles):
            
                Particle1.Acc = Vector(0,0,0)
                
                # Coriolis Acceleration
                Particle1.Acc = Particle1.Acc - 2.0*(self.Omega @ Particle1.Vel)
                
                # Centrifugal Acceleration
                Particle1.Acc = Particle1.Acc - self.Omega @ (self.Omega @ Particle1.Pos)
                
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
        plt.plot(self.Time,self.Potential,label="Potential Energy")
        plt.plot(self.Time,self.Kinetic,label="Kinetic Energy")
        plt.legend(loc="best")
        plt.show()
        
        plt.clf()
        plt.title("Angular Momentum")
        plt.xlabel("yr")
        plt.ylabel("MSun AU^2 yr^-1")
        plt.plot(self.Time,self.Angular,label="Lz")
        plt.legend(loc="best")
        plt.show()
