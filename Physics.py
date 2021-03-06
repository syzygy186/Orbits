
import math, os
import matplotlib.pyplot as plt
import numpy as np 
import functools
from Video import *

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

    def GetTheta(self):
        return np.arctan2(self.y,self.x)

class Particle:
    
    '''
    Particles for the Simulator
    They hold kinematical variables: Position, Velocity and Acceleration
    And the trajectory of the particle, the time should be managed externally.
    
    '''
    
    def __init__(self,Name,Mass,Position,Velocity,Color):
        self.Nam = Name
        self.Mas = Mass
        self.Pos = Position
        self.Vel = Velocity
        self.Acc = Vector(0.,0.,0.)
        self.Col = Color
        self.Tra = [[],[],[]]

        self.IniPos = self.Pos
        self.IniVel = self.Vel

    def RemoveLastTrack(self):
        
        for i in range(3):
            del self.Tra[i][-1]

    def ClearTrajectory(self):
        self.Tra = [[],[],[]]
        self.Pos = self.IniPos
        self.Vel = self.IniVel
        
    def Evolve(self,dt,Omega=Vector(0.,0.,0.)):
        
        self.Vel = self.Vel + self.Acc * dt
        self.Pos = self.Pos + self.Vel * dt
        
        # Record the Trajectory
        self.Tra[0].append(self.Pos.x)
        self.Tra[1].append(self.Pos.y)
        self.Tra[2].append(self.Pos.z)
  
##################################

class Simulator:
    
    '''
    
    Deals with Time management and Interactions
    Democratic treatment of all particles

    
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
        self.Theta     = []
    
    # Go to the Center of mass
    def GoToCM(self):

        print("Going to the CM")
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
            print(Particle.Nam,Particle.Pos,Particle.Vel)
    
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

    def GeneratePlots(self,mode="Show"):
        
        plt.clf()
        plt.title("Orbits")
        plt.xlabel("AU")
        plt.ylabel("AU")
        for Particle in self.Particles:
            plt.plot(Particle.Tra[0],Particle.Tra[1],label=Particle.Nam)
        plt.legend(loc="best")
        if(mode=="Show"):
            plt.show()
        elif (mode=="Save"):
            plt.savefig("Orbits.png")
        
        plt.clf()
        plt.title("Energies")
        plt.xlabel("yr")
        plt.ylabel("MSun AU^2 yr^-2")
        plt.plot(self.Time,self.Energy,label="Total Energy")
        plt.plot(self.Time,self.Potential,label="Potential Energy")
        plt.plot(self.Time,self.Kinetic,label="Kinetic Energy")
        plt.legend(loc="best")
        if(mode=="Show"):
            plt.show()
        elif (mode=="Save"):
            plt.savefig("Energies.png")
        
        plt.clf()
        plt.title("Jacobi")
        plt.xlabel("yr")
        plt.ylabel("AU^2 yr^-2")
        plt.plot(self.Time,self.Jacobi,label="Jacobi Integral")
        plt.legend(loc="best")
        if(mode=="Show"):
            plt.show()
        elif (mode=="Save"):
            plt.savefig("Jacobi.png")
        
        plt.clf()
        plt.title("Distance Between Stars")
        plt.xlabel("yr")
        plt.ylabel("AU")
        plt.plot(self.Time,self.Distance,label="r")
        plt.legend(loc="best")
        if(mode=="Show"):
            plt.show()
        elif (mode=="Save"):
            plt.savefig("Distace between Satrs.png")
        
        plt.clf()
        plt.title("Angular Momentum")
        plt.xlabel("yr")
        plt.ylabel("MSun AU^2 yr^-1")
        plt.plot(self.Time,self.Angular,label="Lz")
        plt.legend(loc="best")
        if(mode=="Show"):
            plt.show()
        elif (mode=="Save"):
            plt.savefig("Total Angular Momentum.png")

    def GenerateVideo(self,Name):
        VG = VideoGenerator(Name,self.Particles)
        VG.Generate()

    def ClearTrajectories(self):
        for Particle in self.Particles:
            Particle.ClearTrajectory()

##################################

# Differentiates between two classes of particles 

class StaticSimulator:

    ''' 
        
        There are two kinds of particles: Sources and Particles

            Sources are used to generate force fields 
            Particles eveolve in the generated potential interction-free

    '''

    def __init__(self, particles, forceField):
        self.Methods = {"Euler":self.SimulateEuler,"RK4":self.SimulateRK4}
        self.forceField = forceField
        self.Particles  = particles

    def GoToCM(self):

        '''
            Shifts into the Center of mass of the Sources and then 
            shifts all particles into this frame
        '''
        # print("Going into the CM of the Force Field")
        self.forceField.GoToCM()
        self.Position = self.forceField.Position
        self.Velocity = self.forceField.Velocity
        for Particle in self.Particles:
            Particle.Pos = Particle.Pos - self.Position
            Particle.Vel = Particle.Vel - self.Velocity
            # print(Particle.Nam,Particle.Pos,Particle.Vel)

    def Simulate(self,Time,Steps,Method):

        '''
            Top Level Simulator Selector
        '''

        self.Methods[Method](Time,Steps)

    def SimulateEuler(self,Time,Steps):

        '''
            Simualtion using Euler's Method for the evolution of the 
            Particles in the force field

        '''
        
        t = 0
        dt = Time/Steps
        
        while t < Time: 
           
            t += dt
            self.forceField.Evolve(dt,"Euler")

            for Particle in self.Particles:
                Particle.Acc = self.forceField(Particle.Pos)
                Particle.Evolve(dt)

    def SimulateRK4(self,Time,Steps):

        '''
            Simualtion using Runge-Kutta of fourth degree for the evolution 
            of the particles in the force field
            Arguments:

                Time  : Time in Years to be Simulated
                Steps : Number of Steps to be used to reach the endpoint time

        '''

        t = 0
        dt = Time/Steps
        
        while t < Time: 
           
            t += dt

            for Particle in self.Particles:
                
                k1  = self.forceField(Particle.Pos)
                vk1 = (dt/2)*k1
                xk1 = Particle.Pos + (dt/2)*vk1

                # Here we step on the Force Field but delete the track 
                # to keep the trajectories of sources and particles of equal length

                self.forceField.Evolve(dt/2,"RK4")
                self.forceField.RemoveLastTrack()

                k2  = self.forceField(xk1)
                vk2 = (dt/2)*k2
                xk2 = Particle.Pos + (dt/2)*vk2

                k3  = self.forceField(xk2)
                vk3 = (dt/2)*k3
                xk3 = Particle.Pos + (dt/2)*vk3
                
                # Second step on the Force Field, now we keep the track

                self.forceField.Evolve(dt/2,"RK4")

                k4  = self.forceField(xk3)
                vk4 = (dt)*k4
                xk4 = Particle.Pos + (dt)*vk4
                
                Acceleration  = (1./6.)*(k1+(2.*k2)+(2.*k3)+k4)
                Particle.Acc  = Acceleration

                Particle.Evolve(dt,"RK4")

    def GeneratePlots(self,args):
        plt.clf()
        plt.title("Orbits")
        plt.xlabel("AU")
        plt.ylabel("AU")
        for Particle in self.Particles:
            plt.plot(Particle.Tra[0],Particle.Tra[1],label=Particle.Nam)
        plt.legend(loc="best")
        if(args["mode"]=="Show"):
            plt.show()
        elif (args["mode"]=="Save"):
            try:
                filename = args["filename"]
            except:
                print("Error: No filename given!")
                return
            
            plt.savefig(filename+".png")

    def GenerateVideo(self,Name):
        VG = VideoGenerator(Name,self.Particles)
        VG.Generate()

    def ClearTrajectories(self):
        self.forceField.ClearTrajectories()
        for Particle in self.Particles:
            Particle.ClearTrajectory()

# Any cofiguration of stars that we wish to study
# they all interract with themselves but not with 
# the testers from the Satatic Simulator

class ForceFieldGenerator:

    '''
        Holds a set of particles to be used to generate a 
        potential that can be used by a Simulator as a 
        Force Field Source

            Holds:
                Particles : All sources which evolve through Gravitational Forces 

    '''

    def __init__(self,Particles,force = lambda r : (G/pow(r*r,1.5))*r):
        self.Methods = {"Euler":self.EvolveEuler, "RK4":self.EvolveRK4}
        self.Particles = Particles
        self.Force     = force
        
    def GravitationalForce(self,r):

        '''
           Wrwapper to be able to modify the underlying force field

        '''

        return self.Force(r)
    
    def RemoveLastTrack(self):

        '''
            Recursively remove last track from all undelying particles
        '''

        for Particle in self.Particles:
            Particle.RemoveLastTrack()

    def GoToCM(self):

        '''
            Shift to the Center of Mass System
            Here rCM = 0 and (d/dt)rCM = 0

            We store the transformation values to be able to go back

        '''

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
            # print(Particle.Nam,Particle.Pos,Particle.Vel)

    def GeneratePlots(self,mode="Save"):
        plt.clf()
        plt.title("Force Field Orbits")
        plt.xlabel("AU")
        plt.ylabel("AU")
        for Particle in self.Particles:
            plt.plot(Particle.Tra[0],Particle.Tra[1],label=Particle.Nam)
        plt.legend(loc="best")
        if(mode=="Show"):
            plt.show()
        elif (mode=="Save"):
            plt.savefig("Force Field Orbits.png")

    def Evolve(self,dt,Method="RK4"):
        self.Methods[Method](dt)
        
    def EvolveRK4(self,dt):
            for i,Particle1 in enumerate(self.Particles):
                Particle1.Acc = Vector(0,0,0)
                for j,Particle2 in enumerate(self.Particles):
                    if i == j:
                        continue
                    r  = Particle2.Pos - Particle1.Pos
    
                    k1  = self.GravitationalForce(r)
                    vk1 = (dt/2)*k1
                    xk1 = r + (dt/2)*vk1
    
                    k2  = self.GravitationalForce(xk1)
                    vk2 = (dt/2)*k2
                    xk2 = r + (dt/2)*vk2
    
                    k3  = self.GravitationalForce(xk2)
                    vk3 = (dt/2)*k3
                    xk3 = r + (dt/2)*vk3
                    
                    k4  = self.GravitationalForce(xk3)
                    vk4 = (dt)*k4
                    xk4 = r + (dt)*vk4
                    
                    Acceleration  = (1./6.)*(k1+(2.*k2)+(2.*k3)+k4)
                    Particle1.Acc  = Particle1.Acc + Acceleration 

            for Particle in self.Particles:
                Particle.Evolve(dt)

    def EvolveEuler(self,dt):
        for i,Particle1 in enumerate(self.Particles):
            Acceleration = Vector(0,0,0)
            for j,Particle2 in enumerate(self.Particles):
                if i == j:
                    continue
                r  = Particle2.Pos - Particle1.Pos
                Acceleration = Acceleration + Particle2.Mas*self.GravitationalForce(r)
            Particle1.Acc = Acceleration
        for Particle in self.Particles:
            Particle.Evolve(dt)
        
    def __call__(self,Position):
        TotalForce = Vector(0,0,0)
        for i,Particle1 in enumerate(self.Particles):
            r = Particle1.Pos - Position
            TotalForce = TotalForce + Particle1.Mas*self.GravitationalForce(r)
        return TotalForce

    def FindRadius(self):
        Radius = 0
        for i,Particle1 in enumerate(self.Particles):
            for j,Particle2 in enumerate(self.Particles):
                if i==j:
                    continue
                r = Particle1.Pos - Particle2.Pos
                r = math.sqrt(r*r)
                if (r > Radius):
                    Radius = r
        return r

    def ClearTrajectories(self):
        for Particle in self.Particles:
            Particle.ClearTrajectory()

##################################

# Currently Unused
class Tracker:
   
    '''
        Kimenatical Tracker

    '''

    def __init__(self,Title,XTitle,YTitle,XUnits,YUnits):
        
        self.Title    = Name
        self.Lambda  = Lambda
        
        self.XTitle  = XTitle
        self.XUnits  = XUnits
        self.XAxis   = []
        
        self.YTitle  = YTitle
        self.YUnits  = YUnits
        self.YAxis   = []
    
    def AddPoint(self,X,Y):
        XAxis.append(X)
        YAxis.append(Y)

    def GeneratePlot(self,Mode="Save"):
        plt.clf()
        plt.title(self.Name)
        plt.xlabel(self.XTitle+" ("+self.XUnits+")")
        plt.xlabel(self.YTitle+" ("+self.YUnits+")")
        for Particle in self.Particles:
            plt.plot(Particle.Tra[0],Particle.Tra[1],label=Particle.Nam)
        plt.legend(loc="best")
        if(Mode=="Show"):
            plt.show()
        elif (Mode=="Save"):
            plt.savefig(self.Name+".png")
        else:
            print("Error: Unrecognized tracker plot generator mode:",Mode)

# 
class Simulator2:

    def __init__(self,Name,Particles=[]):
        self.Name = Name
        self.Velocity  = Vector(0,0,0)
        self.Position  = Vector(0,0,0)
        self.Omega     = Vector(0,0,0)
        self.DOmega    = Vector(0,0,0)
        self.Velocity  = Vector(0,0,0)
        self.Particles = Particles
    
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
        

        GravitationalField = lambda r1,r2 : (-G/math.pow((r1-r2)*(r1-r2),1.5))*(r1-r2)
        NonInertialField   = lambda r,v   : -2.0*(self.Omega @ v) - self.Omega @ (self.Omega @ r) - self.DOmega @ r 

        while t < Time: 

            t += dt
        
            for i, Particle1 in enumerate(self.Particles):
            
                r1 = Particle1.Pos
                v1 = Particle1.Vel 
                Particle1.Acc = Vector(0,0,0)

                Particle1.Acc = Particle1.Acc + NonInertialField(r1,v1)

                # k1  = NonInertialField(r1,v1)
                # vk1 = (dt/2)*k1
                # xk1 = r1 + (dt/2)*vk1
    
                # k2  = NonInertialField(xk1,v1+vk1)
                # vk2 = (dt/2)*k2
                # xk2 = r1 + (dt/2)*vk2
    
                # k3  = NonInertialField(xk2,v1)
                # vk3 = (dt/2)*k3
                # xk3 = r1 + (dt/2)*vk3
                    
                # k4  = NonInertialField(xk3,v1)

                # Acceleration  = (1./6.)*(k1+(2.*k2)+(2.*k3)+k4)

                # Particle1.Acc = Particle1.Acc + Acceleration
                
                # Coriolis Acceleration
                # Particle1.Acc = Particle1.Acc - 2.0*(self.Omega @ Particle1.Vel)
                
                # Centrifugal Acceleration
                # Particle1.Acc = Particle1.Acc - self.Omega @ (self.Omega @ Particle1.Pos)

                #Euler Force
                # Particle1.Acc = Particle1.Acc - self.DOmega @ Particle1.Pos

                for j,Particle2 in enumerate(self.Particles):
                    if i == j:
                        continue

                    r2 = Particle2.Pos
                    
                    # k1  = GravitationalField(r1,r2)
                    # vk1 = (dt/2)*k1
                    # xk1 = r1 + (dt/2)*vk1
    
                    # k2  = GravitationalField(xk1,r2)
                    # vk2 = (dt/2)*k2
                    # xk2 = r1 + (dt/2)*vk2
    
                    # k3  = GravitationalField(xk2,r2)
                    # vk3 = (dt/2)*k3
                    # xk3 = r1 + (dt/2)*vk3
                    
                    # k4  = GravitationalField(xk3,r2)
                    
                    # Acceleration  = (1./6.)*(k1+(2.*k2)+(2.*k3)+k4)
                    
                    # Particle1.Acc = Particle1.Acc + Particle2.Mas*Acceleration

                    Particle1.Acc = Particle1.Acc + Particle2.Mas*GravitationalField(r1,r2)

           
            # Time Evolution
        
            for Particle in self.Particles:
                Particle.Evolve(dt,self.Omega)

    def GeneratePlots(self,mode="Save"):
        plt.clf()
        plt.title(self.Name+" Orbits")
        plt.xlabel("AU")
        plt.ylabel("AU")
        for Particle in self.Particles:
            plt.plot(Particle.Tra[0],Particle.Tra[1],label=Particle.Nam)
        plt.legend(loc="best")
        if(mode=="Show"):
            plt.show()
        elif (mode=="Save"):
            plt.savefig(self.Name+" Orbits.png")

        plt.clf()
        plt.title(self.Name+" Distances")
        plt.xlabel("dt")
        plt.ylabel("AU")
        xAxis = [ i for i in range(len(self.Particles[0].Tra[0]))]
        for i,Particle1 in enumerate(self.Particles):
            for j,Particle2 in enumerate(self.Particles):
                if i>=j:
                    continue
                yAxis = []
                for t in range(len(xAxis)):
                    r = Vector(Particle1.Tra[0][t],Particle1.Tra[1][t],Particle1.Tra[2][t])
                    r = r - Vector(Particle2.Tra[0][t],Particle2.Tra[1][t],Particle2.Tra[2][t])
                    yAxis.append(math.sqrt(r*r))
                plt.plot(xAxis,yAxis,label=Particle1.Nam+"-"+Particle2.Nam+" Distance")
        plt.legend(loc="best")
        if(mode=="Show"):
            plt.show()
        elif (mode=="Save"):
            plt.savefig(self.Name+" Distances.png")

    def GenerateVideo(self):
        VG = VideoGenerator(self.Name+" Video",self.Particles)
        VG.Generate()

##################################

class OrbitFinder(StaticSimulator):

    def __init__(self,Particles,ForceFieldGenerator):
        super().__init__(Particles,ForceFieldGenerator)
        self.Scale  = 1.0
        self.Factor = 10.

    def FindScale(self,Time,dt):
        Distance = 0.0
        t = 0

        while(t < Time):
            self.forceField.Evolve(dt,"Euler")
            R = self.forceField.FindRadius()
            if (R > Distance):
                Distance = R
            t += dt

        self.Scale = Distance

