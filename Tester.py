#!/usr/bin/env python3

import numpy as np
from Physics import *
from Tools import *

def main():    

    Sun1 = Particle("Sun 1",3.0,Vector(1,0,0), Vector(0,2*np.arccos(-1.),0),RED)
    Sun2 = Particle("Sun 2",1.0,Vector(-3,0,0),Vector(0,0,0),BLUE)

    BinaryStarPotential  = ForceFieldGenerator([Sun1,Sun2])
    
    Particles = []
    for i in range(2):
        Position = Vector(Random(1,2),0,0)
        Velocity = Vector(0,Random(-12,-10),0)
        Mass = 1E-6   # Irrelevant unless interactions are turned-on
        Particles.append(Particle("Planet "+str(i+1),Mass,Position,Velocity,WHITE))
    
    Test_Static = True
    Time  = 0.2
    Steps = 750

    if Test_Static == True:

        BinaryStarSimulation = StaticSimulator(Particles,BinaryStarPotential)
        BinaryStarSimulation.GoToCM()
        BinaryStarSimulation.Simulate(Time,Steps,"RK4")
        BinaryStarSimulation.GeneratePlots(mode="Save")
        BinaryStarSimulation.GenerateVideo("Binary System 2")
    
    else:

        BasicSimulator = Simulator(Particles+[Sun1,Sun2])
        BasicSimulator.GoToCM()
        BasicSimulator.Simulate(Time,Steps)
        BasicSimulator.GeneratePlots(mode="Save")
        BasicSimulator.GenerateVideo("Binary System 1")

if __name__ == '__main__':
    main()