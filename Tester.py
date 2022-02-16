#!/usr/bin/env python3

import numpy as np
from Physics import *
from Tools import *

def DimFrame(Frame,Factor):
    Dimensions = Frame.shape
    for y in range(Dimensions[0]):
        for x in range(Dimensions[1]):
            for rbg in range(3):
                Frame[y][x][rbg] = int((Frame[y][x][rbg]*Factor)//1)

def main():

    Particle1 = Particle("Sun 1",3.0,Vector(1,0,0),Vector(0,2*np.arccos(-1.),0),RED)
    Particle2 = Particle("Sun 2",1.0,Vector(-3,0,0),Vector(0,0,0),BLUE)
    
    Particles = [Particle1,Particle2]
    for i in range(1):
        Position = Vector(1.5,0,0)
        Velocity = Vector(0,-10,0)
        Mass = 1E-3 # Irrelevant unless interactions are turned-on
        Particles.append(Particle("Planet "+str(i+1),Mass,Position,Velocity,WHITE))
        
    S = Simulator(Particles)
    S.GoToCM()
    S.Simulate(1.75,1000)
    S.GenerateVideo("Binary System")
    S.GeneratePlots(mode="Save")


if __name__ == '__main__':
    main()