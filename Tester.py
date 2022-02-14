#!/usr/bin/env python3

from Physics import *

def main():

    Particle1 = StaticParticle("Earth",0.0001,Vector(2.1,0,0),Vector(0,1,0))

    SS = StaticSimulator([Particle1])
    # SS.SimulateEuler(20,10000)
    SS.SimulateRungeKutta(20,10000)
    SS.GeneratePlots()

if __name__ == '__main__':
    main()