# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 12:23:24 2019

@author: xl7817
"""


import ballclass as bc
import numpy as np
from numpy import linalg as LA
import scipy.constants as const
import scipy.stats as stats
import pylab as pl
import random as rnd
import operator as op
import matplotlib.pyplot as plt



class Simulation(bc.Ball):

    """
    
    A class for the simulation of balls and a container.
    
    """

    def __init__(self, num):  
        
        self._cont = bc.Ball(m=1e15, r=10, taipu = "container")
        self._ball = [self._cont]
        self._dt = 0    # This is a variable added for the calculation of
                        # pressure acting on the container walls.
        
        if num < 1:
            raise ValueError("Please insert a number >= 1!")
        
        k = 500  # change this variable to easily change the range
                 # of velocities of the balls.
        
        for x in range(num):
            self._ball.append( bc.Ball( m = (32 * 1.66e-27), r = 0.1,
                                    
                               pos = [6.5 - (x%27 * 0.5), 
                                      6.5 - (int(x/27) * 0.5)],   
#                               pos = [0.8 - (x%4 * 0.5),
#                                      0.8 - int(x/4) * 0.5],
                               
                               # the two position vectors are the ones used to 
                               # plot the values for the graphs.
                               
                               v = [rnd.uniform(-k,k), 
                                    rnd.uniform(-k,k)]) )
           
                    # This initialises the balls to be in a 5x5 grid.
                    # The balls are of equal densities with radius and mass
                    # ranging between 0.5 and 1.5.
                    # With the uniform probability, average velocity is 0.
        
    def next_collision(self):
        
        """

        This code works by creating a list of positive times and selecting 
        the next two balls that collide. The time is then selected through 
        sorting the list and all balls are moved to the next collision.

    	"""
        
        timelist = []
        
        for x in range (len(self._ball)):
            for y in range (x+1, len(self._ball)):
                t = self._ball[x].time_to_collision(self._ball[y])
                if t is not None and t > 1e-9:
                    timelist.append([t, x, y])
                    
        mint = op.itemgetter(0)
        
        nxt = sorted(timelist, key = mint)
        
        self._dt = nxt[0][0]       # time
        a = nxt[0][1]       # ball 1
        b = nxt[0][2]       # ball 2
        
        for k in range (len(self._ball)):
            self._ball[k].move(self._dt)   
        
        self._ball[a].collide(self._ball[b])
       

    def kesys(self):
        
        """
        
        This code returns the total kinetic energy in each frame.
        It is included in run() so it prints the KE after every collision. 
        
        """
        
        ke = 0
        
        for x in range (1, len(self._ball)):
            ke += self._ball[x].kenergy()
        
        return ke
       
    def run(self, num_frames, animate=False):
        
        v = [0,0]
        fr = []           
        ke = []
        
        if animate:
            f = pl.figure()
            ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10), aspect='equal')
            
            ax.add_artist(self._cont.get_patch())
            
            for z in range (len(self._ball)):
                ax.add_patch(self._ball[z].get_patch())
                v = np.add(v, self._ball[z].vel())
                
            """
            
            This part moves the code by the time to next collision per frame,
            so the output is not always smooth as some collisions occur 
            almost instantly after each other, while some takes longer.
            
            This is ideal for large numbers of balls as the collision frequency
            is much higher than that with little balls.
            
            """
                
            for frame in range(num_frames):
                fr.append(frame+1)
                ke.append(self.kesys())
                self.next_collision()                
                if animate:
                    pl.pause(0.01)
        
        if animate:
            pl.show()
        
        plt.title("KE of system in each frame (KE check)")
        plt.xlabel("Frame number")
        plt.ylabel("Kinetic energy")
        plt.plot(fr, ke)
        plt.show()
        
                    
    def fakerun(self, num_frames):
        
        """
        
        This allows the system to run without animation as it causes blue 
        screens when animating 729 balls. 
            
        """
        
        fr = []           
        ke = []
        
        for frame in range(num_frames):
            fr.append(frame+1)
            ke.append(self.kesys())
            self.next_collision()
        
        # KE check, should just be a flat line        
        
        plt.title("KE of system in each frame (KE check)")
        plt.xlabel("Frame number")
        plt.ylabel("Kinetic energy")
        plt.plot(fr, ke)
        plt.show()
        

    def hist(self):
               
        fromcenter = []     # list of distances from the center
        
        for z in range (1, len(self._ball)):
            A = self._ball[z].pos()
            fromcenter.append(LA.norm(A))
            
        dballs = []         # list of distances between each ball in simulation
        
        for x in range (1, len(self._ball)):
            for y in range (x+1, len(self._ball)):
                k = self._ball[x].pos() - self._ball[y].pos()
                dballs.append(LA.norm(k))
                
        plt.subplot(1,2,1)
        plt.title("Distance of each ball from center")
        plt.hist(fromcenter, bins = np.linspace(0,10,51), ec = 'black')
        print("Mean distance from center:", np.mean(fromcenter))
        
        plt.subplot(1,2,2)
        plt.title("Distance between each ball")
        plt.hist(dballs, bins = np.linspace(0,20,101), ec = 'black')
        print("Mean distance from each ball:", np.mean(dballs))

        # The more balls there are, the more plots there are in the second
        # histogram as there are n! values.        
        
    def boltzmann(self):
        vels = []
        for x in range (1, len(self._ball)):
            vels.append(self._ball[x].speed())
        
        plt.figure()
        x = np.linspace(0,1000,1000)
        plt.title("Velocity of balls against Boltzmann distribution")         
        plt.hist(vels, bins = 200, density = 1, ec = 'black')     
        params = stats.maxwell.fit(vels, floc = 0)     
        plt.plot(x, stats.maxwell.pdf(x, *params), lw=3)     
        plt.xlabel("Velocity")     
        plt.ylabel("Proportion of Particles")   
        plt.grid = True      
        plt.show()  
        
        
    def frameke(self):
        ke = []
        
        for x in range (1, len(self._ball)):
            ke.append(self._ball[x].kenergy())
        
        plt.figure()
        plt.title("KE of each ball in final frame")
        plt.hist(ke, bins = 200, ec = 'black')
        print("Mean KE of each ball:", np.mean(ke))
    
    def temp(self):
        
        """
        
        The equation used in this function is E = NKbT,
        where E is the sum of KE in the system.
        T = E/NKb
        
        """
        
        E = self.kesys()
        N = len(self._ball)
        
        T = E/(N*const.k)
        print('Temperature is:', T, 'K')
        
        return T
    
    def pressure(self):
        
        # As this is 2D, area is the circumference of the container
        area = 2*np.pi*self._ball[0].radius()
        mom = []
        
        for x in range (1, len(self._ball)):
            p = self._ball[x].momentum()
            mom.append(p)
        
        totmom = np.sum(mom)
        
        F = totmom/self._dt
        P = F / area
        print('Pressure is:', P, 'F/m')
        
        return P
    
                
                
                
                
                
                                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        