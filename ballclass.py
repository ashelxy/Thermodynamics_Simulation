# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 15:21:19 2019

@author: xl7817
"""

import numpy as np
from numpy import linalg as LA
import pylab as pl

class Ball:

    """
    
    Class to create a ball or container for the simulation
    
    """
    
    ball_counter = 0
    
    def __init__(self, m=1, r=1, pos=[0,0], v=[0,0], taipu="ball"):

        self.__mass = m                   # in kgs
        self.__radius = r                 # in m
        self.__position = np.array(pos)   # each in m
        self.__velocity = np.array(v)     # in ms-1
        self.__type = taipu
        
        if self.__type == 'ball':	 
            Ball.ball_counter += 1    
        

        """

    	The following lines are simply to ensure that the ball can be 
        initialised properly and are consistent for the 
        remainder of the project.
         
        """
        
        if type(m) != int and type(m) != float:
            raise ValueError("Mass described isn't valid!")
        if type(r) != int and type(r) != float:
            raise ValueError("Radius of ball isn't valid!")
        if len(self.__position) != 2:
            raise ValueError("Please insert a 2 var array for your position!")
        if len(self.__velocity) != 2:
            raise ValueError("Please insert a 2 vare array for your velocity!")
        if str(self.__type) != "ball" and str(self.__type) != "container":
            raise ValueError("Please insert a ball or container!")
       
        # Initialisation of the images the ball will created here
        
        if self.__type == "ball":
            self.circle = pl.Circle(self.__position, 
                                    self.__radius, 
                                    fc='cyan', ec='black')
            
        if self.__type == "container":
            self.circle = pl.Circle(self.__position, 
                                    self.__radius, ec='black', 
                                    fill=False, ls='solid')
                                   
    def radius(self):
        
        return self.__radius
    
    def pos(self):
        
        return self.__position
    
    def vel(self):
        
        return  self.__velocity
    
    def speed(self):
        
        """
        
        This returns the magnitude of the velocity, so LA.norm() doesn't 
        have to be used each time. Keeps the code tidier.
        
        """
        
        speed = LA.norm(self.__velocity)
        return speed
    
    def momentum(self):

        momentum = self.__mass * self.speed()
        return momentum
        
    def move(self,dt): #dt is in seconds
        
        self.__position = np.add(self.__position, dt*self.__velocity)
        self.circle.center = self.__position    # updates the position patch
        return self.__position
    
    def time_to_collision(self, other):
        
        r12 = np.subtract(self.__position, other.__position)
        v12 = np.subtract(self.__velocity, other.__velocity)
        
        if self.__type == 'ball' and other.__type == 'ball':
            R12 = self.__radius + other.__radius
        else: 
            R12 = self.__radius - other.__radius

        d = (np.dot(r12, v12)**2) - \
            (np.dot(v12, v12)*(np.dot(r12, r12) - R12**2))

        if d < 0:
            t1 = None
            return t1

        else:
            t1 = -(np.dot(r12, v12) / np.dot(v12, v12)) + \
                 (np.sqrt(d) / np.dot(v12, v12)) - 1e-10
                  
            t2 = -(np.dot(r12, v12) / np.dot(v12, v12)) - \
                 (np.sqrt(d) / np.dot(v12, v12)) - 1e-10
                 
            # The additional (-1e-10) at the back of the times are to 
            # take into account the weird rounding thing Python does.
            
            """
            
            This keeps times tidy by only returning the next time 
            to collision, returns an arbitrary negative number 
            if balls have already collided.
            
            """
            
            if t1 > 0 and t2 > 0:
                return min(t1, t2)
            elif t1 > 0 or t2 > 0:
                return max(t1, t2)
            else:
                return -10
            
            
    def kenergy(self):
        
        """
        
        This function determines the kinetic energy of a ball at any time 
        in the container.
        
        """
        
        KE = abs(self.__mass * (LA.norm(self.__velocity)**2) / 2)
        
        return KE
    
    def collide(self, other):

        """

    	These equations are derived from the conservation of momentum
        and the conservation of kinetic energy.

    	"""
        
        if self.__type == "container" and other.__type == "container":
            raise ValueError ("Two containers can't collide!")
            
        m1 = self.__mass
        m2 = other.__mass
        M  = m1 + m2
        
        r1r2 = self.__position - other.__position     
        r2r1 = other.__position - self.__position
        
        r12 = np.dot(r1r2, r1r2)
        r21 = np.dot(r2r1, r2r1)
        
        upar1 = ((np.dot(self.__velocity, r1r2) / r12) * r1r2)      
        upar2 = ((np.dot(other.__velocity, r2r1) / r21) * r2r1)    
        
        uper1 = self.__velocity - upar1         		           
        uper2 = other.__velocity - upar2               	          
            
        if self.__type == "container":
            vpar1 = upar1
            vpar2 = -upar2
                    
        elif other.__type == "container":
            vpar1 = -upar1
            vpar2 = upar2
       
        else:
            vpar1 = (upar1 * (m1 - m2) + upar2  *(2*m2)) / M
            vpar2 = (upar1 * (2*m1) + upar2 * (m2 - m1)) / M
                
        self.__velocity = np.add(uper1, vpar1)
        other.__velocity = np.add(uper2, vpar2)
                                  
        return self.__velocity, other.__velocity
            
                    
    def get_patch(self):
        
        return self.circle
            
    def ball_count(self):
        
        return self.ball_counter
    

        
        
        
    
    
