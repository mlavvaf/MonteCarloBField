#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 16:06:37 2022

@author: maedeh
"""

import numpy as np

def pos(a): 
    
    x = [-550]
    y = [-440, -200, 0, 200, 440]
    z = [-400, 0, 400]
    
    # a = input('Enter increments in mm:') #increments (mm)
    # a = int(a)
    
    
    for i in x:
        if i < 1100:
            # print(i+a)
            x.append(i + a)
    
    pos = []
    for i in x:
        for j in y:
            for k in z:
                pos.append([i, j, k])
                # print([i, j, k])
                
    
    
    np.savetxt("positions" + str(a) + "mm increments.csv",pos, delimiter =",",  fmt ='% s', header = 'x (mm),y (mm),z (mm)', comments='')
    
    return

# pos(120)