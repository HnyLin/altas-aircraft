# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 17:28:11 2023

@author: airfo
"""

import numpy as np

v = 590                             #ft/s
gamma =  1.4
T = 418.77                          #deg R at 28k
slug_lbm = 32.17
R = 53.35                           #ft*lbf / (lbm R)

AR_w = 10.06133                     #Verify
AR_h = 5.6
nu = 0.97

nu_h = 0.9                          #Martins

lam = 0.25                          #1/4c assumption (Unused)

c_root = 13.59                    #ft (Unused)
c_bar = 9.5                         # ft (mean chord)

sweep_w = 20 * np.pi / 180          #Wing Sweep Radians
sweep_h = 20 * np.pi / 180          #Tail Sweep Radians

body_length = 81                    #ft
wing_tip_loc = 34  
horz_tip_loc = 71              #ft

body_length_frac = lam * c_root / body_length

K_f = 0.115                         #Because 1/4c as a % of length <0.1 used smallest provided Kf from Martins
w_f = 10                            #Fuselage Width ft
L_f = body_length
S_w = 805.063                       #Verify ft^2
S_h = 202.18                        #ft^2
c_h = 9.614                         #ft
l_h = (lam * c_h + horz_tip_loc) - (lam * c_root + wing_tip_loc)          #ft



# MAC 

t_root= 3.968

t = c_root / t_root

MAC = c_root * (2/3) * (( 1 + t + t**2 ) / ( 1 + t ))

mean_ac = 2 / 3 * c_root * (1 + lam + lam ** 2) / (1 + lam)         #Cancels out

#Calculation

M = v / np.sqrt(gamma * R * T * slug_lbm)

dC_lw = 2 * np.pi * AR_w / (2 + np.sqrt( AR_w / nu ) ** 2 * (1 + np.tan(sweep_w) ** 2 - M **2) + 4) 

dC_lh0 = 2 * np.pi * AR_h / (2 + np.sqrt( AR_h / nu ) ** 2 * (1 + np.tan(sweep_h) ** 2 - M **2) + 4) 

d_epsilon  = 2 / (np.pi * AR_w) * dC_lw

dC_lh = dC_lh0 * (1 - d_epsilon) * nu_h

term = K_f * w_f ** 2 * L_f / (S_w * mean_ac * dC_lw)

x_np = (l_h * S_h / ( mean_ac * S_w ) * dC_lh / dC_lw - term) * mean_ac

#SM 

x_np_ac = x_np

x_np = x_np + (lam * c_root + wing_tip_loc) 

x_cg = 489/12

SM = (x_np - x_cg) /mean_ac


 