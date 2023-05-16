# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 12:12:29 2023

@author: henry
"""
import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLines

#A3 Code
def tw_wp( WS, nu_p, tw, C_L, rho ):
    lbf_lbm = 32.17
    V_cruise = np.sqrt ( 2 / ( rho * C_L ) * WS * lbf_lbm)
    wp = nu_p / V_cruise * (1 / tw)
    return wp

def get_PS_Plot(weight, precision, lbf_lbm, hp, WS, nu_p, AR, w_to, w_L, BFL, rho_SL, rho, 
                  C_lmax_TO, k_1, k_2, C_d0_TO, mu_g, S_FL, S_L, S_a, C_lmax_L, 
                  rho_cruise, V_cruise, C_D0_cruise, e_cruise, G_ceiling, e, 
                  rho_ceiling, C_Lmax, G_climbvals, k_svals, C_D0vals, 
                  C_Lmax_vals, labels, e_vals, w_vals, N_correction):
    
    plt.figure(figsize = (16,10) )
    plt.xlabel("S (ft^2)")
    plt.ylabel("P (hp)")

    #Takeoff
    #Sea Level 30ft Obsticle 
    TOP = BFL / 37.5                            #Check Units
    tw_TO = WS / ( (rho_SL / rho_SL) * C_lmax_TO * TOP)
    
    wp_TOSL30 = hp * tw_wp( WS, nu_p, tw = tw_TO, C_L = C_lmax_TO, rho = rho_SL)
    
    plt.plot(weight / WS, weight / wp_TOSL30, label = 'Takeoff Field Length SL 30ft Obsticle')
    
    #Sea Level 50ft Obsticle(SL+ISA + 18 deg day)
    
    tw_TO = 1 / ( k_2 / ( ( k_1 * WS / ( BFL * rho_SL) + 0.72 * C_d0_TO ) / C_lmax_TO + mu_g ))
    
    wp_TOSL50 = hp * tw_wp( WS, nu_p, tw = tw_TO, C_L = C_lmax_TO, rho = rho_SL)
    plt.plot(weight / WS, weight / wp_TOSL50, label = 'Takeoff Field Length SL+ISA 50ft Obsticle')
    
    #5000 + Sea Level 50ft Obsticle(5000+ISA + 18 deg day)
    tw_TO = 1 / ( k_2 / ( ( k_1 * WS / ( BFL * rho) + 0.72 * C_d0_TO ) / C_lmax_TO + mu_g ))
    
    wp_TO5K50 = hp * tw_wp( WS, nu_p, tw = tw_TO, C_L = C_lmax_TO, rho = rho_SL)
    plt.plot(weight / WS, weight / wp_TO5K50, label = 'Takeoff Field Length 5000+ISA 50ft Obsticle')
    
    #Landing
    #Sea Level 30ft Obsticle
    ws_L = ( ( rho_SL / rho_SL) * C_lmax_L ) / 80 * (S_L - S_a)
    ws_L = ws_L / (w_L / w_to)
    plt.plot(weight / ws_L * np.ones(precision), np.linspace(0, 50000, precision), label = "Landing Field Length")
    
    #Cruise
    q = rho_cruise * V_cruise ** 2 / 2
    
    tw_cruise = q * C_D0_cruise / WS + WS / ( q * np.pi * AR * e_cruise)
    
    wp_cruise = nu_p / V_cruise * 1 / tw_cruise * lbf_lbm * hp
    
    plt.plot( weight / WS, weight / wp_cruise, label = "Cruise" )
    
    #Stall
    wing_area = 805.1
    V_takeoff = np.sqrt( 2 * w_to * lbf_lbm / (C_lmax_TO * wing_area * rho_SL) )
    V_stall = 1/ 1.1 * V_takeoff
    print(V_stall)
    ws_stall = 1 / 2 * rho * V_stall ** 2 * C_Lmax
    
    plt.plot(weight / ws_stall * np.ones(precision), np.linspace(0, 10000, precision), label = "Stall Requirement")

    #FAR 25.121 Climb Requirements
    temp_vals = []
    for i in range(len(G_climbvals)):
        
        tw_C_uncorrected = k_svals[i] ** 2 / C_Lmax_vals[i] * C_D0vals[i] + C_Lmax_vals[i] / ( k_svals[i] ** 2 * np.pi * e_vals[i] * AR ) + G_climbvals[i] 
        
        tw_C = ( 1 / 0.8 ) * ( 1 / 0.94 ) * N_correction[i] * w_vals[i] / w_to * tw_C_uncorrected
            
        wp_C = hp * tw_wp( WS, nu_p, tw = tw_C, C_L = C_Lmax_vals[i], rho = rho_SL )
        
        temp_vals = np.append(temp_vals, wp_C)
        
        plt.plot(weight / WS, weight / wp_C, label = labels[i])
    
    wp_takeoffclimb = temp_vals[0 : precision]
    wp_transegmentclimb = temp_vals[precision : 2 * precision]
    wp_secsegmentclimb = temp_vals[2 * precision : 3 * precision]
    wp_enrouteclimb = temp_vals[3 * precision : 4 * precision]
    wp_balkedAEO = temp_vals[4 * precision : 5 * precision]
    wp_balkedOEI = temp_vals[5 * precision : ]

    #Ceiling
    
    tw_celing = (1 / (rho_ceiling / rho_SL) ** 0.6) * (2 * np.sqrt(C_D0_cruise / (np.pi * AR * e)) + G_ceiling)
    
    wp_celing = hp * tw_wp( WS, nu_p, tw = tw_celing, C_L = C_Lmax, rho = rho_ceiling )
    
    plt.plot(weight / WS, weight / wp_celing, label = "Ceiling")
    
    plt.title("P (hp) - S ($ft^2$)")
    
    #plt.legend(loc = 'best')
    plt.ylim(0, 40000)
    plt.xlim(400, 1600)
    
    labelLines(plt.gca().get_lines(), zorder=2.5)
    
    plt.fill_between( weight/WS,  weight / wp_balkedOEI, y2 = 40000, where = (weight / WS >= weight / ws_L), interpolate=True, color = "green", alpha = 0.2)
    
    plt.scatter(800, 14750, marker = "*", color = "red", s = 1000)
    
    plt.show()
    
    return wp_TOSL30, wp_TOSL50, wp_TO5K50, wp_cruise, wp_takeoffclimb, \
        wp_transegmentclimb, wp_secsegmentclimb, wp_enrouteclimb, \
            wp_balkedAEO, wp_balkedOEI, wp_celing, ws_L

#Inputs
lbf_lbm = 32.17
hp = 550                #550 ft*lbf / s

precision = 900
WS = np.linspace(1, 300, precision)


nu_p = 0.8                                  #From Report
AR = 10.06                                  #From Design
w_to = 82561.1                              #lbs from weight estimate
w_L = 73118                                 #lbs from weight estimate

BFL = 4500                                  #From RFP

#Inputs (Takeoff)
rho_SL = 0.0765                             #lbm/ft^3
rho = 0.0659                                #Air density at 5000'
C_lmax_TO = 1.7                             #Roskam Vol 1 Table 3.1
BPR = 50                                    #Verify BPR for Turbofan
k_2 = 0.75 * (5 + BPR) / (4 + BPR)
k_1 = 0.0447
C_d0_TO = 0.045                             #From Drag Polars
mu_g = 0.025                                #From Roskam

#Inputs (Landing)
S_FL = 4500
S_L = S_FL                                  #Additional Safety Factor deleted (times 0.6)
S_a = 1000                                  #Metabook eq 4.46 (ft)
C_lmax_L = 3.3                              #Roskam Vol 1 Table 3.1

#Inputs(Cruise)
rho_cruise = 10.66*10**-4  * lbf_lbm        #Air density at 28000'
V_cruise = 350 * 1.688                      #Declared in the beginning for W/P
C_D0_cruise = 0.025                         #Drag Polar Estimate
e_cruise = 0.75                             #Drag Polar Estimate

#Inputs (Ceiling)
G_ceiling = 0.001
e = 0.8
rho_ceiling = 0.0287                        #air density at 30,000   

#FAR 25.121 Climb Inputs
#Takeoff Climb
#Transition Segment Climb
#Second Segment Climb
#Enroute Climb
#Balked Landing Climb
C_Lmax = 1.7
G_climbvals = [0.012, 0, 0.024, 0.012, 0.032, 0.021]
k_svals = [1.2, 1.15, 1.2, 1.25, 1.3, 1.5]
C_D0vals = [0.045 + 0.025, 0.045 + 0.025, 0.045, 0.025, 0.125 + 0.1, (0.125 + 0.1) / 2]
C_Lmax_vals = [C_lmax_TO, C_lmax_TO, C_lmax_TO, C_Lmax, C_lmax_L, 0.85 * C_lmax_L]
labels = ["Takeoff Climb", "Transition Segment Climb", "Second Segment Climb", 
          "Enroute Climb", "Balked Landing Climb AEO", "Balked Landing Climb OEI"]
e_vals = [0.75, 0.75, 0.75, 0.8, 0.7, 0.7]
w_vals = [0.95 * w_to,  0.94 * w_to, 0.93 * w_to, 0.92 * w_to, 1 * w_to, 1 * w_to]
N_correction = [1, 1, 1, 1, 1, 2]

weight = w_to                          

wp_TOSL30, wp_TOSL50, wp_TO5K50, wp_cruise, wp_takeoffclimb, wp_transegmentclimb, \
    wp_secsegmentclimb, wp_enrouteclimb, wp_balkedAEO, wp_balkedOEI, \
        wp_celing, ws_L = get_PS_Plot(weight, precision, lbf_lbm, hp, WS, nu_p, AR, w_to, w_L, BFL, rho_SL, rho, 
                  C_lmax_TO, k_1, k_2, C_d0_TO, mu_g, S_FL, S_L, S_a, C_lmax_L, 
                  rho_cruise, V_cruise, C_D0_cruise, e_cruise, G_ceiling, e, 
                  rho_ceiling, C_Lmax, G_climbvals, k_svals, C_D0vals, 
                  C_Lmax_vals, labels, e_vals, w_vals, N_correction)