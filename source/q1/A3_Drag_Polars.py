# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 10:03:37 2023

@author: henry
"""
import numpy as np
import matplotlib.pyplot as plt
 
def Get_Drag_Polar(AR, Span, Wing_area, Max_Takeoff_W, c_f, c, d, phase):
    
    phase_dict = {"Clean" : 0,
                  "Takeoff flaps" : 1,
                  "Landing flaps" : 2,
                  "Landing gear" : 3
                  }
    i = phase_dict[phase]
    
    #Values from Roskam Vol 1 Table 3.6
    delta_Cd0_vals = (0, 0.02, 0.075, 0.025)
    e_vals = (0.8, 0.75, 0.7, 0.88)
    
    #Calulating Wetted Area
    S_wet = 10 ** c * Max_Takeoff_W ** d
    
    #Calulating Zero Lift Drag
    c_D0 = c_f * S_wet / Wing_area
    
    #C_L Values
    CL_vals = (np.linspace(-2,2,50), np.linspace(-2,2,30), np.linspace(-2.6,2.6,50), np.linspace(-2.6,2.6,50))
    
    C_L = CL_vals[i]
    
    #Calculating Drag
    C_D = ( c_D0 + delta_Cd0_vals[i] ) + C_L ** 2 / (np.pi * AR * e_vals[i])
    
    C_D0 = c_D0 + delta_Cd0_vals[i]
    
    return C_D, C_L, C_D0

#Calculating C_D and C_L values and plotting (AT-V1)
def plot_Drag_Polar(AR, Span, Wing_area, Max_Takeoff_W, c_f, c, d, title):
    plt.figure(figsize = (8,4) )
    plt.title(title)
    plt.xlabel("$C_D$")
    plt.ylabel("$C_L$")
    plt.xlim([0, 0.7])
    C_D_Clean, C_L_Clean, C_D0_Clean = Get_Drag_Polar(AR, Span, Wing_area, Max_Takeoff_W, c_f, c, d, phase = "Clean")
    plt.plot(C_D_Clean, C_L_Clean, label = "Clean", linestyle = "-", linewidth = 2)
    
    C_D_Takeoff, C_L_Takeoff, C_D0_Takeoff = Get_Drag_Polar(AR, Span, Wing_area, Max_Takeoff_W, c_f, c, d, phase = "Takeoff flaps")
    plt.plot(C_D_Takeoff, C_L_Takeoff, label = "Takeoff", linestyle = "-", linewidth = 2)
    
    C_D_Landing_flaps, C_L_Landing_flaps, C_D0_Landing_flaps = Get_Drag_Polar(AR, Span, Wing_area, Max_Takeoff_W, c_f, c, d, phase = "Landing flaps")
    plt.plot(C_D_Landing_flaps, C_L_Landing_flaps, label = "Landing Flap", linestyle = "-", linewidth = 2)
    
    C_D_Landing_gear, C_L_Landing_gear, C_D0_Landing_gear = Get_Drag_Polar(AR, Span, Wing_area, Max_Takeoff_W, c_f, c, d, phase = "Landing gear")
    plt.plot(C_D_Landing_gear, C_L_Landing_gear, label = "Landing Gear", linestyle = "-", linewidth = 2)
    
    plt.legend(loc='best')
    plt.grid()
    plt.show()
    
    print("---------------------------------------------")
    print("C_D0 Clean: ", round(C_D0_Clean, 3))
    print("C_D0 Takeoff: ", round(C_D0_Takeoff, 3))
    print("C_D0 Landing Flaps: ", round(C_D0_Landing_flaps, 3))
    print("C_D0 Landing Gear: ", round(C_D0_Landing_gear, 3))
    print("---------------------------------------------")
    
c = -0.0866                     #Roskam Vol 1 Table 3.5 (For a regional Turboprop)
d = 0.8099                      #Roskam Vol 1 Table 3.5 (For a regional Turboprop)
c_f = 0.0026                    #Raymer 2012 Table 12.3

# Setting Variables From OpenVSP (VT-V1)
AR = 10.06133                   #Aspect Ratio
Span = 96.428                   #Wing Span (ft)
Wing_area = 805.06             #Wing Area (ft^2)
Max_Takeoff_W = 82561.08        #Max Takeoff Weight (lbs)
title = "AT-V1"
print("C_D0 for ", title)
plot_Drag_Polar(AR, Span, Wing_area, Max_Takeoff_W, c_f, c, d, title)

# Setting Variables From OpenVSP (VT-V2)
AR = 7.78                   #Aspect Ratio
Span = 65                #Wing Span (ft)
Wing_area = 542.86             #Wing Area (ft^2)
Max_Takeoff_W = 100963        #Max Takeoff Weight (lbs)
title = "AT-V2"
print("C_D0 for ", title)
plot_Drag_Polar(AR, Span, Wing_area, Max_Takeoff_W, c_f, c, d, title)

# Setting Variables From OpenVSP (VT-V3)
AR = 6.28                   #Aspect Ratio
Span = 90                #Wing Span (ft)
Wing_area = 290             #Wing Area (ft^2)
Max_Takeoff_W = 100963        #Max Takeoff Weight (lbs)
title = "AT-V3"
print("C_D0 for ", title)
plot_Drag_Polar(AR, Span, Wing_area, Max_Takeoff_W, c_f, c, d, title)