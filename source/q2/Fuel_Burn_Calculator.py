#By Henry Lin
#4/8/2023
#Fuel Burn Calculator

import numpy as np
from scipy.optimize import least_squares
import pandas as pd

#Drag Polar Calculation
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
    
    #Calculating K (Induced Drag Model)
    K_vals = 1 / (np.pi * AR * np.array(e_vals))

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
    
    return C_D, C_L, C_D0, K_vals

#Calculating C_D and C_L values and plotting (AT-V1)
def get_Drag_Coeffiecents(AR, Span, Wing_area, Max_Takeoff_W, c_f, c, d):

    C_D_Clean, C_L_Clean, C_D0_Clean, K_Clean = Get_Drag_Polar(AR, Span, Wing_area, Max_Takeoff_W, c_f, c, d, phase = "Clean") 

    C_D_Takeoff, C_L_Takeoff, C_D0_Takeoff, K_Takeoff = Get_Drag_Polar(AR, Span, Wing_area, Max_Takeoff_W, c_f, c, d, phase = "Takeoff flaps")
    
    C_D_Landing_flaps, C_L_Landing_flaps, C_D0_Landing_flaps, K_Landing_flaps = Get_Drag_Polar(AR, Span, Wing_area, Max_Takeoff_W, c_f, c, d, phase = "Landing flaps")
    
    C_D_Landing_gear, C_L_Landing_gear, C_D0_Landing_gear, K_Landing_gear = Get_Drag_Polar(AR, Span, Wing_area, Max_Takeoff_W, c_f, c, d, phase = "Landing gear")
    
    return C_D0_Clean, C_D_Takeoff, C_D0_Landing_flaps, C_D0_Landing_gear, C_D0_Landing_flaps, K_Clean, K_Takeoff, K_Landing_flaps, K_Landing_gear

 #Fuel Fraction Calculator
def Fuel_Fraction_Calculator(MTOW, MPOW, SFC, L_D, R, segments, V_cruise, C_D0, K, AR, e):

    #Calculating Start, Warm-up, and Taxi Fuel Burn
    #Based Upon Assumption of Idling for 15 minutes w/ ideal being 5% of Max Power
    #SFC units lbm/(hp*hr)
    idle_POW = 0.05 * MPOW
    SWT_fuel_burn = SFC * 15/60 * idle_POW          #Units lbm
    SWT_fuel_weight = SWT_fuel_burn * 32.17

    #Calculating Cruise Fuel Fraction
    range_intervals = np.linspace(0, R, segments)


    #Calculating Descent and Landing (Historical Data)
    ff_descent = 0.990
    ff_landing = 0.995

c = -0.0866                     #Roskam Vol 1 Table 3.5 (For a regional Turboprop)
d = 0.8099                      #Roskam Vol 1 Table 3.5 (For a regional Turboprop)
c_f = 0.0026                    #Raymer 2012 Table 12.3

# Setting Variables From OpenVSP (VT-V1)
AR = 10.06133                   #Aspect Ratio
Span = 96.428                   #Wing Span (ft)
Wing_area = 805.06             #Wing Area (ft^2)
Max_Takeoff_W = 82561.08        #Max Takeoff Weight (lbs)
C_D0_Clean, C_D_Takeoff, C_D0_Landing_flaps, C_D0_Landing_gear, C_D0_Landing_flaps, K_Clean, K_Takeoff, K_Landing_flaps, K_Landing_gear \
      = get_Drag_Coeffiecents(AR, Span, Wing_area, Max_Takeoff_W, c_f, c, d)


