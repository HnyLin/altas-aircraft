#By Henry Lin
#4/8/2023
#Fuel Burn Calculator

import numpy as np
from scipy.optimize import least_squares
import pandas as pd

#Drag Polar Calculation

def Get_Drag_Polar(AR, Span, Wing_area, MTOW, c_f, c, d, phase):
    
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
    K = K_vals[i]

    #Calulating Wetted Area
    S_wet = 10 ** c * MTOW ** d
    
    #Calulating Zero Lift Drag
    c_D0 = c_f * S_wet / Wing_area
    
    #C_L Values
    CL_vals = (np.linspace(-2,2,50), np.linspace(-2,2,30), np.linspace(-2.6,2.6,50), np.linspace(-2.6,2.6,50))
    
    C_L = CL_vals[i]
    
    #Calculating Drag
    C_D = ( c_D0 + delta_Cd0_vals[i] ) + C_L ** 2 / (np.pi * AR * e_vals[i])
    
    C_D0 = c_D0 + delta_Cd0_vals[i]
    
    return C_D, C_L, C_D0, K

def get_Drag_Coeffiecents(AR, Span, Wing_area, MTOW, c_f, c, d):

    C_D_Clean, C_L_Clean, C_D0_Clean, K_Clean = Get_Drag_Polar(AR, Span, Wing_area, MTOW, c_f, c, d, phase = "Clean") 

    C_D_Takeoff, C_L_Takeoff, C_D0_Takeoff, K_Takeoff = Get_Drag_Polar(AR, Span, Wing_area, MTOW, c_f, c, d, phase = "Takeoff flaps")
    
    C_D_Landing_flaps, C_L_Landing_flaps, C_D0_Landing_flaps, K_Landing_flaps = Get_Drag_Polar(AR, Span, Wing_area, MTOW, c_f, c, d, phase = "Landing flaps")
    
    C_D_Landing_gear, C_L_Landing_gear, C_D0_Landing_gear, K_Landing_gear = Get_Drag_Polar(AR, Span, Wing_area, MTOW, c_f, c, d, phase = "Landing gear")
    
    return C_D0_Clean, K_Clean

#Fuel Fraction Calculator
def Fuel_Fraction_Calculator(MTOW, MPOW, SFC, R, segments, eta, h_cruise, V_cruise):
    
    #Calculating Start, Warm-up, and Taxi Fuel Burn
    #Based Upon Assumption of Idling for 15 minutes w/ ideal being 5% of Max Power
    #SFC units lbm/(hp*hr)
    idle_POW = 0.05 * MPOW
    SWT_fuel_burn = (1 - hybridization_factors[0]) * SFC * 15/60 * idle_POW          #Units lbm
    SWT_fuel_weight = SWT_fuel_burn * 32.17
    print("SWT Fuel Burn (lbf): ", SWT_fuel_weight)

    W_SWT = MTOW - SWT_fuel_weight
    print("W_SWT (lbf): ", W_SWT)

    #Calculating Takeoff Fuel Fraction
    #Assuming 1min at Max Power
    ff_takeoff = 1 - 1 / 60 * SFC / eta * ( MPOW / W_SWT ) * 32.17 #eta is prop efficency
    print("ff_takeoff", ff_takeoff)

    W_Takeoff = W_SWT * ff_takeoff * (1 - hybridization_factors[1])
    print("W_Takeoff (lbf): ", W_Takeoff)

    Takeoff_fuel_weight = W_SWT - W_Takeoff
    print("Takeoff Fuel Weight (lbf): ", Takeoff_fuel_weight)

    #Calculating Climb Fuel Fractions (Multi-Segment Approach)
    ff_vals_climb = np.ones(segments)

    weight_vals_climb = np.ones(segments)
    weight_vals_climb[0] = W_Takeoff
    
    thrust_weight_vals_climb = np.ones(segments)

    velocity_vals_climb = np.ones(segments)

    C_L_vals_climb = np.ones(segments)

    C_D_vals_climb = np.ones(segments)

    D_vals_climb = np.ones(segments)

    delta_h_e = np.ones(segments - 1)

    temp = np.ones(segments)

    rho_interp = [0.0765, 0.0565, 0.0408, 0.0287, 0.0189]
    h_interp = [0, 10000, 20000, 30000, 40000]

    h_vals = np.linspace(0, h_cruise, segments)
    rho_vals = np.interp(h_vals, h_interp, rho_interp)

    for i in range(len(ff_vals_climb)):
        thrust_weight_vals_climb[i] = eta / V_cruise * MPOW / weight_vals_climb[i]

        MTOW = weight_vals_climb[i]
        C_D0_Clean, K_Clean = get_Drag_Coeffiecents(AR, Span, Wing_area, MTOW, c_f, c, d)

        velocity_vals_climb[i] = np.sqrt( weight_vals_climb[i] / Wing_area / ( 3 * rho_vals[i] * C_D0_Clean) * ( thrust_weight_vals_climb[i] + np.sqrt( thrust_weight_vals_climb[i]**2 + 12 * C_D0_Clean * K_Clean ) ) )

        C_L_vals_climb[i] = 2 * weight_vals_climb[i] / ( rho_vals[i] * velocity_vals_climb[i] ** 2 * Wing_area )

        C_D_vals_climb[i] = C_D0_Clean + K_Clean * C_L_vals_climb[i] ** 2

        D_vals_climb[i] = rho_vals[i] * velocity_vals_climb[i] ** 2 / 2 * Wing_area * C_D_vals_climb[i]

        temp[i] = h_vals[i] + velocity_vals_climb[i]**2 / (2 * 32.17)


    print("Climb Velcoity Vals: ", velocity_vals_climb)
    print("CL Values Climb: ", C_L_vals_climb)

    #Calculating Cruise Fuel Fraction
    range_intervals = np.linspace(0, R, segments)

    #Calculating coeffficent of lift
    #CL_seg = 2 * W_cruise[i]

    #Calculating Lift to Drag Ratio

    #Calculating Descent and Landing (Historical Data)
    ff_descent = 0.990
    ff_landing = 0.995

    return

c = -0.0866                     #Roskam Vol 1 Table 3.5 (For a regional Turboprop)
d = 0.8099                      #Roskam Vol 1 Table 3.5 (For a regional Turboprop)
c_f = 0.0026                    #Raymer 2012 Table 12.3

SFC = 0.4                       #Metabook (Mattingly 1996 Fig 1.17b) lbm / (hp * hr)
eta = 0.9                       #Propeller Efficency?

# Setting Variables From OpenVSP (VT-V1)
AR = 10.06133                   #Aspect Ratio
Span = 96.428                   #Wing Span (ft)
Wing_area = 805.06              #Wing Area (ft^2)

MTOW = 82561.08                 #Max Takeoff Weight (lbs)
MPOW = 7000                     #Hp Check Value!!!!!!!!!!!
R = 500 * 6076.12               #Range (ft)
h_cruise = 25000                #Cruising Altitude (ft)!!!!!!
V_cruise = 350 * 1.688 

segments = 20
hybridization_factors = (0.25, 0, 0, 0.5, 0.5)

# C_D0_Clean, C_D_Takeoff, C_D0_Landing_flaps, C_D0_Landing_gear, C_D0_Landing_flaps, K_Clean, K_Takeoff, K_Landing_flaps, K_Landing_gear \
#       = get_Drag_Coeffiecents(AR, Span, Wing_area, MTOW, c_f, c, d)

Fuel_Fraction_Calculator(MTOW, MPOW, SFC, R, segments, eta, h_cruise, V_cruise)