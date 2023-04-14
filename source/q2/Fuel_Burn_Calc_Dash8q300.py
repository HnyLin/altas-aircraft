#By: Henry Lin
#4/14/2023

import numpy as np
import Fuel_Burn_Calculator_Standalone as FBC

#Inputs for a Dash 8-q300
AR = 13.39
Span = 90
Wing_area = 604.9
MTOW = 43000
MPOW = 4760
R = 500 * 6076.12 
h_cruise = 25000
V_cruise = 287 * 1.688
hybridization_factors = [0, 0, 0, 0, 0, 0]

c = -0.0866                     #Roskam Vol 1 Table 3.5 (For a regional Turboprop)
d = 0.8099                      #Roskam Vol 1 Table 3.5 (For a regional Turboprop)
c_f = 0.0026                    #Raymer 2012 Table 12.3

SFC = 0.4                       #Metabook (Mattingly 1996 Fig 1.17b) lbm / (hp * hr)
eta = 0.9                       #Propeller Efficency

segments = 20

SWT_fuel_burn, Takeoff_fuel_burn, climb_fuel_burn, cruise_fuel_burn, desecent_fuel_burn, landing_fuel_burn, D8total_fuel_burn, D8total_battery_weight, D8total_hybrid_weight = \
    FBC.Fuel_Fraction_Calculator(AR, Wing_area, c_f, c, d, MTOW, MPOW, SFC, R, segments, eta, h_cruise, V_cruise, hybridization_factors)

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
h_cruise = 20000                #Cruising Altitude (ft)!!!!!!
V_cruise = 350 * 1.688 

segments = 20

#Start Warmup Taxi, Takeoff, Climb, Cruise, Descent, Landing (Loitter Unavaliable)
hybridization_factors = (0.2, 0.2, 0, 0, 0.5, 0.5)

SWT_fuel_burn, Takeoff_fuel_burn, climb_fuel_burn, cruise_fuel_burn, desecent_fuel_burn, landing_fuel_burn, total_fuel_burn, total_battery_weight, total_hybrid_weight = \
    FBC.Fuel_Fraction_Calculator(AR, Wing_area, c_f, c, d, MTOW, MPOW, SFC, R, segments, eta, h_cruise, V_cruise, hybridization_factors)

print("-----------------------------------------")
print("Dash 8 Fuel Burn: ", D8total_fuel_burn)
print("ATLAS Fuel Burn: ", total_fuel_burn)