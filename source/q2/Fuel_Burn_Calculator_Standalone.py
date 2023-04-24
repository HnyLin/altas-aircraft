#By Henry Lin
#4/8/2023
#Fuel Burn Calculator V2

'''
Description:
Inputs:
AR: Aspect Ratio (int)
Wing_area: Wing Area (int) [ft^2]
c_f: Roskam Skin Friction Coefficent (0.0026 for a Reg Turboprop)
c: Roskam Historical Regression Vol 1 Table 3.5 (-0.0866 for a Reg Turboprop)
d: Roskam Historical Regression Vol 1 Table 3.5 (0.8099 for a Reg Turboprop)
MTOW: Maximum Takeoff Weight (int) [lbf]
MPOW: Maximum Power (int) [hp]
SFC: Specific Fuel Consuption  (int) [lbm / (hp*hr)] Recommended Value: 0.4 (Metabook, Mattinglu 1996 Fig 1.17b)
R: Range (int) [ft]
segments: Number of Segments (increase for increased accuracy) (int)
eta: Prop Efficency (int)
h_cruise: Cruising Altitude (int) [ft]
V_cruise: Cruising Velocity (int) [ft/s]
hybridization_factors: Takes a List of 6 Values 0 (gas) 1 (electric) order: Start Warmup Taxi, Takeoff, Climb, Cruise, Descent, Landing

Outputs:
6 Values for Fuel Burn [lbs], Function Prints Fuel Burn for each Phase
'''

import numpy as np

#Drag Polar Calculation
def HFCA_to_Battery(fuel_weight):
    fuel_mass = fuel_weight / 32.17
    SED_JetA1 = 43.1 * 429.9                    #Btu/lb
    SED_Battery = 500 * 3600 / 10**6 * 429.9    #Btu/lb
    battery_mass = fuel_mass * SED_JetA1 / SED_Battery
    battery_weight = battery_mass * 32.17

    return battery_weight

def Get_Drag_Polar(AR, Wing_area, MTOW, c_f, c, d, phase):
    
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

def get_Drag_Coeffiecents(AR,  Wing_area, MTOW, c_f, c, d):

    C_D_Clean, C_L_Clean, C_D0_Clean, K_Clean = Get_Drag_Polar(AR, Wing_area, MTOW, c_f, c, d, phase = "Clean") 

    C_D_Takeoff, C_L_Takeoff, C_D0_Takeoff, K_Takeoff = Get_Drag_Polar(AR, Wing_area, MTOW, c_f, c, d, phase = "Takeoff flaps")
    
    C_D_Landing_flaps, C_L_Landing_flaps, C_D0_Landing_flaps, K_Landing_flaps = Get_Drag_Polar(AR, Wing_area, MTOW, c_f, c, d, phase = "Landing flaps")
    
    C_D_Landing_gear, C_L_Landing_gear, C_D0_Landing_gear, K_Landing_gear = Get_Drag_Polar(AR, Wing_area, MTOW, c_f, c, d, phase = "Landing gear")
    
    return C_D0_Clean, K_Clean

#Fuel Fraction Calculator
def Fuel_Fraction_Calculator(AR, Wing_area, c_f, c, d, MTOW, MPOW, SFC, R, segments, eta, h_cruise, V_cruise, hybridization_factors):
    
    #Calculating Start, Warm-up, and Taxi Fuel Burn
    #Based Upon Assumption of Idling for 15 minutes w/ ideal being 5% of Max Power
    #SFC units lbm/(hp*hr)
    idle_POW = 0.05 * MPOW
    SWT_fuel_mass = (1 - hybridization_factors[0]) * SFC * 15/60 * idle_POW         #Units lbm

    SWT_hybrid_weight = hybridization_factors[0] * SFC * 15/60 * idle_POW * 32.17   #Units lbf
    SWT_battery_weight = HFCA_to_Battery(SWT_hybrid_weight)

    SWT_fuel_burn = SWT_fuel_mass * 32.17
    print("SWT Fuel Burn (lbf): ", SWT_fuel_burn)
    print("SWT Battery Weight (lbf): ", SWT_battery_weight)
    W_SWT = MTOW - SWT_fuel_burn

    #Calculating Takeoff Fuel Fraction
    #Assuming 1min at Max Power
    ff_takeoff = 1 - 1 / 60 * SFC / eta * ( MPOW / W_SWT ) * 32.17 #eta is prop efficency

    W_Takeoff = W_SWT * ff_takeoff 

    Takeoff_fuel_burn = (W_SWT - W_Takeoff) * (1-hybridization_factors[1])

    Takeoff_hybrid_weight = (W_SWT - W_Takeoff) * (hybridization_factors[1])
    Takeoff_battery_weight = HFCA_to_Battery(Takeoff_hybrid_weight)

    print("Takeoff Fuel Weight (lbf): ", Takeoff_fuel_burn)
    print("Takeoff Battery Weight (lbf): ", Takeoff_battery_weight)

    #Calculating Climb Fuel Fractions (Multi-Segment Approach)
    ff_vals_climb = np.ones(segments-1)

    weight_vals_climb = np.ones(segments)
    weight_vals_climb[0] = W_Takeoff
    weight_vals_climb[1] = weight_vals_climb[0]
    
    thrust_weight_vals_climb = np.ones(segments)

    velocity_vals_climb = np.ones(segments)

    C_L_vals_climb = np.ones(segments)

    C_D_vals_climb = np.ones(segments)

    D_vals_climb = np.ones(segments)

    he_vals_climb = np.ones(segments)

    delta_he_vals_climb = np.ones(segments - 1)

    rho_interp = [0.0765, 0.0565, 0.0408, 0.0287, 0.0189]
    h_interp = [0, 10000, 20000, 30000, 40000]

    h_vals = np.linspace(0, h_cruise, segments)
    rho_vals = np.interp(h_vals, h_interp, rho_interp)

    c_t = 0.4                                       #Thrust Specific Fuel Consumption

    #Calculating Intial Condition
    i = 0
    thrust_weight_vals_climb[i] = eta / V_cruise * MPOW / weight_vals_climb[i] * 550        #Horsepower Conversion

    MTOW = weight_vals_climb[i]
    C_D0_Clean, K_Clean = get_Drag_Coeffiecents(AR,  Wing_area, MTOW, c_f, c, d)

    velocity_vals_climb[i] = np.sqrt( 32.17 * weight_vals_climb[i] / Wing_area / ( 3 * rho_vals[i] * C_D0_Clean) * ( thrust_weight_vals_climb[i] + np.sqrt( thrust_weight_vals_climb[i]**2 + 12 * C_D0_Clean * K_Clean ) ) )

    C_L_vals_climb[i] = 2 * weight_vals_climb[i] / ( rho_vals[i] * velocity_vals_climb[i] ** 2 * Wing_area ) * 32.17        #lbf_lbm converison

    C_D_vals_climb[i] = C_D0_Clean + K_Clean * C_L_vals_climb[i] ** 2

    D_vals_climb[i] = rho_vals[i] * velocity_vals_climb[i] ** 2 / 2 * Wing_area * C_D_vals_climb[i] / 32.17                 #lbf_lbm conversion

    he_vals_climb[i] = h_vals[i] + velocity_vals_climb[i]**2 / (2 * 32.17)

    for i in range(1, segments-1):
        thrust_weight_vals_climb[i] = eta / V_cruise * MPOW / weight_vals_climb[i] * 550    #Horsepower Conversion

        MTOW = weight_vals_climb[i]
        C_D0_Clean, K_Clean = get_Drag_Coeffiecents(AR,  Wing_area, MTOW, c_f, c, d)

        velocity_vals_climb[i] = np.sqrt( 32.17 * weight_vals_climb[i] / Wing_area / ( 3 * rho_vals[i] * C_D0_Clean) * ( thrust_weight_vals_climb[i] + np.sqrt( thrust_weight_vals_climb[i]**2 + 12 * C_D0_Clean * K_Clean ) ) )

        C_L_vals_climb[i] = 2 * weight_vals_climb[i] / ( rho_vals[i] * velocity_vals_climb[i] ** 2 * Wing_area ) * 32.17    #lbf_lbm conversion

        C_D_vals_climb[i] = C_D0_Clean + K_Clean * C_L_vals_climb[i] ** 2

        D_vals_climb[i] = rho_vals[i] * velocity_vals_climb[i] ** 2 / 2 * Wing_area * C_D_vals_climb[i] / 32.17             #lbf_lbm conversion

        he_vals_climb[i] = h_vals[i] + velocity_vals_climb[i]**2 / (2 * 32.17)

        delta_he_vals_climb[i] = he_vals_climb[i] - he_vals_climb[i-1]

        ff_vals_climb[i] = np.exp( -c_t / 3600 * delta_he_vals_climb[i] / ( velocity_vals_climb[i] * ( 1 - D_vals_climb[i] / ( weight_vals_climb[i] * thrust_weight_vals_climb[i] ) ) ) )
        
        weight_vals_climb[i+1] = weight_vals_climb[i] * ff_vals_climb[i]

    weight_climb = weight_vals_climb[-1]                        #Weight of Plane After Climb (lbf)
    climb_fuel_burn = weight_vals_climb[0] - weight_climb       #Weight of Fuel Burned During Climb (lbf)

    hybrid_fuel_weight = climb_fuel_burn * hybridization_factors[2]
    climb_battery_weight = HFCA_to_Battery(hybrid_fuel_weight)

    climb_fuel_burn = climb_fuel_burn * (1 - hybridization_factors[2])
    # print("Climb Exit Velocity (ft/s): ", velocity_vals_climb[-2])
    print("Climb Fuel Burn (lbf): ", climb_fuel_burn)
    print("Climb Battery Weight (lbf): ", climb_battery_weight)

    #Calculating Cruise Fuel Fraction
    #Allocating Space
    weight_vals_cruise = np.ones(segments)
    weight_vals_cruise[0] = weight_climb
    CL_vals_cruise = np.ones(segments)
    LD_vals_cruise = np.ones(segments)
    ff_vals_cruise = np.ones(segments)

    rho_cruise = np.interp(28000, h_interp, rho_interp)

    #Calculating coeffficent of lift
    for i in range(segments-1):
        CL_vals_cruise[i] = 2 * weight_vals_cruise[i] / ( rho_cruise * V_cruise**2 * Wing_area ) * 32.17     #lbm_lbf conversion

        MTOW = weight_vals_cruise[i]
        C_D0_Clean, K_Clean = get_Drag_Coeffiecents(AR,  Wing_area, MTOW, c_f, c, d)

        LD_vals_cruise[i] = CL_vals_cruise[i] / ( C_D0_Clean + K_Clean * CL_vals_cruise[i]**2 )

        ff_vals_cruise[i] = np.exp( -R / segments * c_t / 3600 / ( V_cruise * LD_vals_cruise[i] ) )

        weight_vals_cruise[i+1] = ff_vals_cruise[i] * weight_vals_cruise[i]

    cruise_fuel_burn = weight_vals_cruise[0] - weight_vals_cruise[-1]

    hybrid_fuel_weight = cruise_fuel_burn * hybridization_factors[3]
    cruise_battery_weight = HFCA_to_Battery(hybrid_fuel_weight)
    cruise_fuel_burn = cruise_fuel_burn * (1 - hybridization_factors[3])
    print("Cruise Fuel Burn (lbf): ", cruise_fuel_burn)
    print("Cruise Battery Weight (lbf): ", cruise_battery_weight)

    #Calculating Descent and Landing (Historical Data)
    weight_descent_entry = weight_vals_cruise[-1]
    ff_descent = 0.990
    weight_descent_exit = ff_descent * weight_descent_entry
    desecent_fuel_burn = weight_descent_entry - weight_descent_exit

    hybrid_fuel_weight = desecent_fuel_burn * hybridization_factors[4]
    descent_battery_weight = HFCA_to_Battery(hybrid_fuel_weight)

    desecent_fuel_burn = desecent_fuel_burn * (1 - hybridization_factors[4])
    print("Descent Fuel Burn (lbf): ", desecent_fuel_burn)
    print("Descent Fuel Weight (lbf): ", descent_battery_weight)

    ff_landing = 0.995
    weight_landing_exit = ff_landing * weight_descent_exit
    landing_fuel_burn = weight_descent_exit - weight_landing_exit

    hybrid_fuel_weight = landing_fuel_burn * hybridization_factors[5]
    landing_battery_weight = HFCA_to_Battery(hybrid_fuel_weight)

    landing_fuel_burn = landing_fuel_burn * (1 - hybridization_factors[5])
    print("Landing Fuel Burn (lbf): ", landing_fuel_burn)
    print("Landing Battery Weight (lbf): ", landing_battery_weight)

    total_fuel_burn = SWT_fuel_burn + Takeoff_fuel_burn + climb_fuel_burn + cruise_fuel_burn + desecent_fuel_burn + landing_fuel_burn
    print("Total Fuel Burn (lbf): ", total_fuel_burn)

    total_battery_weight = SWT_battery_weight + Takeoff_battery_weight + climb_battery_weight + cruise_battery_weight + descent_battery_weight + landing_battery_weight
    print("Total Battery Weight (lbf): ", total_battery_weight)

    total_hybrid_weight = total_battery_weight + total_fuel_burn
    print("Total Hybrid Weight (lbf): ", total_hybrid_weight)

    return SWT_fuel_burn, Takeoff_fuel_burn, climb_fuel_burn, cruise_fuel_burn, desecent_fuel_burn, landing_fuel_burn, total_fuel_burn, total_battery_weight, total_hybrid_weight