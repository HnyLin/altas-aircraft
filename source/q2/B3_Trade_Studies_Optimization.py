#B3 Trade Studies
'''
Base code is from "B1_MTOW_Refined.py"
Idea is to investigate changes to MTOW and Block Fuel Burn is different geometric designs were pursued.
In the current iteration of this code, we are investigating varying AR.
Variables to be changed in Maks code:
A, t_c_root, need to go back through B1
'''
import numpy as np
import matplotlib.pyplot as plt
import Fuel_Burn_Calculator_Standalone as FBC
import scipy.optimize as optimize
'''
I apologize to whoever gets to read this, go to the bottom to see loop. This script combines 
Henry's fuel/battery function and my empty weight function  to calculate, and iterate on, the MTOW.
Our takeoff power also iterates with each loop too.
'''
def calcEmptyWeight(W_TO, P_rshp, AR, t_c_root, S_w):
    '''
    Hard coding all constants in this function so you don't need to call them. May or may not be a good idea. 
    If one wants to change aspects to the design, they would need to edit the code within this function
    to reflect that design

    INPUTS:
    W_TO - Takeoff Weight (lbs)
    P_rshp - Shaft HP at takeoff (hp)
    AR - Wing Aspect Ratio
    t_c_root - thickness to chord ratio for wing
    S_w - wing ref area

    OUTPUTS:
    Component Weights
    '''
    W_dg = W_TO             #design weight, lbs
    N_z = 1.5*3.5           #ult load factor, raymor table 14.2
    #S_w = 805.06257         #trapezoidal wing area, ft^2
    #AR = 10.0613            #Wing Aspect Ratio
    #t_c_root = 0.15450      #Wing thickness-to-chord ratio 
    lam = .7525             #Taper ratio  
    Lam = 20                #Wing Sweep degrees, check rads
    S_csw = 2*((4.254*9.024*2) + (6.295*(5.449 + 5.802)))        #wing control surface area, upper and lower, 295.202474 ft^2

    #HT
    K_uht = 1.143           #coeff, NEED CHANGE TO 1.0, SEE RAYMER 
    F_w = 10
    B_h = 31.6
    S_ht = 202.18
    L_t = 24.5        #estimated length, wing MAC to tail MAC
    K_y = 0.3*L_t
    Lam_ht = 20
    A_h = 5.6
    S_e = 2*((1.350+3.020)*0.5*8.478)   #elevator area, 37.04886

    #VT
    Ht_Hv = 1
    S_vt = 190.79
    K_z = L_t
    Lam_vt = 40 #About a 40deg sweep
    A_v = 70.4
    t_c_root_vt =  0.16512 #check with tre

    #fuselage
    B_w = 96.4
    K_door = 1.12
    K_Lg = 1.12
    L = 73    #measured aribtrary length, excludes nose cap and radome cowling 
    S_f = 2093
    df =10*12   #diamaeter of fuselage, in
    D = 0.02*df + 1     # fuselage structrual depth in

    #LG_m
    K_mp = 1.15
    W_l = 0.043*W_dg    #Raymer approximation, Table 15.2
    N_gear = 3
    N_l = 1.5*N_gear
    L_m = 3.68*12   #open vsp measure, (in)
    N_mw = 4
    N_mss = 3
    V_stall = 205.913

    #LG_n
    K_np = 1.15
    L_n = 3.771*12  #open vsp measure, (in)
    N_nw = 2

    #enginecontrol
    N_en = 2
    L_ec = 31

    #flightcontrol
    N_f = 5
    N_m = 1
    S_r = 2*(2*(1.033+2.288)*0.5*10.691)    #rudder area, 71.009622
    S_cs = S_csw + S_r + S_e    #total control surface area, 403.260956
    I_yaw = 2494172.151 #simulated using const density in openvsp

    #fuelsys
    V_t = 128.604   #total volume
    V_i = 128.604   #integral fuel tank volume (wings)
    V_p = 128.604   #self-sealing, "protexted tank volume"
    N_t = 4         #number of tanks (4)

    #avionics
    W_uav = 1100 #keep for now, look into typical uninstalled avionics 

    #instruments
    K_r = 1.0
    K_tp = 0.793
    N_c = 4
    N_en = 2
    L_f = 81

    #hydraulics
    #already defined

    #Air Conditioning   #readdress later
    N_p = 53
    #V_pr =         #need internal pressurized volume, too complex for now

    #nacelle group
    K_ng = 1.15
    N_Lt = 14.167
    N_w = 2.354     #ft, EXTREMELY SMALL, MAY NEED RESIZING
    S_n = 52.454

    W_engine = P_rshp**(0.9306)*10**(-0.1205)  #Weight of engine, lbs

    K_p = 1.4
    K_tr = 1.0


    #=======================
    #STRUCTURAL WEIGHT
    #=======================

    #Wing Weight
    W_wing = 0.0051*(W_dg*N_z)**0.557*S_w**0.649*AR**0.5*(t_c_root)**-0.4*(1 + lam)**0.1*np.cos(Lam/180.0*np.pi)**-1.0*S_csw**0.1
    #Horizontal Tail Weight
    W_HT = 0.0379*K_uht*(1 + F_w/B_h)**-0.25*W_dg**0.639*N_z**0.10*S_ht**0.75*L_t**-1.0*K_y**0.704*np.cos(Lam_ht/180.0*np.pi)**-1.0*A_h**0.166*(1 + S_e/S_ht)**0.1
    #Vertical Tail Weight
    W_VT = 0.0026*(1 + Ht_Hv)**0.225*W_dg**0.556*N_z**0.536*L_t**(-0.5)*S_vt**0.5*K_z**0.875*np.cos(Lam_vt/180*np.pi)**(-1.0)*A_v**0.35*(t_c_root_vt)**(-0.5)

    #K_ws needed for W_fuse
    K_ws = 0.75*((1 + 2*lam)/(1 + lam))*(B_w*np.tan(Lam/L))

    #Fuselage Weight
    W_fuse = 0.3280*K_door*K_Lg*(W_dg*N_z)**0.5*L**0.25*S_f**0.302*(1 + K_ws)**0.04*(L/D)**0.10
    #Main Landing Gear Weight
    W_lg_main = 0.0106*K_mp*W_l**0.888*N_l**0.25*L_m**0.4*N_mw**0.321*N_mss**-0.5*V_stall**0.1
    #Nose Landing Gear Weight
    W_lg_nose = 0.032*K_np*W_l**0.646*N_l**0.2*L_n**0.5*N_nw**0.45

    #=======================
    #CONTROL SYSTEM WEIGHT
    #=======================
    #Engine Controls Weight
    W_encl = 5.0*N_en + 0.80*L_ec

    #Flight Controls Weight
    W_fc = 145.9*N_f**0.554*(1 + N_m/N_f)**-1.0*S_cs**0.20*(I_yaw*10**-6)**0.07

    #=======================
    #Systems Weight
    #=======================
    #Fuel System Weight
    W_fs = 2.405*V_t**0.606*(1 + V_i/V_t)**-1.0*(1 + V_p/V_t)*N_t**0.5

    #Avionics Weight
    W_av = 1.73*W_uav**0.983

    #Instruments Weight
    W_instr = 4.509*K_r*K_tp*N_c**0.541*N_en*(L_f + B_w)**0.5

    #Hydraulics or electronics Weight (tbd)
    W_hyd = 0.2673*N_f*(L_f + B_w)**0.937

    #Anti-icing
    W_ai = 0.002*W_dg

    #Air Conditionting
    #W_ac = 63.36*N_p**0.25*(V_pr/1000)**0.604*W_uav**0.10      #holding off for now, need a way to find pressurized volume

    #=======================
    #Propulsion Weight
    #=======================
    #Needed for nacelle group
    W_ec = 2.331*W_engine**0.901*K_p*K_tr
    W_em = 2*245.8 #weight, lbs, of like EMs,if time premits, look into max  Power requirement for our EM, and look for irl EMs that can satisfy https://skiesmag.com/news/electric-motor-manufacturer-magnix-set-to-conquer-aviation-market/

    #Nacelle Group Weight
    W_ng = .6724*K_ng*N_Lt**0.10*N_w**0.294*N_z**0.119*W_ec**0.611*N_en**0.984*S_n**0.224 + W_em #check later, may be only counting as single nacelle


    W_empty = W_ng + W_ai + W_hyd + W_instr + W_av + W_fs + W_fc + W_encl + W_lg_nose + W_lg_main + W_fuse +  W_HT + W_wing + W_VT


    #print('=============================')
    #print('Summary of Weight Results')
    #print('=============================')
    #print('Wing: %0.3f' % W_wing)
    #print('Vertical Tail: %0.3f' % W_VT)
    #print('Horizontal Tail: %0.3f' % W_HT)
    #print('Fuselage: %0.3f' % W_fuse)
    #print('Main Landing Gear: %0.3f' % W_lg_main)
    #print('Nose Landing Gear: %0.3f' % W_lg_nose)
    #print('Engine Controls: %0.3f' % W_encl)
    #print('Flight Controls: %0.3f' % W_fc)
    #print('Fuel System: %0.3f' % W_fs)
    #print('Avionics: %0.3f' % W_av)
    #print('Instruments: %0.3f' % W_instr)
    #print('Hydraulics: %0.3f' % W_hyd)
    #print('Anti-icing: %0.3f' % W_ai)
    #print('Nacelle Group: %0.3f' % W_ng)
    #print('=============================')
    #print('NOW HENRYS RESULTS!')
    #print('=============================')
    return W_empty
#====================================================
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
    #print("SWT Fuel Burn (lbf): ", SWT_fuel_burn)
    #print("SWT Battery Weight (lbf): ", SWT_battery_weight)
    W_SWT = MTOW - SWT_fuel_burn

    #Calculating Takeoff Fuel Fraction
    #Assuming 1min at Max Power
    ff_takeoff = 1 - 1 / 60 * SFC / eta * ( MPOW / W_SWT ) * 32.17 #eta is prop efficency

    W_Takeoff = W_SWT * ff_takeoff 

    Takeoff_fuel_burn = (W_SWT - W_Takeoff) * (1-hybridization_factors[1])

    Takeoff_hybrid_weight = (W_SWT - W_Takeoff) * (hybridization_factors[1])
    Takeoff_battery_weight = HFCA_to_Battery(Takeoff_hybrid_weight)

    #print("Takeoff Fuel Weight (lbf): ", Takeoff_fuel_burn)
    #print("Takeoff Battery Weight (lbf): ", Takeoff_battery_weight)

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
    # #print("Climb Exit Velocity (ft/s): ", velocity_vals_climb[-2])
    #print("Climb Fuel Burn (lbf): ", climb_fuel_burn)
    #print("Climb Battery Weight (lbf): ", climb_battery_weight)

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
    #print("Cruise Fuel Burn (lbf): ", cruise_fuel_burn)
    #print("Cruise Battery Weight (lbf): ", cruise_battery_weight)

    #Calculating Descent and Landing (Historical Data)
    weight_descent_entry = weight_vals_cruise[-1]
    ff_descent = 0.990
    weight_descent_exit = ff_descent * weight_descent_entry
    desecent_fuel_burn = weight_descent_entry - weight_descent_exit

    hybrid_fuel_weight = desecent_fuel_burn * hybridization_factors[4]
    descent_battery_weight = HFCA_to_Battery(hybrid_fuel_weight)

    desecent_fuel_burn = desecent_fuel_burn * (1 - hybridization_factors[4])
    #print("Descent Fuel Burn (lbf): ", desecent_fuel_burn)
    #print("Descent Fuel Weight (lbf): ", descent_battery_weight)

    ff_landing = 0.995
    weight_landing_exit = ff_landing * weight_descent_exit
    landing_fuel_burn = weight_descent_exit - weight_landing_exit

    hybrid_fuel_weight = landing_fuel_burn * hybridization_factors[5]
    landing_battery_weight = HFCA_to_Battery(hybrid_fuel_weight)

    landing_fuel_burn = landing_fuel_burn * (1 - hybridization_factors[5])
    #print("Landing Fuel Burn (lbf): ", landing_fuel_burn)
    #print("Landing Battery Weight (lbf): ", landing_battery_weight)

    total_fuel_burn = SWT_fuel_burn + Takeoff_fuel_burn + climb_fuel_burn + cruise_fuel_burn + desecent_fuel_burn + landing_fuel_burn
    #print("Total Fuel Burn (lbf): ", total_fuel_burn)

    total_battery_weight = SWT_battery_weight + Takeoff_battery_weight + climb_battery_weight + cruise_battery_weight + descent_battery_weight + landing_battery_weight
    #print("Total Battery Weight (lbf): ", total_battery_weight)

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Changed Battery Weight to Be Greatest Battery Weight of All Phases!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # total_hybrid_weight = total_battery_weight + total_fuel_burn
    # #print("Total Hybrid Weight (lbf): ", total_hybrid_weight)
    Battery_Weights = [SWT_battery_weight, Takeoff_battery_weight, climb_battery_weight, cruise_battery_weight, descent_battery_weight, landing_battery_weight]
    Battery_Weights.sort()

    total_hybrid_weight = total_fuel_burn + Battery_Weights[-1] * 1.2

    return SWT_fuel_burn, Takeoff_fuel_burn, climb_fuel_burn, cruise_fuel_burn, desecent_fuel_burn, landing_fuel_burn, total_fuel_burn, total_battery_weight, total_hybrid_weight
#================================================================================================================
'''
#Dash 8
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

segments = 30

SWT_fuel_burn, Takeoff_fuel_burn, climb_fuel_burn, cruise_fuel_burn, desecent_fuel_burn, landing_fuel_burn, D8total_fuel_burn, D8total_battery_weight, D8total_hybrid_weight = \
    FBC.Fuel_Fraction_Calculator(AR, Wing_area, c_f, c, d, MTOW, MPOW, SFC, R, segments, eta, h_cruise, V_cruise, hybridization_factors)

c = -0.0866                     #Roskam Vol 1 Table 3.5 (For a regional Turboprop)
d = 0.8099                      #Roskam Vol 1 Table 3.5 (For a regional Turboprop)
c_f = 0.0026                    #Raymer 2012 Table 12.3

SFC = 0.4                       #Metabook (Mattingly 1996 Fig 1.17b) lbm / (hp * hr)
eta = 0.9                       #Propeller Efficency?
'''

def tradeStudies(AR, t_c_root, Wing_area, V_cruise, h1, h2, h3, h4):
    '''
    Trade Studies Loop. Slight modification of B1_MTOW_Refined.py
    Takes input variables (for )
    INPUTS:
    AR - Wing Aspect Ratio
    t_c_root - maximum thickness-to-chord ratio (constant along wing)
    Wing_area - total wing area, ft^2
    V_cruise - cruise speed, knts
    h1, h2, h3, h4 - hybrid. factors for Warmup Taxi, Takeoff, Descent, Landing (climb and cruise are zero hybrid)

    OUTPUTS:
    Weight Breakdown
    Fuel Breakdown


    COMMENTS:
    1) May want to suppress some of the printed results, can be a lot when rerunning code

    2) Tried but commented out code that tries to evaluate change in one variable (ex. AR vs MTOW) for a more simple trade study
    
    3) Made wing span variable with in code, which is calucalated from AR and wing area
    '''
    c = -0.0866                     #Roskam Vol 1 Table 3.5 (For a regional Turboprop)
    d = 0.8099                      #Roskam Vol 1 Table 3.5 (For a regional Turboprop)
    c_f = 0.0026                    #Raymer 2012 Table 12.3

    SFC = 0.4                       #Metabook (Mattingly 1996 Fig 1.17b) lbm / (hp * hr)
    eta = 0.9                       #Propeller Efficency?

    # Setting Variables From OpenVSP (VT-V1)
    #AR = 10.06133                   #Aspect Ratio
    #Wing_area = 805.06              #Wing Area (ft^2)

    #Span = 96.428                   #Wing Span (ft)
    Span = np.sqrt(AR*Wing_area)

    MTOW = 82561.08                 #Max Takeoff Weight (lbs)
    MPOW = 7000                     #Hp Check Value!!!!!!!!!!!
    R = 500 * 6076.12               #Range (ft)
    h_cruise = 28000                #Cruising Altitude (ft)!!!!!!
    V_cruise = V_cruise * 1.688     #Convert V_cruise to ft/s

    segments = 20

    #Start Warmup Taxi, Takeoff, Climb, Cruise, Descent, Landing (Loitter Unavaliable)
    hybridization_factors = (h1, h2, 0, 0, h3, h4)

    #OTHER VARIABLES FOR LOOP
    W_P = 11.25     #lbf/hp
    W_crew_and_payload = 12660      #weight of crew, passengers, and payload, lbs

    #Loop setup
    tol = 1e-6
    dif = 1
    p = 0
    #MTOW_plot = MTOW
    while dif > tol:
    #while p <50:
        p = p+1

        W_empty= calcEmptyWeight(MTOW, MPOW, AR, t_c_root, Wing_area)

        SWT_fuel_burn, Takeoff_fuel_burn, climb_fuel_burn, cruise_fuel_burn, desecent_fuel_burn, landing_fuel_burn, total_fuel_burn, total_battery_weight, total_hybrid_weight = Fuel_Fraction_Calculator(AR, Wing_area, c_f, c, d, MTOW, MPOW, SFC, R, segments, eta, h_cruise, V_cruise, hybridization_factors)

        MTOW_new = W_empty + total_hybrid_weight + W_crew_and_payload
        dif = abs(MTOW_new - MTOW)

        #MTOW_plot[p] = MTOW_new
        MTOW = MTOW_new

        MPOW = MTOW/W_P

    #print('New MTOW is: ', MTOW_new)
    #print('New Power Req is:', MPOW)
    #print('Difference is: ', dif)
    #print('Iterations: ;', p)
    
    return MTOW_new, MPOW, total_fuel_burn
#================================================================================================================

#Test function
#tradeStudies(AR, t_c_root, Wing_area, V_cruise, h1, h2, h3, h4)

passengers = 50


#SCIPY Optimzation
def objective_function(params):
    AR, t_c_root, Wing_area, V_cruise, h1, h2, h3, h4 = params
    '''
    Trade Studies Loop. Slight modification of B1_MTOW_Refined.py
    Takes input variables (for )
    INPUTS:
    AR - Wing Aspect Ratio
    t_c_root - maximum thickness-to-chord ratio (constant along wing)
    Wing_area - total wing area, ft^2
    V_cruise - cruise speed, knts
    h1, h2, h3, h4 - hybrid. factors for Warmup Taxi, Takeoff, Descent, Landing (climb and cruise are zero hybrid)

    OUTPUTS:
    Weight Breakdown
    Fuel Breakdown


    COMMENTS:
    1) May want to suppress some of the printed results, can be a lot when rerunning code

    2) Tried but commented out code that tries to evaluate change in one variable (ex. AR vs MTOW) for a more simple trade study
    
    3) Made wing span variable with in code, which is calucalated from AR and wing area
    '''
    c = -0.0866                     #Roskam Vol 1 Table 3.5 (For a regional Turboprop)
    d = 0.8099                      #Roskam Vol 1 Table 3.5 (For a regional Turboprop)
    c_f = 0.0026                    #Raymer 2012 Table 12.3

    SFC = 0.4                       #Metabook (Mattingly 1996 Fig 1.17b) lbm / (hp * hr)
    eta = 0.9                       #Propeller Efficency?

    # Setting Variables From OpenVSP (VT-V1)
    #AR = 10.06133                   #Aspect Ratio
    #Wing_area = 805.06              #Wing Area (ft^2)

    #Span = 96.428                   #Wing Span (ft)
    Span = np.sqrt(AR*Wing_area)

    MTOW = 82561.08                 #Max Takeoff Weight (lbs)
    MPOW = 7000                     #Hp Check Value!!!!!!!!!!!
    R = 500 * 6076.12               #Range (ft)
    h_cruise = 28000                #Cruising Altitude (ft)!!!!!!
    V_cruise = V_cruise * 1.688     #Convert V_cruise to ft/s

    segments = 20

    #Start Warmup Taxi, Takeoff, Climb, Cruise, Descent, Landing (Loitter Unavaliable)
    hybridization_factors = (h1, h2, 0, 0, h3, h4)

    #OTHER VARIABLES FOR LOOP
    W_P = 11.25     #lbf/hp
    W_crew_and_payload = 12660      #weight of crew, passengers, and payload, lbs

    #Loop setup
    tol = 1e-6
    dif = 1
    p = 0
    #MTOW_plot = MTOW
    while dif > tol:
    #while p <50:
        p = p+1

        W_empty= calcEmptyWeight(MTOW, MPOW, AR, t_c_root, Wing_area)

        SWT_fuel_burn, Takeoff_fuel_burn, climb_fuel_burn, cruise_fuel_burn, desecent_fuel_burn, landing_fuel_burn, total_fuel_burn, total_battery_weight, total_hybrid_weight = Fuel_Fraction_Calculator(AR, Wing_area, c_f, c, d, MTOW, MPOW, SFC, R, segments, eta, h_cruise, V_cruise, hybridization_factors)

        MTOW_new = W_empty + total_hybrid_weight + W_crew_and_payload
        dif = abs(MTOW_new - MTOW)

        #MTOW_plot[p] = MTOW_new
        MTOW = MTOW_new

        MPOW = MTOW/W_P

    #print('New MTOW is: ', MTOW_new)
    #print('New Power Req is:', MPOW)
    #print('Difference is: ', dif)
    #print('Iterations: ;', p)
    
    return total_fuel_burn
#Setting Initial Guess
initial_guess = [10.06, 0.15450, 805.06, 350, 0.25, 0.25, 0.25, 0.25]

#Setting Bounds
bound_vals = ((8,20), (0.1, 1), (600, 1000), (280, 450), (0, 1), (0, 1), (0, 1), (0, 1))

#Optimize
result = optimize.minimize(objective_function, x0 = initial_guess, bounds = bound_vals, options= {'disp': True} )
print("Optimized Values")
print(result.x)

print("Optimum Fuel Weight (lbf): ", result.fun)
print("Optimum Fuel Weight Per Passenger (lbf): ", result.fun/50)
#================================================================================================================


#Trade Study One AR Sweep
resolution = 20
MTOW_vals = np.zeros(resolution)
MPOW_vals = np.zeros(resolution)
total_fuel_burn_vals = np.zeros(resolution)

AR_vals = np.linspace(8, 15, resolution)

counter = 0
for AR in AR_vals:
    MTOW_new, MPOW, total_fuel_burn = tradeStudies(AR, 0.15450, 805.06, 350, 0.2, 0.2, 0.5, 0.5)       #Will Produce results from Weights Refinement
    MTOW_vals[counter] = MTOW_new
    MPOW_vals[counter] = MPOW
    total_fuel_burn_vals[counter] = total_fuel_burn
    counter = counter + 1
    print("Loop Count:", counter)

#Plotting
plt.figure(figsize=(8,12))
plt.subplot(2, 1, 1)
plt.plot(AR_vals, MTOW_vals)
plt.xlabel("Aspect Ratio")
plt.ylabel("Maximum Takeoff Weight (lbf)")
plt.title("AR vs. Maximum Takeoff Weight")

plt.subplot(2, 1, 2)
plt.plot(AR_vals, total_fuel_burn_vals / passengers)
plt.xlabel("Aspect Ratio")
plt.ylabel("Fuel Burn per Passenger Weight (lbf)")
plt.title("AR vs. Fuel Burn per Passenger")
plt.show()

#Trade Study Two t_c_root Sweep
MTOW_vals = np.zeros(resolution)
MPOW_vals = np.zeros(resolution)
total_fuel_burn_vals = np.zeros(resolution)

t_c_root_vals = np.linspace(0.1, 0.4, resolution)

counter = 0
for t_c_root  in t_c_root_vals:
    MTOW_new, MPOW, total_fuel_burn = tradeStudies(10.06, t_c_root, 805.06, 350, 0.2, 0.2, 0.5, 0.5)       #Will Produce results from Weights Refinement
    MTOW_vals[counter] = MTOW_new
    MPOW_vals[counter] = MPOW
    total_fuel_burn_vals[counter] = total_fuel_burn
    counter = counter + 1
    print("Loop Count:", counter)

#Plotting
plt.figure(figsize=(8,12))
plt.subplot(2, 1, 1)
plt.plot(t_c_root_vals, MTOW_vals)
plt.xlabel("Root Thickness / Chord")
plt.ylabel("Maximum Takeoff Weight (lbf)")
plt.title("t/c vs. Maximum Takeoff Weight")

plt.subplot(2, 1, 2)
plt.plot(t_c_root_vals, total_fuel_burn_vals / passengers)
plt.xlabel("Root Thickness / Chord")
plt.ylabel("Fuel Burn per Passenger Weight (lbf)")
plt.title("t/c vs. Fuel Burn per Passenger")

plt.show()

#Trade Study Three Wing Area Sweep
MTOW_vals = np.zeros(resolution)
MPOW_vals = np.zeros(resolution)
total_fuel_burn_vals = np.zeros(resolution)

wing_area_vals = np.linspace(600, 1000, resolution)
counter = 0 

for Wing_area in wing_area_vals:
    MTOW_new, MPOW, total_fuel_burn = tradeStudies(10.06, 0.15450, Wing_area, 350, 0.2, 0.2, 0.5, 0.5)       #Will Produce results from Weights Refinement
    MTOW_vals[counter] = MTOW_new
    MPOW_vals[counter] = MPOW
    total_fuel_burn_vals[counter] = total_fuel_burn
    counter = counter + 1
    print("Loop Count:", counter)

#Plotting
plt.figure(figsize=(8,12))
plt.subplot(2, 1, 1)
plt.plot(wing_area_vals, MTOW_vals)
plt.xlabel("Wing Area (ft^2)")
plt.ylabel("Maximum Takeoff Weight (lbf)")
plt.title("Wing Area vs. Maximum Takeoff Weight")

plt.subplot(2, 1, 2)
plt.plot(wing_area_vals, total_fuel_burn_vals / passengers)
plt.xlabel("Wing Area (ft^2)")
plt.ylabel("Fuel Burn per Passenger Weight (lbf)")
plt.title("Wing Area vs. Fuel Burn per Passenger")

plt.show()