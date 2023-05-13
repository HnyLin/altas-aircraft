import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import time
import pandas as pd

#Functions

#Optimization Code / Weight Estimate / Fuel Burn Estimate
#================================================================================================================
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

#================================================================================================================

#Refined Drag Polars
'''
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Requires Induced_Drag_Data.xlsx (AVL Data) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Account for:
Zero Lift Drag
Drag due to Flaps
Trim Drag
Lift Induced Drag

Required Drag Curves
Clean
Takeoff Flaps, Gear Up
Takeoff Flaps, Gear Down
Landing Flaps, Gear Up
Landing Flaps, gear Down
'''
def getMach(altitude, velocity):

    #Fixed Values
    gamma = 1.4
    R = 53.35                                                   #ft*lbf/(lbm * R)

    #Atmophereic Data
    t_interp = np.array([59, 23.36, -12.26, -47.83, -69.70])
    t_interp = t_interp + 459.67                                #Rankine Conversion
    h_interp = [0, 10000, 20000, 30000, 40000]                  #ft

    T = np.interp(altitude, h_interp, t_interp)

    speed_of_sound = np.sqrt(R * T * gamma * 32.174)

    print("Speed of Sound: ", speed_of_sound)

    Mach_num = velocity / speed_of_sound

    return Mach_num

#Calculating Skin Friction Coefficent
def get_C_f(altitude, velocity, char_length_vals, percent_lam_flow_vals):

    #Fixed Values
    gamma = 1.4
    R = 53.35               #ft*lbf/(lbm * R)

    #Interpolation Data Bank
    rho_interp = [0.0765, 0.0565, 0.0408, 0.0287, 0.0189]       #lbm/ft^3
    t_interp = np.array([59, 23.36, -12.26, -47.83, -69.70])
    t_interp = t_interp + 459.67                                #Rankine Conversion
    mu_interp = np.array([3.737, 3.534, 3.324, 3.107, 2.969])
    mu_interp = mu_interp * 10 ** (-7)                          #slug / (ft*s)

    h_interp = [0, 10000, 20000, 30000, 40000]

    #Interpolates Enviromental Values from Altitude
    rho = np.interp(altitude, h_interp, rho_interp)
    T = np.interp(altitude, h_interp, t_interp)
    viscosity = np.interp(altitude, h_interp, mu_interp)

    Reynolds_component = rho * velocity * char_length_vals / viscosity * 32.174

    C_f_laminar_vals = 1.328 / np.sqrt(Reynolds_component)

    speed_of_sound = np.sqrt(R * T * gamma * 32.174)

    Mach_num = velocity / speed_of_sound

    C_f_turbulent_vals = 0.455 / ( np.log10(Reynolds_component) ** 2.58 * ( 1 + 0.144 * Mach_num ** 2) ** 0.65 )

    C_f_vals = C_f_laminar_vals * (percent_lam_flow_vals) + C_f_turbulent_vals * ( 1 - percent_lam_flow_vals)

    return C_f_vals

#Calculating Zero Lift Drag
def get_CD_0(S_ref, drag_area_vals, skin_friction_coefficent_vals, form_factor_vals, interference_factor_vals, wetted_area_vals):

    #Miscellaneous Form Drag
    CD_miss = 1 / S_ref * np.sum(drag_area_vals)

    CD_0 = 1 / S_ref * np.sum(skin_friction_coefficent_vals * form_factor_vals * interference_factor_vals * wetted_area_vals) + CD_miss

    #Leak and Proturbance Drag (Est 5 - 10% of total parasite drag)
    CD_LP = 0.075 * CD_0

    CD_0 = CD_0 + CD_LP

    return CD_0

#Calculating Flap Drag
def get_flap_drag(flap_length, chord, flapped_area, S_ref, flap_angle, slat_angle, slat_length, slatted_area):

    #Flap Angle Degrees to Rad
    flap_angle = flap_angle * np.pi / 180
    slat_angle = slat_angle * np.pi / 180
    
    #For slotted flaps
    delta_CD_flap = 0.9 * ( flap_length / chord ) ** 1.38 * (flapped_area / S_ref) * np.sin(flap_angle) ** 2

    #For slotted slats
    delta_CD_slat = 0.9 * (slat_length / chord) ** 1.38 * (slatted_area / S_ref) * np.sin(slat_angle) ** 2

    delta_CD_flap_slat = delta_CD_flap + delta_CD_slat

    return delta_CD_flap_slat

#Calculating Trim Drag
def get_CD_trim(length_wingac_to_tailac, length_wingac_cg, CL_w, CM_ac_minus_t, tail_area, S_ref, mean_chord, AR_tail):

    V_HT = length_wingac_to_tailac * tail_area / ( S_ref * mean_chord )

    CL_t = ( CL_w * length_wingac_cg / mean_chord + CM_ac_minus_t ) * length_wingac_to_tailac / ( length_wingac_to_tailac - length_wingac_cg ) * 1 / V_HT

    oswald_eff = 1.78 * ( 1 - 0.045 * AR_tail ** 0.68 ) - 0.64 

    CD_trim = CL_t ** 2 / ( np.pi * oswald_eff * AR_tail ) * ( tail_area / S_ref )

    return CD_trim

#Induced Drag (From AVL)
#Landing
df = pd.read_excel(r'C:\Users\henry\OneDrive\Documents\EAE130B\atlas-aircraft\source\q2\Induced_Drag_Data.xlsx', sheet_name='Landing', )
CD_i_landing_vals = df['CD_i']
Cl_max_landing_vals = df['Clmax']

CD_i_landing_vals = CD_i_landing_vals.to_numpy()
Cl_max_landing_vals = Cl_max_landing_vals.to_numpy()

#Takeoff
df = pd.read_excel(r'source\q2\Induced_Drag_Data.xlsx', sheet_name='Takeoff', )
CD_i_takeoff_vals = df['CD_i']
Cl_max_takeoff_vals = df['Clmax']

CD_i_takeoff_vals = CD_i_takeoff_vals.to_numpy()
Cl_max_takeoff_vals = Cl_max_takeoff_vals.to_numpy()

#Clean
df = pd.read_excel(r'C:\Users\henry\OneDrive\Documents\EAE130B\atlas-aircraft\source\q2\Induced_Drag_Data.xlsx', sheet_name='Clean', )
CD_i_clean_vals = df['CD_i']
Cl_max_clean_vals = df['Clmax']

CD_i_clean_vals = CD_i_clean_vals.to_numpy()
Cl_max_clean_vals = Cl_max_clean_vals.to_numpy()

#================================================================================================================

#================================================================================================================

#================================================================================================================

print("========================================================================")
print("Optimization")
print("========================================================================")

#Calculating Fuel Burn Per Passenger for Dash 8-q300
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
    Fuel_Fraction_Calculator(AR, Wing_area, c_f, c, d, MTOW, MPOW, SFC, R, segments, eta, h_cruise, V_cruise, hybridization_factors)

D8fuel_burn_per_pass = D8total_fuel_burn / 50
print("Dash 8-q300 Fuel Weight Per Passenger 500 nmi range (lbf): ", round(D8fuel_burn_per_pass, 2))

#================================================================================================================


#For Calulating Optimium Aircraft Parameters (Commented Out Due to Long Run Time)
'''
#Setting Initial Guess
initial_guess = [12.06, 0.15450, 800, 350, 0.25, 0.25, 0.25, 0.25]

#Setting Bounds
bound_vals = ((10, 13.14), (0.1, 0.25), (600, 1000), (280, 450), (0, 1), (0, 1), (0, 1), (0, 1))

#Optimize
start_time = time.time()
result = optimize.minimize(objective_function, x0 = initial_guess, bounds = bound_vals, options= {'disp': True}, tol = 10 ** -8 )
end_time = time.time()
print("Elapsed Timed (min): ", (end_time - start_time)/60)
print("Optimized Values")
print(result.x)

print("Optimum Fuel Weight (lbf): ", result.fun)
print("Optimum Fuel Weight Per Passenger (lbf): ", result.fun/50)

'''
#Optimization Results
optimized_fuel_weight = 3465                #lbf
MTOW = 40265                                #lbf
AR = 12.19
t_c_root = 0.25
S_ref = 800                                 #ft^2
V_cruise = 350                              #ktas
h1 = 0.24                                   #Start Warmup Taxi
h2 = 0.16                                   #Takeoff
h3 = 0.37                                   #Descent
h4 = 0.37                                   #Landing

optimized_fuel_weight_per_pass = optimized_fuel_weight / 50
print("Optimized Fuel Weight Per Passenger 500 nmi Range (lbf): ", round(optimized_fuel_weight_per_pass, 2))

#Calculating Percent Difference (Block Fuel 500 nmi)
percent_difference = (optimized_fuel_weight - D8total_fuel_burn) / abs( D8total_fuel_burn ) * 100
print("Percent Difference in Fuel Burn: ", round(percent_difference, 1))
#================================================================================================================

#Refined Drag Polars
print("========================================================================")
print("Refined Drag Polars")
print("========================================================================")
Cl_max = 3.3
rho_takeoff_landing = 0.0659                #lbm/ft^3

V_stall = np.sqrt( ( 2 * MTOW * 32.174 ) / ( Cl_max * rho_takeoff_landing * S_ref) )

print("Stall Velcocity (ft/s): ", V_stall)

V_takeoff_landing = 1.3 * V_stall   #ft/s
h_takeoff_landing = 5000            #ft
h_cruise = 28000                    #ft

Mach_takeoff_landing = getMach(5000, V_takeoff_landing)
Mach_cruise = getMach(28000, V_cruise)

'''
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Revise Component Data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'''

#Component Data
#[Wing Section 1, Wing Section 2, V Tail, H Tail, Winglet, Nacelle, Fueselage]

char_length_vals = np.array([11.904, 7.0956, 18.32, 5.584, 2.3196, 14.167, 81])

percent_lam_flow_vals = np.array( [0.5, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25] )

wetted_area_vals = np.array( [883.718, 839.026, 458.405, 389.813, 46.496, 52.454, 2093] )      

interference_factor_vals = np.array( [1, 1, 1, 1, 1, 1.5, 1] )               

form_factor_vals_takeoff_landing = np.array( [1.33018616, 1.325467794, 1.261023316, 1.289911252, 1.259013234, 1.058329216, 1.133150585] )

form_factor_vals_cruise = np.array( [1.425822647, 1.419737634, 1.33662719, 1.373882348, 1.3340349, 1.058329216, 1.133150585] )

drag_area_vals_geardown = np.array([(0.139+0.419*(Mach_takeoff_landing - 0.161)**2), 0.15, 0.15, 0.25]) * 91

drag_area_vals_gearup = np.array([(0.139+0.419*(Mach_takeoff_landing - 0.161)**2) * 91])

#Calculating Coefifecent of Friction Values

altitude = h_takeoff_landing
velocity = V_takeoff_landing
Cf_vals_takeoff_landing = get_C_f(altitude, velocity, char_length_vals, percent_lam_flow_vals)
# print("Cf Takeoff Landing:", Cf_vals_takeoff_landing)

altitude = h_cruise
velocity = V_cruise
Cf_vals_cruise = get_C_f(altitude, velocity, char_length_vals, percent_lam_flow_vals)
# print("Cf Cruise: ", Cf_vals_cruise)

#Zero Lift Drag (Takeoff Landing (Gearup) )
skin_friction_coefficent_vals = Cf_vals_takeoff_landing
form_factor_vals = form_factor_vals_takeoff_landing
drag_area_vals = drag_area_vals_gearup
CD_0_takeoff_landing_gearup = get_CD_0(S_ref, drag_area_vals, skin_friction_coefficent_vals, form_factor_vals, interference_factor_vals, wetted_area_vals)
# print("CD_0 Landing Takeoff (Gearup): ", CD_0_takeoff_landing_gearup)

#Zero Lift Drag (Takeoff Landing (Gear Down) )
skin_friction_coefficent_vals = Cf_vals_takeoff_landing
form_factor_vals = form_factor_vals_takeoff_landing
drag_area_vals = drag_area_vals_geardown
CD_0_takeoff_landing_geardown = get_CD_0(S_ref, drag_area_vals, skin_friction_coefficent_vals, form_factor_vals, interference_factor_vals, wetted_area_vals)
# print("CD_0 Landing Takeoff (Gear Down): ", CD_0_takeoff_landing_geardown)

#Zero Lift Drag (Cruise)
skin_friction_coefficent_vals = Cf_vals_cruise
form_factor_vals = form_factor_vals_cruise
drag_area_vals = drag_area_vals_gearup
CD_0_cruise = get_CD_0(S_ref, drag_area_vals, skin_friction_coefficent_vals, form_factor_vals, interference_factor_vals, wetted_area_vals)
# print("CD_0 Cruise: ", CD_0_cruise)

#Calculating Flap Drag (Takeoff)
flap_angle_takeoff = 30             #degrees
flap_length = (5.187 + 4.016) /2    #ft
chord = (12.829 + 10.997) /2        #ft
flapped_area = 122.73               #ft^2
flap_angle = flap_angle_takeoff
slat_angle = 0                          #No Slats
slat_length = 0                         #No Slats
slatted_area = 0                        #No Slats
delta_CD_flap_slat_takeoff = get_flap_drag(flap_length, chord, flapped_area, S_ref, flap_angle, slat_angle, slat_length, slatted_area)

# print("Delta C_D Flaps and Slats (Takeoff): ", delta_CD_flap_slat_takeoff)

#Calculating Flap Drag (Landing)
flap_angle_landing = 70     #degrees
delta_CD_flap_slat_landing = get_flap_drag(flap_length, chord, flapped_area, S_ref, flap_angle, slat_angle, slat_length, slatted_area)

# print("Delta C_D Flaps and Slats (Landing): ", delta_CD_flap_slat_landing)

#Clean (Cruise)
CL_clean = Cl_max_clean_vals
CD_clean = CD_i_clean_vals + CD_0_cruise

#Takeoff Flaps, Gear Up
CL_takeoff_gearup = Cl_max_takeoff_vals
CD_takeoff_gearup = CD_i_takeoff_vals + CD_0_takeoff_landing_gearup + delta_CD_flap_slat_takeoff

#Takeoff Flaps, Gear Down
CL_takeoff_geardown = Cl_max_takeoff_vals
CD_takeoff_geardown = CD_i_takeoff_vals + CD_0_takeoff_landing_geardown + delta_CD_flap_slat_takeoff

#Landing Flaps, Gear Up
CL_landing_gearup = Cl_max_landing_vals
CD_landing_gearup = CD_i_landing_vals + CD_0_takeoff_landing_gearup + delta_CD_flap_slat_landing

#Landing Flaps, gear Down
CL_landing_geardown = Cl_max_landing_vals
CD_landing_geardown = CD_i_landing_vals + CD_0_takeoff_landing_geardown + delta_CD_flap_slat_landing


#Plotting
#Takeoff Gear Up
plt.figure(figsize=(12, 12))
plt.plot(CD_takeoff_gearup, CL_takeoff_gearup, label = "Takeoff Flaps, Gear Up", marker = ".", markersize = 10)
plt.plot(CD_takeoff_geardown, CL_takeoff_geardown, label = "Takeoff Flaps, Gear Down", marker = ".", markersize = 10)
plt.plot(CD_landing_gearup, CL_landing_gearup, label = "Landing Flaps, Gear Up", marker = ".", markersize = 5)
plt.plot(CD_landing_geardown, CL_landing_geardown, label = "Landing Flaps, Gear Down", marker = ".", markersize = 5)
plt.plot(CD_clean, CL_clean, label = "Clean", marker = ".", markersize = 10)
plt.ylabel("$C_L$")
plt.xlabel("$C_D$")
plt.title("$C_L$ vs $C_D$ Drag Polars")
plt.legend()
plt.show()

#================================================================================================================
