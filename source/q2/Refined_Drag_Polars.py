'''
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

import numpy as np
import matplotlib.pyplot as plt

#Calculating Skin Friction Coefficent
def get_C_f(altitude, velocity, char_length, percent_lam_flow):

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

    Reynolds_component = rho * velocity * char_length / viscosity * 32.174

    C_f_laminar = 1.328 / np.sqrt(Reynolds_component)

    speed_of_sound = np.sqrt(R * T * gamma * 32.174)

    Mach_num = velocity / speed_of_sound

    C_f_turbulent = 0.455 / ( np.log10(Reynolds_component) ** 2.58 * ( 1 + 0.144 * Mach_num ** 2) ** 0.65 )

    C_f = C_f_laminar * (percent_lam_flow) + C_f_turbulent * ( 1 - percent_lam_flow)

    return C_f

#Calculating Zero Lift Drag
def get_CD_0(S_ref, drag_area_vals, skin_friction_coefficent_vals, form_factor_vals, interference_factor_vals, wetted_area_vals):

    #Miscellaneous Form Drag
    CD_miss = 1 / S_ref * sum(drag_area_vals)

    CD_0 = 1 / S_ref * sum(skin_friction_coefficent_vals * form_factor_vals * interference_factor_vals * wetted_area_vals) + CD_miss

    #Leak and Proturbance Drag (Est 5 - 10% of total parasite drag)
    CD_LP = 0.075 * CD_0

    CD_0 = CD_0 + CD_LP

    return CD_0

#Calculating Flap Drag
def get_flap_drag(flap_length, chord, flapped_area, S_ref, flap_angle):
    
    #For a plain and split flap
    delta_CD_flap = 1.7 * ( flap_length / chord ) ** 1.38 * (flapped_area / S_ref) * np.sin(flap_angle) ** 2

    return delta_CD_flap

#Calculating Trim Drag
def get_CD_trim(length_wingac_to_tailac, length_wingac_cg, CL_w, CM_ac_minus_t, tail_area, S_ref, mean_chord, AR_tail):

    V_HT = length_wingac_to_tailac * tail_area / ( S_ref * mean_chord )

    CL_t = ( CL_w * length_wingac_cg / mean_chord + CM_ac_minus_t ) * length_wingac_to_tailac / ( length_wingac_to_tailac - length_wingac_cg ) * 1 / V_HT

    oswald_eff = 1.78 * ( 1 - 0.045 * AR_tail ** 0.68 ) - 0.64 

    CD_trim = CL_t ** 2 / ( np.pi * oswald_eff * AR_tail ) * ( tail_area / S_ref )

    return CD_trim

#Induced Drag (From AVL)