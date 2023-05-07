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
import pandas as pd

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




#Aircraft Geometry

#Gneral Parameters
MTOW = 67551                        #lbf
S_ref = 805.06                      #ft^2

V_cruise = 350 * 1.688      #ft/s
V_stall = 121                       #ft/s
V_takeoff_landing = 1.3 * V_stall   #ft/s
h_takeoff_landing = 5000            #ft
h_cruise = 28000                    #ft

Mach_takeoff_landing = getMach(5000, V_takeoff_landing)
print("Mach Takeoff Landing: ", Mach_takeoff_landing)

Mach_cruise = getMach(28000, V_cruise)
print('Mach Cruise: ', Mach_cruise)

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
print("Cf Takeoff Landing:", Cf_vals_takeoff_landing)

altitude = h_cruise
velocity = V_cruise
Cf_vals_cruise = get_C_f(altitude, velocity, char_length_vals, percent_lam_flow_vals)
print("Cf Cruise: ", Cf_vals_cruise)

#Zero Lift Drag (Takeoff Landing (Gearup) )
skin_friction_coefficent_vals = Cf_vals_takeoff_landing
form_factor_vals = form_factor_vals_takeoff_landing
drag_area_vals = drag_area_vals_gearup
CD_0_takeoff_landing_gearup = get_CD_0(S_ref, drag_area_vals, skin_friction_coefficent_vals, form_factor_vals, interference_factor_vals, wetted_area_vals)
print("CD_0 Landing Takeoff (Gearup): ", CD_0_takeoff_landing_gearup)

#Zero Lift Drag (Takeoff Landing (Gear Down) )
skin_friction_coefficent_vals = Cf_vals_takeoff_landing
form_factor_vals = form_factor_vals_takeoff_landing
drag_area_vals = drag_area_vals_geardown
CD_0_takeoff_landing_geardown = get_CD_0(S_ref, drag_area_vals, skin_friction_coefficent_vals, form_factor_vals, interference_factor_vals, wetted_area_vals)
print("CD_0 Landing Takeoff (Gear Down): ", CD_0_takeoff_landing_geardown)

#Zero Lift Drag (Cruise)
skin_friction_coefficent_vals = Cf_vals_cruise
form_factor_vals = form_factor_vals_cruise
drag_area_vals = drag_area_vals_gearup
CD_0_cruise = get_CD_0(S_ref, drag_area_vals, skin_friction_coefficent_vals, form_factor_vals, interference_factor_vals, wetted_area_vals)
print("CD_0 Cruise: ", CD_0_cruise)

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

print("Delta C_D Flaps and Slats (Takeoff): ", delta_CD_flap_slat_takeoff)

#Calculating Flap Drag (Landing)
flap_angle_landing = 70     #degrees
delta_CD_flap_slat_landing = get_flap_drag(flap_length, chord, flapped_area, S_ref, flap_angle, slat_angle, slat_length, slatted_area)

print("Delta C_D Flaps and Slats (Landing): ", delta_CD_flap_slat_landing)

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
plt.plot(CD_takeoff_gearup, CL_takeoff_gearup)
plt.ylabel("$C_L$")
plt.xlabel("$C_D$")
plt.title("$C_L$ vs $C_D$ Takeoff Gear Up")
plt.show()

#Takeoff Gear Down
plt.figure(figsize=(12, 12))
plt.plot(CD_takeoff_geardown, CL_takeoff_geardown)
plt.ylabel("$C_L$")
plt.xlabel("$C_D$")
plt.title("$C_L$ vs $C_D$ Takeoff Gear Down")
plt.show()

#Landing Gear Up
plt.figure(figsize=(12, 12))
plt.plot(CD_landing_gearup, CL_landing_gearup)
plt.ylabel("$C_L$")
plt.xlabel("$C_D$")
plt.title("$C_L$ vs $C_D$ Landing Gear Up")
plt.show()

#Landing Gear Down
plt.figure(figsize=(12, 12))
plt.plot(CD_landing_geardown, CL_landing_geardown)
plt.ylabel("$C_L$")
plt.xlabel("$C_D$")
plt.title("$C_L$ vs $C_D$ Landing Gear Down")
plt.show()

#Clean
plt.figure(figsize=(12, 12))
plt.plot(CD_clean, CL_clean)
plt.ylabel("$C_L$")
plt.xlabel("$C_D$")
plt.title("$C_L$ vs $C_D$ Cruise Clean")
plt.show()