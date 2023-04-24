'''
Account for:
Zero Lift Drag
Drag due to Flaps
Lift Induced Drag
Trim Drag

Required Drag Curves
Clean
Takeoff Flaps, Gear Up
Takeoff Flaps, Gear Down
Landing Flaps, Gear Up
Landing Flaps, gear Down
'''

import numpy as np

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
