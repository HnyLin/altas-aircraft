# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 14:34:17 2023

@author: Henry Lin
"""
"""
Description:
    To calculate the weight parameters for the designed aircraft you must provide:
        1) passengers (The Number of Passangers)
            Note: Assumes that for every 50 passengers there are 2 pilots and 1 cabin member
        2) mission (Supply a list of strings describing the phases of the mission)
            Possible Entries: Warmup, Taxi, Takeoff, Climb, Cruise, Loitter, Descent, Landing
        3) hybridization_factors (List of idential size to "mission" of values between 
                                  0 (pure gas power) and 1 (pure electric power))
            Note: The values of the mission stage and hybridization factor are indexed together, list must
            be in order of the mission
        4) R (Range in nautical miles (nmi))
        5) c (Fuel Burn of the gas turboprop (lbm / (lbf * hr)))
        6) nu (Battery Efficency (Dimensionless Quanitity between (0 and 1)))
        7) Cl_Cd (Lift / Drag Ratio)
        8) E (Endurance in hours (Used for loitter fuel weight ratio calculations))
        9) V_inf (Cruise speed in knots)
        10) e_star (Battery Energy Density (wS/kg))
        11) A (Empty Weight Factor Coefficent 1)
        12) C (Emply Weight Factor Coefficent 2)
            Note: Empty Weight Factor We / W0 = A * W0 ** C

Important: If the output is infinity it means that the weight of the battery exceeds 
the weight of the aircraft, therefore the plane is impossible.

Contact me if you have any questions -Henry                               
"""

import numpy as np
from scipy.optimize import least_squares
import pandas as pd

def prelim_pass_payload(passengers):
    
    #Calculating Crew and Crew Payload Weight lbs
    num_pilots = (int(passengers/50)+1)*2
    num_crew = int(passengers/50)+1
    W_crew = (num_crew+num_pilots)*190
    W_crew_bag = num_pilots*30
    
    #Calculating Passenger and Passenger Weight Weight lbs
    W_passengers = passengers*200
    W_passengers_bag = passengers*40
    
    W_crew_pass = W_crew + W_passengers
    W_payload = W_crew_bag + W_passengers_bag
    
    return (W_crew_pass, W_payload)

class FuelFactorCalculator:
    def __init__(self, R, c, nu, Cl_Cd, E, V_inf, e_star, mission, hybridization_factors):
        
        # hybridization factor of 1 = pure electric, 0 = pure gas
                
        self.mission = mission
        self.hybridization_factors = hybridization_factors
        
        self.gas_interval_fractions = {"Climb" : 0.985,
                                       "Descent": 0.990,
                                       "Cruise": np.exp(-R*c / (V_inf * Cl_Cd)),
                                       "Loitter": np.exp(-E*c / (Cl_Cd)),
                                       "Landing": 0.995,
                                       "Start and Takeoff": 0.970
                                       }
    
    def get_gas_fuel_fraction(self):
        
        stored_interval_fuel_fractions = []
        
        total_weight_fraction = 1
        for (hybridization, stage) in zip(self.hybridization_factors, self.mission):
            interval_fraction = self.gas_interval_fractions[stage]
            hybrid_interval_fraction = 1 - (1 - interval_fraction) * (1 - hybridization)
            
            stored_interval_fuel_fractions.append(hybrid_interval_fraction)
            
            total_weight_fraction *= hybrid_interval_fraction

        fuel_fraction = 1 - total_weight_fraction
        fuel_fraction = fuel_fraction * 1.06     #6% reserves and trapped fuel

        return fuel_fraction, stored_interval_fuel_fractions
    
    def get_electric_fuel_fraction(self):
        
        total_weight_fraction = 1
        for stage in self.mission:
            interval_fraction = self.gas_interval_fractions[stage]
            total_weight_fraction *= interval_fraction
         
        #Calculating the weight of cruise fuel
        temp_ff = 1 - total_weight_fraction
        temp_cruise_interval_ff = self.gas_interval_fractions["Cruise"]
        cruise_weight = temp_ff * temp_cruise_interval_ff
        
        
        temp_weight_fraction = np.zeros(len(self.mission))
        temp_weight_fraction[0] = 1
        
        mission_gas_interval_fraction = np.array([self.gas_interval_fractions[segment] for segment in self.mission])
        
        
        for i in range(1, len(self.mission)):
            temp_weight_fraction[i] = temp_weight_fraction[i-1] * mission_gas_interval_fraction[i]
        
        temp_weight_fraction = temp_weight_fraction / cruise_weight
        
        electric_fuel_fraction = (1 - mission_gas_interval_fraction) * temp_weight_fraction
        
        cruise_segments = [segment == "Cruise" for segment in self.mission]
        
        electric_fuel_fraction[cruise_segments] = 1
        
        energy_fraction = np.sum(np.array(self.hybridization_factors) * electric_fuel_fraction)
        
        return energy_fraction
    

def Calc_Estimated_Weight(passengers, mission, hybridization_factors, R, c, nu, Cl_Cd, E, V_inf, e_star, A, C):
 
    #Count the number of crusise intervals
    count = pd.Series(mission).value_counts()
    num_cruise = count._get_value("Cruise")
    R = R / num_cruise
    
    #Solve For Passenger/Crew/Baggage Weight    
    
    W_crew_pass, W_payload = prelim_pass_payload(passengers)
    
    #Solve for battery mass
    m_battery = R * 1852 / (nu * e_star * Cl_Cd)    
    
    #Iterative Solver
    
    W_guess = 10000             #lbs
    g = 9.81                    #m/s
    
    #Fuel Factor Calculation
    
    f_factors = FuelFactorCalculator(R, c, nu, Cl_Cd, E, V_inf, e_star, mission, hybridization_factors)
    
    fuel_fraction, stored_interval_fuel_fractions = f_factors.get_gas_fuel_fraction()
    
    def minimize_function(W0):
        
        gas_FF, stored_interval_fuel_fractions = f_factors.get_gas_fuel_fraction()
        
        empty_weight_fraction = A * W_guess**C
        
        return np.abs(W0 - (W_crew_pass + W_payload) / (1 - empty_weight_fraction - gas_FF 
         - (m_battery * g) * f_factors.get_electric_fuel_fraction()))
    
    solution = least_squares(minimize_function, x0 = W_guess, bounds = (0, np.inf))
    solution.x = np.inf if solution.x < 100 else solution.x
    
    EWF = A * solution.x ** C
    EW = EWF * solution.x
    DTW = solution.x - W_crew_pass - W_payload - EW
    DTWF = DTW / solution.x   
    TOGW = solution.x
    
    CW = TOGW
    current_weight = []
    for i in range(len(stored_interval_fuel_fractions)):
        CW_new = CW * stored_interval_fuel_fractions[i]
        CW = CW_new
        current_weight.append(CW_new)

    current_weight.insert(0, TOGW) 
    current_weight = np.array(current_weight) 
    interval_fuel_burn = np.empty(len(current_weight) - 1)
   
    for i in range(len(interval_fuel_burn)):
        interval_fuel_burn[i] = current_weight[i] - current_weight[i+1]
    
    return TOGW, EW, EWF, DTW, DTWF, stored_interval_fuel_fractions, interval_fuel_burn
