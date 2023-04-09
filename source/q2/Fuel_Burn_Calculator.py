#By Henry Lin
#4/8/2023
#Fuel Burn Calculator
import numpy as np
from scipy.optimize import least_squares
import pandas as pd

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