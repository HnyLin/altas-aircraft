import Weight_Estimate_Calculator_A2_Henry as WE
import pandas as pd

#Concept 1 Inputs
passengers = 50  #Number of Passengers
mission = ("Start and Takeoff", "Climb", "Cruise", "Descent", "Loitter", "Landing")  # Mission Profile of Design Aircraft
hybridization_factors = (0.25, 0, 0, 0.5, 0.5, 0.5)
R = 1000         #nmi                        Range Per Cruise Interval
c = 0.4          #lb/(lbf*hr)                Fuel Burn
nu = 0.73        #Dimensionless              Battery efficency
Cl_Cd = 15.088   #Dimensionless              Lift/Drag
E = 0.025        #hours                      Lottier Endurance
V_inf = 350      #knots                      Cruise Velocity
e_star = 500*3600#wS/kg                      Battery Energy Density
A = 0.96         #EWF Coefficent
C = -0.05        #EWF Coefficent

W_crew_pass, W_payload = WE.prelim_pass_payload(passengers)
TOGW, EW, EWF, DTW, DTWF, stored_interval_fuel_fractions, interval_fuel_burn = WE.Calc_Estimated_Weight(passengers, mission, hybridization_factors, R, c, nu, Cl_Cd, E, V_inf, e_star, A, C)

print("For a Hybrid Electric Aircraft: ")
print(f"Crew and Passenger Weight (lbs): {W_crew_pass}")
print(f"Payload Weight (lbs): {W_payload}")
print("------------------------------------------------")
print("Takeoff Gross Weight Estimate (lbs): %0.2f" % TOGW)
print("Empty Weight (lbs): %0.2f" % EW)
print("Empty Weight Fraction %0.3f" % EWF)
print("Drivetrain Weight (lbs): %0.2f" % DTW )
print("Drivetrain Weight Fraction: %0.3f" % DTWF)
d = {'Mission Phase: ': mission, 'Hybridization Factor: ':hybridization_factors, 'Hybrid Fuel Burn (lbs): ': interval_fuel_burn}
fburntable = pd.DataFrame(data = d)
print(fburntable)
print("--------------------------------------------------------------------------------")
print("Hybrid Fuel Burn (lbs): ", round(sum(interval_fuel_burn)))

print("------------------------------------------------")
print("---------------------------------------------------------")
print("------------------------------------------------")

hybridization_factors = (0, 0, 0, 0, 0, 0)
TOGW, EW, EWF, DTW, DTWF, stored_interval_fuel_fractions, interval_fuel_burn = WE.Calc_Estimated_Weight(passengers, mission, hybridization_factors, R, c, nu, Cl_Cd, E, V_inf, e_star, A, C)
print("For a Gas Only Conventional Aircraft: ")
print(f"Crew and Passenger Weight (lbs): {W_crew_pass}")
print(f"Payload Weight (lbs): {W_payload}")
print("------------------------------------------------")
print("Takeoff Gross Weight Estimate (lbs): %0.2f" % TOGW)
print("Empty Weight (lbs): %0.2f" % EW)
print("Empty Weight Fraction %0.3f" % EWF)
print("Drivetrain Weight (lbs): %0.2f" % DTW )
print("Drivetrain Weight Fraction: %0.3f" % DTWF)
d = {'Mission Phase: ': mission, 'Hybridization Factor: ':hybridization_factors, 'Conventional Fuel Burn (lbs): ': interval_fuel_burn}
fburntable = pd.DataFrame(data = d)
print(fburntable)
print("--------------------------------------------------------------------------------")
print("Conventional Fuel Burn (lbs): ", round(sum(interval_fuel_burn)))

#For Dash 8 Q300
passengers = 50  #Number of Passengers
mission = ("Start and Takeoff", "Climb", "Cruise", "Descent", "Loitter", "Landing")  # Mission Profile of Design Aircraft
hybridization_factors = (0, 0, 0, 0, 0, 0)
R = 1000         #nmi                        Range Per Cruise Interval
c = 0.4          #lb/(lbf*hr)                Fuel Burn
nu = 0.73        #Dimensionless              Battery efficency
Cl_Cd = 15       #Dimensionless              Lift/Drag
E = 0.025        #hours                      Lottier Endurance
V_inf = 278      #knots                      Cruise Velocity
e_star = 500*3600#wS/kg                      Battery Energy Density
A = 0.96         #EWF Coefficent
C = -0.05        #EWF Coefficent

W_crew_pass, W_payload = WE.prelim_pass_payload(passengers)
TOGW, EW, EWF, DTW, DTWF, stored_interval_fuel_fractions, interval_fuel_burn = WE.Calc_Estimated_Weight(passengers, mission, hybridization_factors, R, c, nu, Cl_Cd, E, V_inf, e_star, A, C)

print("For a Hybrid Electric Aircraft: ")
print(f"Crew and Passenger Weight (lbs): {W_crew_pass}")
print(f"Payload Weight (lbs): {W_payload}")
print("------------------------------------------------")
print("Takeoff Gross Weight Estimate (lbs): %0.2f" % TOGW)
print("Empty Weight (lbs): %0.2f" % EW)
print("Empty Weight Fraction %0.3f" % EWF)
print("Drivetrain Weight (lbs): %0.2f" % DTW )
print("Drivetrain Weight Fraction: %0.3f" % DTWF)
d = {'Mission Phase: ': mission, 'Hybridization Factor: ':hybridization_factors, 'Hybrid Fuel Burn (lbs): ': interval_fuel_burn}
fburntable = pd.DataFrame(data = d)
print(fburntable)
print("--------------------------------------------------------------------------------")
print("Hybrid Fuel Burn (lbs): ", round(sum(interval_fuel_burn)))

print("------------------------------------------------")
print("---------------------------------------------------------")
print("------------------------------------------------")

hybridization_factors = (0, 0, 0, 0, 0, 0)
TOGW, EW, EWF, DTW, DTWF, stored_interval_fuel_fractions, interval_fuel_burn = WE.Calc_Estimated_Weight(passengers, mission, hybridization_factors, R, c, nu, Cl_Cd, E, V_inf, e_star, A, C)
print("For a Gas Only Conventional Aircraft: ")
print(f"Crew and Passenger Weight (lbs): {W_crew_pass}")
print(f"Payload Weight (lbs): {W_payload}")
print("------------------------------------------------")
print("Takeoff Gross Weight Estimate (lbs): %0.2f" % TOGW)
print("Empty Weight (lbs): %0.2f" % EW)
print("Empty Weight Fraction %0.3f" % EWF)
print("Drivetrain Weight (lbs): %0.2f" % DTW )
print("Drivetrain Weight Fraction: %0.3f" % DTWF)
d = {'Mission Phase: ': mission, 'Hybridization Factor: ':hybridization_factors, 'Conventional Fuel Burn (lbs): ': interval_fuel_burn}
fburntable = pd.DataFrame(data = d)
print(fburntable)
print("--------------------------------------------------------------------------------")
print("Conventional Fuel Burn (lbs): ", round(sum(interval_fuel_burn)))