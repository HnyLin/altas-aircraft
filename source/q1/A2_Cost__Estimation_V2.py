# Cost Estimation V2
'''
First: This is post-A2 and will only consider the cost of our selected aircraft, AT-V1 (the first design)

Second: All equations/Factors (Variables that start with "F_") come from Finger et al. source, Section 2 [Finger S.2] unless otherwise commented
(or are determined by our design).Many of these equations use the same factors, however the value of those factors change depending on the cost 
one is examining. Because of this, any variable that changes value in subsequent equations will be assigned that new value before that next equation.
This may make the code look messy... oh well

Third: Equations are modified slighty, namely they do not consider F_cert in the original equations. The choice for this is due to a lack of
knowledge in the certification process, as well as we just haven't covered the certifications required (besides FAR 25) in the RFP. Additionally, Prof. Harvey believes
it won't be relevant for this case.

Fourth: This uses W_TO, W_Bat, and P_SHP values that were from initial estimations. W_TO, W_Bat, and P_SHP (P_TO) are calculated among others in A4 as dimensional values,
and should be used if possible. To highlight where these substitutes were used, look for "************************************** " in comments
'''
#Library 
from math import log10, floor
#Function
def round_sig(x, sig=3):
    '''Simple Significant figure calculator, rounds to 3 significant figures'''
    return round(x, sig-int(floor(log10(abs(x))))-1)

#==========================================================================================================
#RTD&E Costs: Engineering, Tooling,  Manufacturing, Development, Flight Testing, Quality Control, Materials
#==========================================================================================================
#Variables
y = 2022        #Year of for analysis
Q = 50          #Number of aircraft produced in 5-year period 
Q_M = 50/60     #"                 " produced in one month (assume consisent buillding over 5 year aircraft)
V_H = 275       #Max level airspeed, KTAS (FOR NOW put cruise speed)
Q_Proto = 2     #Number of prototypes (assume 2)
CPI = 1.24037   #CPI  from 2012 to 2022 (using CPI caluclator for $1000 Jan'12 to $1000 Jan'22 [https://www.bls.gov/data/inflation_calculator.htm])
                #Using this because Finger et al. equations  use CPI from 2012 to present day

W_TO = 82560        #**************************************
W_airframe = 45000  #**************************************

#W_TO = 75000    #testing
#W_airframe = 32500


f = 0.4     #Experimenting with Composites Factor, see lines below for explanation

#CoE (Cost of ... Engineering)
F_CF = 1.03
F_Comp = f*2.00
F_Press = 1.03
F_HyE = 1.66
R_Eng = 2.576*y - 5058      #[Nikolai Fig 24.4]
C_Eng = 0.083*W_airframe**0.791*V_H**1.521*Q**0.183*F_CF*F_Comp*F_Press*F_HyE*R_Eng*CPI

#CoT
F_Taper = 1                 #we have a tapered aircraft
F_CF = 1.02
F_Comp = f*2.00           #Assumes 40% composites in the aircraft (multplying the 100% composites factor by 40%)
F_Press = 1.01
F_HyE = 1.10
R_Tool = 2.883*y - 5666     #[Nikolai Fig 24.4]
C_Tool = 2.1036*W_airframe**0.764*V_H**0.899*Q**0.178*Q_M**0.066*F_Taper*F_CF*F_Comp*F_Press*F_HyE*R_Tool*CPI

#CoMFG
F_CF = 1.01
F_Comp = f*1.25           #Assumes 40% composites in the aircraft
F_HyE = 1.10
R_MFG =2.316*y - 4552       #[Nikolai Fig 24.4]
C_MFG = 20.2588*W_airframe**0.74*V_H**0.543*Q**0.524*F_CF*F_Comp*F_HyE*R_MFG*CPI

#CoD
F_CF = 1.01
F_Comp = f*1.50           #Assumes 40% composites in the aircraft
F_Press = 1.03
F_HyE = 1.05
C_Dev = 0.06458*W_airframe**0.873*V_H**1.89*Q_Proto**0.346*F_CF*F_Comp*F_Press*F_HyE*CPI

#CoFT
F_HyE = 1.50
C_FT = 0.009646*W_airframe**1.16*V_H**1.3718*Q_Proto**1.281*F_HyE*CPI

#CoQC
F_Comp = f*1.50           #Assumes 40% composites in the aircraft
F_HyE = 1.50
C_QC = 0.13*C_MFG*F_Comp*F_HyE

#CoM
F_CF = 1.02
F_Press = 1.01
F_Hye = 1.05
C_Mat = 24.896*W_airframe**0.689*V_H**0.624*Q**0.792*F_CF*F_Press*F_HyE*CPI
#=========================================================================================
#Propulsion Costs: Int. Comb. Engine, Elect. Motor, Power Mngmt. Sys., Battery, Propellers
#=========================================================================================
#Variables
N_engine = 2        #Number of Engine
N_Motor = 2
N_Prop = 2
D_P = 11            #11ft Diameter (from OpenVSP model)



W_bat = W_TO*0.3        #************************************** Should be from A4 Dimensional Values, just estimating now
E_Bat = 0.5*W_bat/2.2   #We assume a battery specific energy of 0.5kWhr/kg, and our battery weight is in lbm

#Calculating Power of motors/comb engine from 
H_P = 0.25
P_SHP = 7338        #**************************************Total Power needed (takeoff) (P_EM + P_ICE)
#P_SHP = 4760            #testing
P_EM = 0.25*P_SHP
P_ICE = (1 - H_P)*P_SHP

P_EM_tot = N_Motor*P_EM

#CoICE (Cost of ...Int. Comb. Engine)
C_ICE = 174*N_engine*P_ICE*CPI

#CoEM
C_EM = 174*N_Motor*P_EM*CPI

#CoPMS
C_PMS = 150*P_EM_tot*CPI

#CoB
C_Bat = 200*E_Bat*CPI

#CoP
C_CSProp = 210*N_Prop*D_P**2*(P_SHP/D_P)**0.12*CPI
#================================
#Sum Everything Up
#===============================
#Unit Cost?
C_TOT = C_CSProp + C_Bat + C_PMS + C_EM + C_ICE + C_Mat + C_QC + C_FT + C_Dev + C_MFG + C_Tool + C_Eng
C_unit_price = C_TOT/Q
#Sales Price
Profit = C_unit_price*1.15 - C_unit_price #15% Profit Margin
C_unit_sales_price = C_unit_price*1.15 

print("Cost of Engineering: %0.3f"% C_Eng)
print("Cost of Tooling: %0.3f"% C_Tool)
print("Cost of Manufacturing: %0.3f"% C_MFG)
print("Cost of Development: %0.3f"% C_Dev)
print("Cost of Flight Testing: %0.3f"% C_FT)
print("Cost of Quality Control: %0.3f"% C_QC)
print("Cost of Materials: %0.3f"% C_Mat)
print("Cost of Combustion Engine: %0.3f"% C_ICE)
print("Cost of Electric Motors: %0.3f"% C_EM)
print("Cost of Power Management System: %0.3f"% C_PMS)
print("Cost of Batteries: %0.3f"% C_Bat)
print("Cost of Propeller: %0.3f"% C_CSProp)
print("===================================")
print("Total Program Cost (USD): %0.3f"% C_TOT)
print("Unit Price (USD): %0.3f"% C_unit_price)
print("===================================")
print("Unit Price (USD): %0.3f"% Profit)
print("Unit Sales Price (USD): %0.3f"% C_unit_sales_price)
print('\n')
print("===================================")
print('PRICE ROUNDED TO 3 SIGNIFICANT FIGURES')
print("===================================")
print('\n')

#Sig Figs: 3 (based off of factors used and  max level flight 275 from RFP and )

print("Cost of Engineering: %0.3f"% round_sig(C_Eng))
print("Cost of Tooling: %0.3f"% round_sig(C_Tool))
print("Cost of Manufacturing: %0.3f"% round_sig(C_MFG))
print("Cost of Development: %0.3f"% round_sig(C_Dev))
print("Cost of Flight Testing: %0.3f"% round_sig(C_FT))
print("Cost of Quality Control: %0.3f"% round_sig(C_QC))
print("Cost of Materials: %0.3f"% round_sig(C_Mat))
print("Cost of Combustion Engine: %0.3f"% round_sig(C_ICE))
print("Cost of Electric Motors: %0.3f"% round_sig(C_EM))
print("Cost of Power Management System: %0.3f"% round_sig(C_PMS))
print("Cost of Batteries: %0.3f"% round_sig(C_Bat))
print("Cost of Propeller: %0.3f"% round_sig(C_CSProp))
print("===================================")
print("Total Program Cost (USD): %0.3f"% round_sig(C_TOT))
print("Unit Price (USD): %0.3f"% round_sig(C_unit_price))
print("===================================")
print("Profit Margin (USD): %0.3f"% round_sig(Profit))
print("Unit Sales Price (USD): %0.3f"% round_sig(C_unit_sales_price))

print(R_Eng, R_MFG,R_Tool)