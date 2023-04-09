#Libraries
import numpy as np
#======================================================================
#STEP 1 Read Inputs
W_P = 11.25     #lbf/hp
W_S = 78        #lbf/sqft
W_TO = 82561.1  #ATV1 GTOW from A2, lbf
Phi = 0.25      #Hybridization Ratio assumed constant for now (for simplicity, constant for entire flight)

#Otras variables
g = 32.17               #gravity ft/s^2
R = 1000*1.15078*5280   #range, ft (from nmi)
A = 0.96                #EWF Coefficent
C = -0.05               #EWF Coefficent
k_pg = 1.35             #factor for multi-engine prop, from Metabook 7.3.3
P_rshp = 2380           #Rated SHP, (from Dash8-300 .xslx)
W_crew = 11140          #weight of crew and passengers, lbf (from A2 weight estimate)
W_payload = 2120        #weight of payload, lbf (from A2 weight estimate)
wingWeight = 10         #typical weight per unit area for transport aircraft, lbf/sqft


L_D = 15.09                     #Lift to Drag ratio of ATV1
W_pt = 500                      #MADEUP, need research from comparable aircraft engine architecture
e_bat = 500*5.604423234e-5                #Battery specific energy (500 J/kg-hr to BTU/slug-s [convert to lbf*ft/s later in code] (Note: originally had 24923*778.2, may be incorrect)
e_fuel = 43.1e6*5.604423234e-5     #Jet-A1 specific energy (from 43.1MJ/kg to BTU/slug-s    (Note: originally had 1.86e4*32.3*778.2, may be incorrect)
                                                    #(source: https://www.skytanking.com/news-info/glossary/jet-fuel/#:~:text=The%20specific%20energy%20of%20the,and%20lots%20of%20hot%20air.)
                                                    #Assume can provide this energy per second (thus original units are MJ/(kg*s)

eta_1p = 0.3*0.96           #overall fuel efficiency for parallel hybrid
eta_2p = 0.99*0.96*0.96     #overall battery efficiency for parallel hybrid
eta_3p = 0.9                #"       " assuming propeller efficiency 0.90 (from selected component efficiecies in Table 4.1 of RFP)



print('Battery energy is (hp/slug*s)')
print(e_bat)
#======================================================================
#Create functions for loop (Step 2)
def calcWPS(W_TO, W_P, W_S):
    '''
    DESCRIPTION:
    Simple calculator for calculating new power and wing area from power loading, wing loading and weight 
    
    INPUT:
        W       -   Weight, lbs
        W_P     -   Power Loading, lbs/hp
        W_S     -   Wing Loading, lbs/sqft
    

    OUTPUT:
        P_new   -   Power, lbf*ft/s
        S_new   -   Wing Area, sqft

    '''
    P_new = W_TO/W_P*550        #hp to lbf*ft/s
    S_new = W_TO/W_S            #sqft
    return P_new, S_new

#======================================================================
#Create functions for loop (Step 3.1) changing major stuff

def calcWeightEverythingElse(P_new, S_new, wingWeight, A, C, P_rshp, k_pg, W_TO):
    '''
    DESCRIPTION:
    Simple calculator for the "everything else" weight in de Vries sources, assuming twin engine aircraft
    
    INPUT:
        P
        S
        wingWeight
        A
        C
        P_rshp
        k_pg
        W_TO

    OUTPUT:
        W_e

    '''
    #Calculate wing weight
    #W_wing = wingWeight*S_new   #lbs

    #Calculate wing everything else prime
    W_e = A*W_TO**C*W_TO       #lbs

    #Ugliness ensues

    #3.1.3 Calculate Powergroup weight (powerplant)

    W_eng = P_rshp**(0.9306)*10**(-0.1205)  #Weight of engine, lbs
    W_pg = k_pg*(W_eng + 0.24*P_new)            #Weight of power group -> assume powerplant weight

    W_ep_plus_W_wing =  W_pg*2 - W_e 

    W_e_new =  W_ep_plus_W_wing + W_pt

    #Calculate weight everything else
    #W_e = W_pt + W_pg*2 + W_wing + W_ep   #lbs

    return W_e_new

#======================================================================
#Create functions for loop (Step 3.2)
def calcWeightEnergy(Phi, g, eta_1, eta_2, eta_3, e_bat, e_fuel, L_D, W_e, W_payload):
    '''
    DESCRIPTION:
    Simple calculator for finding fuel and battery weights of a hybrid electric aircrft
    
    INPUT:
        Phi     -   Hybridization Ratio
        g       -   gravity, ft/s^2
        eta_1   -   overall fuel powertrain efficiency
        eta_2   -   overall battery powertrain efficiency
        eta_3   -   overall fuel powertrain efficiency
        eta_bat -   battery specific energy
        eta_fuel-   fuel specific energy
        L_D     -   Lift-to-Drag ratio

    OUTPUT:
        W_bat   -   
        W_fuel  -

    '''



    A = eta_3*(e_fuel*550/g)*L_D*(eta_1 + eta_2*(Phi/(1 - Phi)))
    B = W_e + W_payload
    Cp1 = (g/(e_bat*550))*(Phi + e_bat/e_fuel*(1 - Phi))
    Cp2 = g/(e_bat*550)

    E_o_tot = (B-B*np.exp(R/A))/(Cp2*np.exp(R/A) - Cp1)

    #Solve for battery/fuel weight
    W_bat = (g/(e_bat*0.8*550))*(Phi*E_o_tot)
    W_fuel = (g/(e_fuel*550))*((1 - Phi)*E_o_tot)      
    
    return W_bat, W_fuel

#======================================================================
#Power, Wing Area, Weight approximation for Parallel Hybrid:
tol = 0.1
res = 1
p = 0
while res > tol:
    p = p+1
    P_new, S_new = calcWPS(W_TO, W_P, W_S)      

    W_e = calcWeightEverythingElse(P_new, S_new, wingWeight, A, C, P_rshp, k_pg, W_TO)

    W_bat, W_fuel = calcWeightEnergy(Phi, g, eta_1p, eta_2p, eta_3p, e_bat, e_fuel, L_D, W_e, W_payload)

    W_TO_new = W_crew + W_payload + W_fuel + W_bat + W_e
    
    res = np.abs(W_TO_new - W_TO)
    print(P_new, S_new)
    print(W_TO)

    print(p)
    W_TO = W_TO_new
P_new, S_new = calcWPS(W_TO, W_P, W_S) 

print('=============================')
print('Summary of Results')
print('=============================')
print('Take off weight (lbs): %0.3f' % W_TO)
print('Take off power (lbf-ft): %0.3f' % P_new)
print('New Wing Area (sqft): %0.3f' % S_new)
