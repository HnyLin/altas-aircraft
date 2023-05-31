import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import time
import pandas as pd
import plotly.graph_objects as go
from labellines import labelLines
import math
from math import log10, floor
from sympy import symbols, Eq, solve

def get_WP(MTOW, Wing_area, AR, h_cruise):

    #Variables for W_P Calculation
    lbf_lbm = 32.17
    hp = 550                #550 ft*lbf / s

    precision = 300
    WS = np.linspace(1, 300, precision)


    nu_p = 0.8                                  #From Report
    #AR = 10.06                                  #From Design
    w_to = MTOW                                 #lbs from weight estimate
    w_L = MTOW                                  #lbs from weight estimate
    wing_area = Wing_area

    BFL = 4500                                  #From RFP

    #Inputs (Takeoff)
    rho_SL = 0.0025332 * lbf_lbm
    rho = 20.48 * 10 ** -4 * lbf_lbm            #Air density at 5000'
    C_lmax_TO = 1.7                             #Roskam Vol 1 Table 3.1
    BPR = 50                                    #Verify BPR for Turbofan
    k_2 = 0.75 * (5 + BPR) / (4 + BPR)
    k_1 = 0.0447
    C_d0_TO = 0.045                             #From Drag Polars
    mu_g = 0.025                                #From Roskam

    #Inputs (Landing)
    S_FL = 4500
    S_L = 0.6 * S_FL
    S_a = 1000                                  #Metabook eq 4.46 (ft)
    C_lmax_L = 3.3                              #Roskam Vol 1 Table 3.1

    #Inputs(Cruise)
    rho_interp = [0.0765, 0.0565, 0.0408, 0.0287, 0.0189]
    h_interp = [0, 10000, 20000, 30000, 40000]
    rho_cruise = np.interp(h_cruise, h_interp, rho_interp)
    V_cruise = 350 * 1.688                      #Declared in the beginning for W/P
    C_D0_cruise = 0.025                         #Drag Polar Estimate
    e_cruise = 0.75                             #Drag Polar Estimate

    #Inputs (Ceiling)
    G_ceiling = 0.001
    e = 0.8
    rho_ceiling = 8.91 * 10 ** -4   * lbf_lbm        #air density at 30,000   

    #FAR 25.121 Climb Inputs
    #Takeoff Climb
    #Transition Segment Climb
    #Second Segment Climb
    #Enroute Climb
    #Balked Landing Climb
    C_Lmax = 1.7
    G_climbvals = [0.012, 0, 0.024, 0.012, 0.032, 0.021]
    k_svals = [1.2, 1.15, 1.2, 1.25, 1.3, 1.5]
    C_D0vals = [0.045 + 0.025, 0.045 + 0.025, 0.045, 0.025, 0.125 + 0.1, (0.125 + 0.1) / 2]
    C_Lmax_vals = [C_lmax_TO, C_lmax_TO, C_lmax_TO, C_Lmax, C_lmax_L, 0.85 * C_lmax_L]
    labels = ["Takeoff Climb", "Transition Segment Climb", "Second Segment Climb", 
            "Enroute Climb", "Balked Landing Climb AEO", "Balked Landing Climb OEI"]
    e_vals = [0.75, 0.75, 0.75, 0.8, 0.7, 0.7]
    w_vals = [0.95 * w_to,  0.94 * w_to, 0.93 * w_to, 0.92 * w_to, 1 * w_to, 1 * w_to]            #Made up Weights
    N_correction = [1, 1, 1, 1, 1, 2]


    def tw_wp( WS, nu_p, tw, C_L, rho ):
        lbf_lbm = 32.17
        V_cruise = np.sqrt ( 2 / ( rho * C_L ) * WS * lbf_lbm)
        wp = nu_p / V_cruise * (1 / tw)
        return wp

    def get_Restraint(precision, lbf_lbm, hp, WS, nu_p, AR, w_to, w_L, BFL, rho_SL, rho, 
                    C_lmax_TO, k_1, k_2, C_d0_TO, mu_g, S_FL, S_L, S_a, C_lmax_L, 
                    rho_cruise, V_cruise, C_D0_cruise, e_cruise, G_ceiling, e, 
                    rho_ceiling, C_Lmax, G_climbvals, k_svals, C_D0vals, 
                    C_Lmax_vals, labels, e_vals, w_vals, N_correction, wing_area):
        
        ##plt.figure(figsize = (16,10) )
        #plt.xlabel("W/S (lbf / ft^2)")
        #plt.ylabel("W/P (lbf / hp)")

        #Takeoff
        #Sea Level 30ft Obsticle 
        TOP = BFL / 37.5                            #Check Units
        tw_TO = WS / ( (rho_SL / rho_SL) * C_lmax_TO * TOP)
        
        wp_TOSL30 = hp * tw_wp( WS, nu_p, tw = tw_TO, C_L = C_lmax_TO, rho = rho_SL)
        
        #plt.plot(WS, wp_TOSL30, label = 'Takeoff Field Length SL 30ft Obsticle')
        
        #Sea Level 50ft Obsticle(SL+ISA + 18 deg day)
        
        tw_TO = 1 / ( k_2 / ( ( k_1 * WS / ( BFL * rho_SL) + 0.72 * C_d0_TO ) / C_lmax_TO + mu_g ))
        
        wp_TOSL50 = hp * tw_wp( WS, nu_p, tw = tw_TO, C_L = C_lmax_TO, rho = rho_SL)
        #plt.plot(WS, wp_TOSL50, label = 'Takeoff Field Length SL+ISA 50ft Obsticle')
        
        #5000 + Sea Level 50ft Obsticle(5000+ISA + 18 deg day)
        tw_TO = 1 / ( k_2 / ( ( k_1 * WS / ( BFL * rho) + 0.72 * C_d0_TO ) / C_lmax_TO + mu_g ))
        
        wp_TO5K50 = hp * tw_wp( WS, nu_p, tw = tw_TO, C_L = C_lmax_TO, rho = rho)
        #plt.plot(WS, wp_TO5K50, label = 'Takeoff Field Length 5000+ISA 50ft Obsticle')
        
        #Landing
        #Sea Level 30ft Obsticle
        ws_L = ( ( rho_SL / rho_SL) * C_lmax_L ) / 80 * (S_L - S_a)
        ws_L = ws_L / (w_L / w_to)
        #plt.plot(ws_L * np.ones(precision), np.linspace(0, precision, precision), label = "Landing Field Length")
        
        #Cruise
        q = rho_cruise * V_cruise ** 2 / 2
        
        tw_cruise = q * C_D0_cruise / WS + WS / ( q * np.pi * AR * e_cruise)
        
        wp_cruise = nu_p / V_cruise * 1 / tw_cruise * lbf_lbm * hp
        
        #plt.plot( WS, wp_cruise, label = "Cruise" )
        
        #Stall
        V_takeoff = np.sqrt( 2 * w_to / (C_lmax_TO * wing_area * rho_SL) )
        V_stall = 1/ 1.1 * V_takeoff
        ws_stall = 1 / 2 * rho * V_stall ** 2 * C_Lmax
        
        #plt.plot(ws_stall * np.ones(precision), np.linspace(0, precision, precision), label = "Stall Requirement")
        
        #FAR 25.121 Climb Requirements
        temp_vals = []
        for i in range(len(G_climbvals)):
            
            tw_C_uncorrected = k_svals[i] ** 2 / C_Lmax_vals[i] * C_D0vals[i] + C_Lmax_vals[i] / ( k_svals[i] ** 2 * np.pi * e_vals[i] * AR ) + G_climbvals[i] 
            
            tw_C = ( 1 / 0.8 ) * ( 1 / 0.94 ) * N_correction[i] * w_vals[i] / w_to * tw_C_uncorrected
                
            wp_C = hp * tw_wp( WS, nu_p, tw = tw_C, C_L = C_Lmax_vals[i], rho = rho_SL )
            
            temp_vals = np.append(temp_vals, wp_C)
            
            #plt.plot(WS, wp_C, label = labels[i])
        
        wp_takeoffclimb = temp_vals[0 : precision]
        wp_transegmentclimb = temp_vals[precision : 2 * precision]
        wp_secsegmentclimb = temp_vals[2 * precision : 3 * precision]
        wp_enrouteclimb = temp_vals[3 * precision : 4 * precision]
        wp_balkedAEO = temp_vals[4 * precision : 5 * precision]
        wp_balkedOEI = temp_vals[5 * precision : ]

        #Ceiling
        
        tw_celing = (1 / (rho_ceiling / rho_SL) ** 0.6) * (2 * np.sqrt(C_D0_cruise / (np.pi * AR * e)) + G_ceiling)
        
        wp_celing = hp * tw_wp( WS, nu_p, tw = tw_celing, C_L = C_Lmax, rho = rho_ceiling )
        
        #plt.plot(WS, wp_celing, label = "Ceiling")
        
        #plt.title("W/P - W/S")
        
        ##plt.legend(loc = 'best')
        #plt.ylim(0, 40)
        #plt.xlim(0, 100)
        
        #labelLines(plt.gca().get_lines(), zorder=2.5)
        
        #loc = np.where( (wp_balkedOEI - wp_cruise) < 0 )
        #loc = max(loc)
        
        
        #plt.fill_between(WS, wp_balkedOEI, where = ( WS < ws_L ) & (wp_balkedOEI < wp_cruise), color = "green", alpha = 0.2)
        #plt.fill_between(WS, wp_cruise, where = WS < WS[min(loc) + 1], color = "green", alpha = 0.2)
        
        #plt.show()

        #Serial Efficency Values
        nu_GT = 0.3
        nu_GB = 0.96                     #None
        nu_EM1 = 0.96
        nu_PM = 0.99
        nu_EM2 = 0.96
        nu_P2 = 0.85

        sshaft_PR = 1                #Guess
        ssupplied_PR = 0.25          #Guess

        #Parallel Efficency Values
        nu_GT = 0.3
        nu_GB = 0.96
        nu_P1 = 0.9
        nu_EM1 = 0.96
        nu_PM = 0.99     

        pshaft_PR = 0            #Guess
        psupplied_PR = 0.25         #Guess


        labels = ["Takeoff Sea Level 30ft", "Takeoff Sea Level 50ft", 
                "Takeoff 5000ft 50ft", "Cruise", "Takeoff Climb", 
                "Transition Segment Climb", "Second Segment Climb", "Enroute Climb", 
                "Balked Landing Climb AEO", "Balked Landing Climb OEI", "Ceiling"]

        A_serial = np.array([[-nu_GT, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                            [0, -nu_GB, 1, 1, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, -nu_P1, 0, 0, 0, 0, 1, 0],
                            [0, 0, -nu_EM1, 0, 1, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, -nu_PM, -nu_PM, 1, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, -nu_EM2, 1, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, -nu_P2, 0, 1],
                            [ssupplied_PR, 0, 0, 0, 0, (ssupplied_PR - 1), 0, 0, 0, 0],
                            [0, 0, 0, sshaft_PR, 0, 0, 0, (sshaft_PR - 1), 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0, 1, 1]
                            ])

        A_parallel = np.array([[-nu_GT, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                            [0, -nu_GB, 1, 1, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, -nu_P1, 0, 0, 0, 0, 1, 0],
                            [0, 0, -nu_EM1, 0, 1, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, -nu_PM, -nu_PM, 1, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, -nu_EM2, 1, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, -nu_P2, 0, 1],
                            [psupplied_PR, 0, 0, 0, 0, (psupplied_PR - 1), 0, 0, 0, 0],
                            [0, 0, 0, pshaft_PR, 0, 0, 0, (pshaft_PR - 1), 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0, 1, 1]
                            ])

        pw = [wp_TOSL30, wp_TOSL50, wp_TO5K50, wp_cruise, wp_takeoffclimb, wp_transegmentclimb, \
            wp_secsegmentclimb, wp_enrouteclimb, wp_balkedAEO, wp_balkedOEI, \
                wp_celing]
            
        pw = 1 / np.array(pw)

        #Serial Calculation
        x_matrix_serial = []
        for i in range(len(labels)):
            b = np.zeros((10,1))
            x_vals = np.ones([10, precision])
            for j in range(precision):
                b[-1] = pw[i, j]
                x = np.linalg.solve(A_serial, b)
                x_vals[:, j] = x.T
            
            x_matrix_serial.append(x_vals)

        #Parallel Calculation
        x_matrix_parallel = []
        for i in range(len(labels)):
            b = np.zeros((10,1))
            x_vals = np.ones([10, precision])
            for j in range(precision):
                b[-1] = pw[i, j]
                x = np.linalg.solve(A_parallel, b)
                x_vals[:, j] = x.T
            
            x_matrix_parallel.append(x_vals)
            
        #Extracting Data for Pgt and Pem1
        gtloc = 1
        em1loc = 4
        batloc = 5

        P_gt_serial = np.ones([len(labels), precision])
        P_em1_serial = np.ones([len(labels), precision])
        P_bat_serial = np.ones([len(labels), precision])

        P_gt_parallel = np.ones([len(labels), precision])
        P_em1_parallel = np.ones([len(labels), precision])
        P_bat_parallel = np.ones([len(labels), precision])

        for i in range(len(labels)):
            x_temp = x_matrix_serial[i]
            P_gt_serial[i, :] = x_temp[gtloc, :]
            P_em1_serial[i, :] = x_temp[em1loc, :]
            P_bat_serial[i, :] = x_temp[batloc, :]
            
            x_temp = x_matrix_parallel[i]
            P_gt_parallel[i, :] = x_temp[gtloc, :]
            P_em1_parallel[i, :] = x_temp[em1loc, :]
            P_bat_parallel[i, :] = x_temp[batloc, :]

        '''
        #Plotting Battery Serial
        plt.figure(figsize = (12, 12))
        plt.xlabel("W/S (lbf / ft^2)")
        plt.ylabel("W/P (lbf / hp)")
        plt.title("P_bat Serial Configuration")
        plt.ylim(0, 40)
        plt.xlim(0,100)

        for i in range(len(labels)):
            plt.plot(WS, 1 / P_bat_serial[i, :], label = labels[i])

        plt.plot(ws_L * np.ones(precision), np.linspace(0, precision, precision), label = "Landing Field Length")

        labelLines(plt.gca().get_lines(), zorder=2.5)

        loc = np.where( (P_bat_serial[3,:] - P_bat_serial[9,:]) < 0 )

        plt.fill_between(WS, 1 / P_bat_serial[3, :], where = WS < min(WS[loc]), color = 'green', alpha = 0.2 )

        plt.fill_between(WS, 1 / P_bat_serial[9,:], where = (WS > min(WS[loc]) - 2 ) & (WS < ws_L), color = 'green', alpha = 0.2  )

        plt.show()

        #Plotting Battery Parallel
        plt.figure(figsize = (12, 12))
        plt.xlabel("W/S (lbf / ft^2)")
        plt.ylabel("W/P (lbf / hp)")
        plt.title("P_bat Parallel Configuration")
        plt.ylim(0, 40)
        plt.xlim(0,100)

        for i in range(len(labels)):
            plt.plot(WS, 1 / P_bat_parallel[i, :], label = labels[i])

        plt.plot(ws_L * np.ones(precision), np.linspace(0, precision, precision), label = "Landing Field Length")

        labelLines(plt.gca().get_lines(), zorder=2.5)

        loc = np.where( (P_bat_parallel[3,:] - P_bat_parallel[9,:]) < 0 )

        plt.fill_between(WS, 1 / P_bat_parallel[3, :], where = WS < min(WS[loc]), color = 'green', alpha = 0.2 )

        plt.fill_between(WS, 1 / P_bat_parallel[9,:], where = (WS > min(WS[loc]) - 2 ) & (WS < ws_L), color = 'green', alpha = 0.2  )

        plt.show()
            
            
        #Plotting Pgt Serial
        plt.figure(figsize = (12, 12))
        plt.xlabel("W/S (lbf / ft^2)")
        plt.ylabel("W/P (lbf / hp)")
        plt.title("P_gt Serial Configuration")
        plt.ylim(0, 40)
        plt.xlim(0,100)

        for i in range(len(labels)):
            plt.plot(WS, 1 / P_gt_serial[i, :], label = labels[i])
            
        plt.plot(ws_L * np.ones(precision), np.linspace(0, precision, precision), label = "Landing Field Length")
            
        labelLines(plt.gca().get_lines(), zorder=2.5)

        loc = np.where( (P_gt_serial[3,:] - P_gt_serial[9,:]) < 0 )

        plt.fill_between(WS, 1 / P_gt_serial[3, :], where = WS < min(WS[loc]), color = 'green', alpha = 0.2 )

        plt.fill_between(WS, 1 / P_gt_serial[9,:], where = (WS > min(WS[loc]) - 2 ) & (WS < ws_L), color = 'green', alpha = 0.2  )

        plt.show()

        #Plotting Pgt Parallel
        plt.figure(figsize = (12, 12))
        plt.xlabel("W/S (lbf / ft^2)")
        plt.ylabel("W/P (lbf / hp)")
        plt.title("P_gt Parallel Configuration")
        plt.ylim(0, 40)
        plt.xlim(0,100)

        for i in range(len(labels)):
            plt.plot(WS, 1 / P_gt_parallel[i, :], label = labels[i])
            
        plt.plot(ws_L * np.ones(precision), np.linspace(0, precision, precision), label = "Landing Field Length")
            
        labelLines(plt.gca().get_lines(), zorder=2.5)

        loc = np.where( (P_gt_parallel[3,:] - P_gt_parallel[9,:]) < 0 )

        plt.fill_between(WS, 1 / P_gt_parallel[3, :], where = WS < min(WS[loc]), color = 'green', alpha = 0.2 )

        plt.fill_between(WS, 1 / P_gt_parallel[9,:], where = (WS > min(WS[loc]) - 2 ) & (WS < ws_L), color = 'green', alpha = 0.2  )


        plt.show()

        #Plotting Pem1 Serial
        plt.figure(figsize = (12, 12))
        plt.xlabel("W/S (lbf / ft^2)")
        plt.ylabel("W/P (lbf / hp)")
        plt.title("P_em1 Serial Configuration")
        plt.ylim(0, 40)
        plt.xlim(0,100)

        for i in range(len(labels)):
            plt.plot(WS, abs(1 / P_em1_serial[i, :]), label = labels[i])
            
        plt.plot(ws_L * np.ones(precision), np.linspace(0, precision, precision), label = "Landing Field Length")

        labelLines(plt.gca().get_lines(), zorder=2.5)

        loc = np.where( (P_em1_serial[3,:] - P_em1_serial[9,:]) < 0 )

        plt.fill_between(WS, 1 / P_em1_serial[3, :], where = WS < min(WS[loc]), color = 'green', alpha = 0.2 )

        plt.fill_between(WS, 1 / P_em1_serial[9,:], where = (WS > min(WS[loc]) - 2 ) & (WS < ws_L), color = 'green', alpha = 0.2  )

        plt.show()

        #Plotting Pem1 Parallel
        plt.figure(figsize = (12, 12))
        plt.xlabel("W/S (lbf / ft^2)")
        plt.ylabel("W/P (lbf / hp)")
        plt.title("P_em1 Parallel Configuration")
        plt.ylim(0, 40)
        plt.xlim(0,100)


        for i in range(len(labels)):
            plt.plot(WS, abs(1 / P_em1_parallel[i, :]), label = labels[i])
            
        plt.plot(ws_L * np.ones(precision), np.linspace(0, precision, precision), label = "Landing Field Length")

        labelLines(plt.gca().get_lines(), zorder=2.5)

        loc = np.where( (abs(P_em1_parallel[3,:]) - abs(P_em1_parallel[9,:])) < 0 )

        plt.fill_between(WS, abs(1 / P_em1_parallel[3, :]), where = WS < min(WS[loc]), color = 'green', alpha = 0.2 )

        plt.fill_between(WS, abs(1 / P_em1_parallel[9,:]), where = (WS > min(WS[loc]) - 2 ) & (WS < ws_L), color = 'green', alpha = 0.2  )

        plt.show()

        '''

        return P_gt_parallel, P_em1_parallel, P_bat_parallel, ws_L
        
    P_gt_parallel, P_em1_parallel, P_bat_parallel, ws_L = get_Restraint(precision, lbf_lbm, hp, WS, nu_p, AR, w_to, w_L, BFL, rho_SL, rho, 
                    C_lmax_TO, k_1, k_2, C_d0_TO, mu_g, S_FL, S_L, S_a, C_lmax_L, 
                    rho_cruise, V_cruise, C_D0_cruise, e_cruise, G_ceiling, e, 
                    rho_ceiling, C_Lmax, G_climbvals, k_svals, C_D0vals, 
                    C_Lmax_vals, labels, e_vals, w_vals, N_correction, wing_area)
    
    OEI_loc = 9

    #Parallel Configuation Chosen
    P_gt_OEI = 1 / P_gt_parallel[ OEI_loc, : ]
    P_em1_OEI = 1 / P_em1_parallel[ OEI_loc, : ]
    P_bat_OEI = 1 / P_bat_parallel[ OEI_loc, : ]
    
    #Interpolate to Determine WP at max WS
    WP_gt = np.interp(ws_L, WS, P_gt_OEI)
    WP_em1 = np.interp(ws_L, WS, P_em1_OEI)
    WP_bat = np.interp(ws_L, WS, P_bat_OEI)

    W_P = max(WP_gt, WP_em1, WP_bat)

    return W_P

MTOW = 82561
Wing_area = 800
h_cruise = 28000
AR = 13
W_P  = get_WP(MTOW, Wing_area, AR, h_cruise)
