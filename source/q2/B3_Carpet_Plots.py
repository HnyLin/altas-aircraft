import numpy as np
import plotly.graph_objects as go

def Get_AR_CD0_Carpet(V_cruise, AR_vals, C_D0_vals, e):
    '''
    AR_vals & C_D0_vals must be same size
    '''
    rho_interp = [0.0765, 0.0565, 0.0408, 0.0287, 0.0189]
    h_interp = [0, 10000, 20000, 30000, 40000]
    rho = np.interp(h_cruise, h_interp, rho_interp)

    V = V_cruise * 1.688
    LD_full = []
    WS_full = []

    for AR in AR_vals:
        LD_vals = []
        WS_vals = []
        
        for C_D0 in C_D0_vals:
            #Carpet Plot AR vs CD_0 Trade Study
            LD_max = 0.5 * np.sqrt(np.pi * AR * e / C_D0 )

            C_L = np.sqrt(C_D0 * np.pi * AR * e)
            WS = 0.5 * rho * V ** 2 * C_L / 32.174

            LD_vals = np.append(LD_vals, LD_max)
            WS_vals = np.append(WS_vals, WS)

        LD_full = np.append(LD_full, LD_vals)  
        WS_full = np.append(WS_full, WS_vals)

    AR_vals = np.repeat(AR_vals, n)
    C_D0_vals = np.tile(C_D0_vals, n)


    fig = go.Figure(go.Carpet(
        a = AR_vals,
        b = C_D0_vals,
        y = LD_full,
        x = WS_full,
        aaxis=dict(
            tickprefix='AR = ',
            smoothing=0.2,
        ),
        baxis=dict(
            tickprefix='CD_0 = ',
            smoothing=0.4,
        )
    ))

    fig.update_layout(
        xaxis=dict(
            tickprefix = 'W/S =',
            showgrid=True,
            showticklabels=True
        ),
        yaxis=dict(
        tickprefix = 'L/D = ',
            showgrid=True,
            showticklabels=True,
        )
    )
    
    fig.show()
    return

#Assuming Cruise Conditions
h_cruise = 28000
e = 0.8
V_cruise = 350
n = 7
AR_vals = np.linspace(5, 20, n)
C_D0_vals = np.linspace(0.005, 0.005*n, n)
Get_AR_CD0_Carpet(V_cruise, AR_vals, C_D0_vals, e)