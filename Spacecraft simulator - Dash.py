#This is a simple spacecraft simulator solution developed for a homework task
#Credit goes to Prof. Dr. B. Dachwald for the main calculations and the task

import numpy as np
import scipy.integrate

from dash import Dash, html, dcc, callback, Output, Input
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

simlen = 2000  #Simulation lenght in h

#Basic values taken from BepiColombo

𝛼 = 90   #control variable in °
#FT = 10  #Thrust in N
#c = 2 #exhaust velocity in km/s

rE = 6378   #Earth radius in km
a0 = 300    #Initial orbit altitude
r0 = a0+rE   #Initial orbit radius in km
v0 = 0   #Initial orbital velocity in km/s (set zero to use Eccentricity)
E = 0.2    #Eccentricity
mD = 2700   #Drymass in kg
mF0 = 1400  #Initial Fuelmass
m0 = mD + mF0
µ = 3.986e5 #Gtravitational parameter in km3/s2

if v0 == 0:
    v0 = np.sqrt(µ/r0) #For a circular orbit
    ai = r0/(1-E) #Semi-major axis
    v0 = (2*µ/r0 - µ/ai)**0.5

def Event_FuelEmpty(t, y, FT, µ, 𝛼, c): #Check if Fuel is empty, used in solve_IVP as  breakpoint when = 0
    r, φ, pr, ρφ, m = y
    return m - mD
Event_FuelEmpty.terminal = True

def motionEQ(t, y, FT, µ, 𝛼, c): #y = [r, φ, pr, ρφ, m]
    r, φ, pr, ρφ, m = y
    r_der = pr
    φ_der = ρφ
    pr_der = r * (ρφ**2) - µ/(r**2) + (FT/m)*np.cos(𝛼)
    ρφ_der =  (1/r) * ((FT/m)*np.sin(𝛼) - 2 * pr * ρφ)
    m_der = - FT/c

    return [r_der, φ_der, pr_der, ρφ_der, m_der]

def updateIVP(FT, c, simlen):
    #Calculates the Flightpath. Terminates at 0 fuel or after simlen
    sol = scipy.integrate.solve_ivp(motionEQ, [0, simlen], [r0, 0, 0, v0/r0, m0], args=(FT, µ, 𝛼, c), t_eval= np.linspace(0, simlen, simlen), events=Event_FuelEmpty)

    xSol = sol.y[0] * np.cos(sol.y[1])
    ySol = sol.y[0] * np.sin(sol.y[1])

    FuelDepleted = 0

    P_FinalOrbit = calcOrbitalElements(sol.y[0][-1], sol.y[2][-1], sol.y[3][-1])[1]

    #Calculate the final orbit
    if sol.t_events[0].size > 0:
        FuelDepleted = sol.t_events[0][0]
        #solBal = scipy.integrate.solve_ivp(motionEQ, [int(sol.t[-1]), simlen], [sol.y[0][-1], sol.y[1][-1], sol.y[2][-1], sol.y[3][-1], sol.y[4][-1]], args=(0, µ, 𝛼, c), t_eval= np.linspace(int(sol.t[-1]), simlen, (simlen-int(sol.t[-1]))))
        solBal = scipy.integrate.solve_ivp(motionEQ, [0, int(P_FinalOrbit)], [sol.y[0][-1], sol.y[1][-1], sol.y[2][-1], sol.y[3][-1], sol.y[4][-1]], args=(0, µ, 𝛼, c), t_eval= np.linspace(0, int(P_FinalOrbit), int(P_FinalOrbit)))
        xSol = np.append(xSol, solBal.y[0] * np.cos(solBal.y[1]))
        ySol = np.append(ySol, solBal.y[0] * np.sin(solBal.y[1]))


    return sol, xSol, ySol, FuelDepleted

def calcOrbitalElements(r, pr, pφ):
    A_FinalOrbit = 1/((2/r) - (pr**2 + (r * pφ)**2)/µ)
    P_FinalOrbit = 2 * np.pi * np.sqrt((A_FinalOrbit**3)/µ)

    return A_FinalOrbit, P_FinalOrbit

#Beginning of App code
app = Dash()

app.layout = [
    html.H1(children='Spacecraft simulator', style={'textAlign':'center'}),
    html.Br(),
    html.Label('Thrust[N]'),
    dcc.Slider(
            min=0,
            max=4,
            marks={i: str(i) for i in range(1, 9)},
            value=0.29,
            id='Thrust-setting',
            updatemode='drag',),
    
    html.Label('Exhaust Velocity[km/s]'),
    dcc.Slider(
            min=0.1,
            max=10,
            value=4.3,
            id='Exhaust-setting',
            updatemode='drag',),

    dcc.Graph(id='graph-content')
]

@callback(
    Output('graph-content', 'figure'),
    Input('Thrust-setting', 'value'),
    Input('Exhaust-setting', 'value'),
)

def update_graph(FT,c):

    maxSimlen = 3600 * 24 *14 #Maximum 14 Tage

    sol, xSol, ySol, FuelDepletedTime = updateIVP(float(FT),float(c),maxSimlen)

    fig = go.Figure()
    fig.add_trace(go.Line(x=xSol[:int(FuelDepletedTime)], y=ySol[:int(FuelDepletedTime)], name='Spacecraft Trajectory'))

    if FuelDepletedTime  != 0:
        #Marker point at which fuel runs out
        fig.add_trace(go.Scatter(
            x=[xSol[int(FuelDepletedTime)-1]],  # X-coordinate of Fuel cutoff point
            y=[ySol[int(FuelDepletedTime)-1]], 
            mode='markers',
            marker=dict(size=12, color='red', symbol='star'),
            name='No Fuel'
            ))
        
        #Final orbit
        radius =  np.sqrt(xSol[int(FuelDepletedTime):]**2 + ySol[int(FuelDepletedTime):]**2)

        fig.add_trace(go.Line(
            x=xSol[int(FuelDepletedTime):], 
            y=ySol[int(FuelDepletedTime):], 
            mode='lines', 
            marker=dict(color='red'),
            name = "Ending Orbit",
            hovertemplate=(
                "<b>X:</b> %{x}<br>"
                "<b>Y:</b> %{y:.2f}<br>"
                "<b>Radius:</b> %{customdata:.2f}<extra></extra>"),
            customdata=radius
            ))
        
        #Anotate with final orbit infomation
        fig.add_annotation(
            x=xSol[int(FuelDepletedTime)-1],  
            y=ySol[int(FuelDepletedTime)-1],  
            text=f"Orbital Period: {calcOrbitalElements(sol.y[0][int(FuelDepletedTime)], sol.y[2][int(FuelDepletedTime)], sol.y[3][int(FuelDepletedTime)])[0]/(3600):.3f}h", 
            showarrow=True,  
            arrowhead=2,  
            ax=50,  # X-offset for the arrow
            ay=-50,  # Y-offset for the arrow
            font=dict(color="black", size=12)
        )
    
    fig.add_shape(type="circle",
        xref="x", yref="y",
        x0=-rE, y0=-rE, x1=rE, y1=rE,
        fillcolor="PaleTurquoise",
        line_color="LightSeaGreen",
    )

    fig.update_layout(
        xaxis=dict(scaleanchor="y"),  # Link the scale of x-axis to y-axis
        yaxis=dict(),
        height=1000,
    )
    
    return fig

if __name__ == '__main__':
    app.run(debug=True)