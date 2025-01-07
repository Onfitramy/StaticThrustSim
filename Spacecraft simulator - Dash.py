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

𝛼 = np.deg2rad(90)   #control variable in °
#FT = 0.29  #Thrust in N
#c = 2 #exhaust velocity in km/s

rE = 6378   #Earth radius in km
a0 = 300    #Initial orbit altitude
r0 = a0+rE   #Initial orbit radius in km
v0 = 0   #Initial orbital velocity in km/s (set zero to use Eccentricity)
E = 0    #Eccentricity
mD = 2700   #Drymass in kg
mF0 = 1400  #Initial Fuelmass
µ = 3.986e5 #Gtravitational parameter in km3/s2

if v0 == 0:
    v0 = np.sqrt(µ/r0) #For a circular orbit
    ai = r0/(1-E) #Semi-major axis
    v0 = (2*µ/r0 - µ/ai)**0.5

def Event_FuelEmpty(t, y, FT, µ, 𝛼, c, mD): #Check if Fuel is empty, used in solve_IVP as  breakpoint when = 0
    r, φ, pr, ρφ, m = y
    return m - mD
Event_FuelEmpty.terminal = True

def motionEQ(t, y, FT, µ, 𝛼, c, mD): #y = [r, φ, pr, ρφ, m]
    r, φ, pr, ρφ, m = y
    #𝛼 = np.arctan((r*ρφ)/pr)
    r_der = pr
    φ_der = ρφ
    pr_der = r * (ρφ**2) - µ/(r**2) + (FT/m)*np.cos(𝛼)
    ρφ_der =  (1/r) * ((FT/m)*np.sin(𝛼) - 2 * pr * ρφ)
    m_der = - FT/c

    return [r_der, φ_der, pr_der, ρφ_der, m_der]

def updateIVP(FT, c, simlen, DryMass, FuelMass):

    m0 = DryMass + FuelMass

    #Calculates the Flightpath. Terminates at 0 fuel or after simlen
    sol = scipy.integrate.solve_ivp(motionEQ, [0, simlen], [r0, 0, 0, v0/r0, m0], args=(FT, µ, 𝛼, c, DryMass), t_eval= np.linspace(0, simlen, simlen), events=Event_FuelEmpty)

    #Convert to 
    xSol = sol.y[0] * np.cos(sol.y[1])
    ySol = sol.y[0] * np.sin(sol.y[1])

    FuelDepleted = 0

    P_FinalOrbit = calcOrbitalElements(sol.y[0][-1], sol.y[2][-1], sol.y[3][-1])[1]

    #Calculate the final orbit
    if sol.t_events[0].size > 0:
        FuelDepleted = sol.t_events[0][0]
        solBal = scipy.integrate.solve_ivp(motionEQ, [0, int(P_FinalOrbit)], [sol.y[0][-1], sol.y[1][-1], sol.y[2][-1], sol.y[3][-1], sol.y[4][-1]], args=(0, µ, 𝛼, c, DryMass), t_eval= np.linspace(0, int(P_FinalOrbit), int(P_FinalOrbit)))
        xSol = np.append(xSol, solBal.y[0] * np.cos(solBal.y[1]))
        ySol = np.append(ySol, solBal.y[0] * np.sin(solBal.y[1]))
    else:
        FuelDepleted = sol.t[-1]


    return sol, xSol, ySol, FuelDepleted

def calcOrbitalElements(r, pr, pφ):
    A_FinalOrbit = 1/((2/r) - (pr**2 + (r * pφ)**2)/µ)
    if A_FinalOrbit > 0:                                                #Orbit is not parabolic or Hyperbolic
        P_FinalOrbit = 2 * np.pi * np.sqrt((A_FinalOrbit**3)/µ)
    else:
        P_FinalOrbit =  3600*24*2   #Last orbit simulation goes for 2 days

    return A_FinalOrbit, P_FinalOrbit

#Beginning of App code
app = Dash()

app.layout = [
    html.H1(children='Spacecraft simulator', style={'textAlign':'center'}),
    html.Br(),
    html.Label('Thrust[N]'),
    dcc.Slider(
            min=0,
            max=0.1,
            marks={i: str(i) for i in range(1, 9)},
            value=0.01,             #!!Startwerte!!
            id='Thrust-setting',
            updatemode='drag',),
    
    html.Label('Exhaust Velocity[km/s]'),
    dcc.Slider(
            min=1,
            max=100,
            value=42.168,           #!!Startwerte!!
            id='Exhaust-setting',
            updatemode='drag',),
    
    html.Label('Dry Mass[kg]'),
    dcc.Slider(
            min=100,
            max=10000,
            value=2700,             #!!Startwerte!!
            id='DryMass-setting',
            updatemode='drag',),
    
    html.Label('Fuel Mass[kg]'),
    dcc.Slider(
            min=100,
            max=10000,
            value=700,              #!!Startwerte!!
            id='FuelMass-setting',
            updatemode='drag',),

    html.Div(id='dV_readout', style={'whiteSpace': 'pre-line'}),

    dcc.Checklist(
    ['Thrust in Velocity Vector'],
    ['Thrust in Velocity Vector'],
    inline=True,
    id="Checklist"
    ),

    dcc.Graph(id='graph-content') 
]

@callback(
    Output('graph-content', 'figure'),
    Output('dV_readout', 'children'),
    Input('Thrust-setting', 'value'),
    Input('Exhaust-setting', 'value'),
    Input('DryMass-setting', 'value'),
    Input('FuelMass-setting', 'value'),
    Input('Checklist', 'value'),
)

def update_graph(FT,c, DryMass, FuelMass, check):

    maxSimlen = 3600 * 24 *31 #Maximum 1 monat

    sol, xSol, ySol, FuelDepletedTime = updateIVP(float(FT),float(c),maxSimlen, DryMass, FuelMass)

    main_radiusTime =  np.stack((np.sqrt(xSol[:int(FuelDepletedTime)]**2 + ySol[:int(FuelDepletedTime)]**2), sol.t), axis=-1)

    fig = go.Figure()
    fig.add_trace(go.Line(x=xSol[:int(FuelDepletedTime)], y=ySol[:int(FuelDepletedTime)], name='Spacecraft Trajectory',
        hovertemplate=(
                "<b>X:</b> %{x}<br>"
                "<b>Y:</b> %{y:.2f}<br>"
                "<b>Radius:</b> %{customdata[0]:.2f}<br>"
                "<b>Time:</b> %{customdata[1]:.2f}<extra></extra>"),
         customdata = main_radiusTime))

    if FuelDepletedTime  != 0:
        #Marker point at which fuel runs out
        fig.add_trace(go.Scatter(
            x=[xSol[int(FuelDepletedTime)-1]],  # X-coordinate of Fuel cutoff point
            y=[ySol[int(FuelDepletedTime)-1]], 
            mode='markers',
            marker=dict(size=12, color='red', symbol='star'),
            name='No Fuel/MaxTime'
            ))
        
        #Final orbit
        final_radius =  np.sqrt(xSol[int(FuelDepletedTime):]**2 + ySol[int(FuelDepletedTime):]**2)

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
            customdata=final_radius
            ))
        
        #Anotate with final orbit infomation
        fig.add_annotation(
            x=xSol[int(FuelDepletedTime)-1],  
            y=ySol[int(FuelDepletedTime)-1],  
            text=f"Orbital Period: {calcOrbitalElements(sol.y[0][int(FuelDepletedTime)-1], sol.y[2][int(FuelDepletedTime)-1], sol.y[3][int(FuelDepletedTime)-1])[0]/(3600):.3f}h", 
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
    
    return fig, f"The craft has a dV of {(-c * np.log(DryMass/(DryMass+FuelMass))):.2f} km/s"

if __name__ == '__main__':
    app.run(debug=True)