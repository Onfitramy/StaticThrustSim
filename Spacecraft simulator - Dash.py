import numpy as np
import scipy.integrate

from dash import Dash, html, dcc, callback, Output, Input
import plotly.express as px
import pandas as pd

simlen = 2000  #Simulation lenght in h

𝛼 = 90   #control variable in °
FT = 10  #Thrust in N
c = 20000 #exhaust velocity in m/s

r0 = 300   #Initial orbit height in km
v0 = 0   #Initial orbital velocity in km/s (set zero to use Eccentricity)
E = 0    #Eccentricity
m0 = 1000   #Initial mass in kg
µ = 3.986e5 #Gtravitational parameter in km3/s2

if v0 == 0:
    v0 = np.sqrt(µ/r0) #For a circular orbit
    ai = r0/(1-E) #Semi-major axis
    v0 = (2*µ/r0 - µ/ai)**0.5

def motionEQ(t, y, FT, µ, 𝛼, c): #y = [r, φ, pr, ρφ, m]
    r, φ, pr, ρφ, m = y
    r_der = pr
    φ_der = ρφ
    pr_der = r * (ρφ**2) - µ/(r**2) + (FT/m)*np.cos(𝛼)
    ρφ_der =  (1/r) * ((FT/m)*np.sin(𝛼) - 2 * pr * ρφ)
    m_der = - FT/c

    return [r_der, φ_der, pr_der, ρφ_der, m_der]

sol = scipy.integrate.solve_ivp(motionEQ, [0, simlen], [r0, 0, 0, v0/r0, m0], args=(FT, µ, 𝛼, c), t_eval= np.linspace(0, simlen, simlen*5))

xSol = sol.y[0] * np.cos(sol.y[1])
ySol = sol.y[0] * np.sin(sol.y[1])

#Beginning of App code
app = Dash()

app.layout = [
    html.H1(children='Spacecraft simulator', style={'textAlign':'center'}),
    html.Br(),
    html.Label('Thrust[N]'),
    dcc.Slider(
            min=0,
            max=10,
            marks={i: str(i) for i in range(1, 9)},
            value=2,
            id='Thrust-setting',
            updatemode='drag',),
    
    html.Label('Exhaust Velocity[m/s]'),
    dcc.Slider(
            min=100,
            max=100000,
            value=20000,
            id='Exhaust-setting',
            updatemode='drag',),

    html.Label('Simulation lenght'),
    dcc.Slider(
            min=100,
            max=10000,
            marks={i: str(i) for i in range(1, 9)},
            value=2000,
            id='Simlen-setting',
            updatemode='drag',),
    dcc.Graph(id='graph-content')
]

@callback(
    Output('graph-content', 'figure'),
    Input('Thrust-setting', 'value'),
    Input('Exhaust-setting', 'value'),
    Input('Simlen-setting', 'value')
)

def update_graph(FT,c,simlen):

    xSol, ySol = updateIVP(float(FT),float(c),simlen)

    fig = px.line(x=xSol, y=ySol, title='Spacecraft Trajectory')

    fig.update_layout(
        xaxis=dict(scaleanchor="y"),  # Link the scale of x-axis to y-axis
        yaxis=dict(),
        height=1000
    )
    
    return fig

def updateIVP(FT, c, simlen):
    sol = scipy.integrate.solve_ivp(motionEQ, [0, simlen], [r0, 0, 0, v0/r0, m0], args=(FT, µ, 𝛼, c), t_eval= np.linspace(0, simlen, simlen*5))

    xSol = sol.y[0] * np.cos(sol.y[1])
    ySol = sol.y[0] * np.sin(sol.y[1])

    return xSol, ySol

if __name__ == '__main__':
    app.run(debug=True)