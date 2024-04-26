import os
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go

# Reading results    
scenario = snakemake.params.scenario
results_folder = snakemake.params.output_folder

color_map = {'Biomass': 'seagreen',
             'Natural Gas': 'crimson',
             'Heavy Fuel Oil': 'gray',
             'Hydropower': 'darkturquoise',
             'Light Fuel Oil': 'lightgray',
             'Coal': 'rgb(50,50,50)',
             'Geothermal': 'purple',
             'Nuclear': 'lightblue',
             'Solar': 'orange',
             'Wind': 'rosybrown'}

scenarios = {'BU': 'Bottom-up no carbon cost',
             'High': 'Top-down high demand, no carbon cost',
             'Low': 'Top-down low demand, no carbon cost',
             'BU_CT_low': 'Bottom-up low cost of carbon',
             'High_CT_low': 'Top-down high demand, low cost of carbon',
             'Low_CT_low': 'Top-down low demand, low cost of carbon',
             'BU_CT_high': 'Bottom-up high cost of carbon',
             'High_CT_high': 'Top-down high demand, high cost of carbon',
             'Low_CT_high': 'Top-down low demand, high cost of carbon'}

# Installed capacity plots

dff = pd.read_csv(snakemake.input.capacity)

os.makedirs(os.path.join(results_folder, 'Capacity', 'PNG'), exist_ok=True)
os.makedirs(os.path.join(results_folder, 'Capacity', 'HTML'), exist_ok=True)
os.makedirs(os.path.join(results_folder, 'Capacity', 'CSV'), exist_ok=True)
for country in dff['Country'].unique():
    fig = px.bar(dff.loc[dff['Country']==country], 
                 x="y", y="Total installed capacity (GW)", 
                 color="Source", 
                 title=f'Installed Capacity - {country}<br><sup>{scenarios[scenario]} scenario</sup>',
                 template='plotly_white',
                 labels={'y': 'Year'},
                 color_discrete_map=color_map,
                )
    fig.write_image(os.path.join(results_folder, 'Capacity', 'PNG', f'{country}.png'), width=700, height=400, scale=4)
    fig.write_html(os.path.join(results_folder, 'Capacity', 'HTML', f'{country}.html'), full_html=False, include_plotlyjs='cdn')
    dff.loc[dff['Country']==country].to_csv(os.path.join(results_folder, 'Capacity', 'CSV', f'{country}.csv'), index=False)
    
# Creating demand plots

dff = pd.read_csv(snakemake.input.demand)
    
os.makedirs(os.path.join(results_folder, 'Demand', 'PNG'), exist_ok=True)
os.makedirs(os.path.join(results_folder, 'Demand', 'HTML'), exist_ok=True)
os.makedirs(os.path.join(results_folder, 'Demand', 'CSV'), exist_ok=True)
for country in dff['Country'].unique():
    fig = px.bar(dff.loc[dff['Country']==country], 
                 x="y", y="Total annual demand (GWh)", 
                 color="Demand type", 
                 title=f'Electricity demand - {country}<br><sup>{scenarios[scenario]} scenario</sup>',
                 template='plotly_white',
                 labels={'y': 'Year'},
                 color_discrete_map={'Current electrified': 'teal',
                                     'New electrified': 'firebrick'},
                )
    fig.write_image(os.path.join(results_folder, 'Demand', 'PNG', f'{country}.png'), width=700, height=400, scale=4)
    fig.write_html(os.path.join(results_folder, 'Demand', 'HTML', f'{country}.html'), full_html=False, include_plotlyjs='cdn')
    dff.loc[dff['Country']==country].to_csv(os.path.join(results_folder, 'Demand', 'CSV', f'{country}.csv'), index=False)

# Creating production plots

dff = pd.read_csv(snakemake.input.generation)

os.makedirs(os.path.join(results_folder, 'Generation', 'PNG'), exist_ok=True)
os.makedirs(os.path.join(results_folder, 'Generation', 'HTML'), exist_ok=True)
os.makedirs(os.path.join(results_folder, 'Generation', 'CSV'), exist_ok=True)
for country in dff['Country'].unique():
    fig = px.bar(dff.loc[dff['Country']==country], 
                 x="y", y="Total annual production (GWh)", 
                 color="Source", 
                 title=f'Electricity generation - {country}<br><sup>{scenarios[scenario]} scenario</sup>',
                 template='plotly_white',
                 labels={'y': 'Year'},
                 color_discrete_map=color_map,
                )
    fig.write_image(os.path.join(results_folder, 'Generation', 'PNG', f'{country}.png'), width=700, height=400, scale=4)
    fig.write_html(os.path.join(results_folder, 'Generation', 'HTML', f'{country}.html'), full_html=False, include_plotlyjs='cdn')
    dff.loc[dff['Country']==country].to_csv(os.path.join(results_folder, 'Generation', 'CSV', f'{country}.csv'), index=False)