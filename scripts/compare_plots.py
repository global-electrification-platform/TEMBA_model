import os
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go

# Reading results    
scenario_names = snakemake.params.scenarios

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

scenarios = {'BU': 'Bottom-up', 
             'High': 'Top-down high', 
             'Low': 'Top-down low',
             'BU_CT_high': 'Bottom-up', 
             'High_CT_high': 'Top-down high', 
             'Low_CT_high': 'Top-down low'}

carbon_tax = {'BU': 'No',
             'High': 'No',
             'Low': 'No',
             'BU_CT_low': 'Low',
             'High_CT_low': 'Low',
             'Low_CT_low': 'Low',
             'BU_CT_high': 'High',
             'High_CT_high': 'High',
             'Low_CT_high': 'High'}

dfs = []
for scenario in scenario_names:
    df = pd.read_csv("results/Scenario_CSVs/{}/InstalledCapacity.csv".format(scenario))
    df['Scenario'] = scenarios[scenario.split('_')[0]]
    df['Carbon Price'] = carbon_tax[scenario]
    dfs.append(df)
df_capacity = pd.concat(dfs, axis=0)

dfs = []
for scenario in scenario_names:
    df = pd.read_csv("results/Scenario_CSVs/{}/ElectricityDemand.csv".format(scenario))
    df['Scenario'] = scenarios[scenario.split('_')[0]]
    df['Carbon Price'] = carbon_tax[scenario]
    dfs.append(df)
df_demand = pd.concat(dfs, axis=0)

dfs = []
for scenario in scenario_names:
    df = pd.read_csv("results/Scenario_CSVs/{}/ElectricityGeneration.csv".format(scenario))
    df['Scenario'] = scenarios[scenario.split('_')[0]]
    df['Carbon Price'] = carbon_tax[scenario]
    dfs.append(df)
df_production = pd.concat(dfs, axis=0)

for country in df['Country'].unique():
    folder = os.path.join('results', 'Main_results', 'Aggregated', country)
    os.makedirs(folder, exist_ok=True)
    
    # Installed capacity graphs
    dff = df_capacity.loc[df_capacity['Country']==country]
    fig_capacity = px.bar(dff, 
                          x="y", y="Total installed capacity (GW)", 
                          color="Source",
                          template='plotly_white',
                          labels={'y': 'Year',
                                  "Total installed capacity (GW)": "Installed capacity (GW)"},
                          color_discrete_map=color_map,
                          facet_col='Scenario',
                          facet_row='Carbon Price',
                          category_orders={"Carbon Price": ["No", "High"]},
                          width=900, height=600
                         )
    fig_capacity.for_each_annotation(lambda a: a.update(text=a.text.replace("=", ': ')))
    fig_capacity.update_layout(margin=dict(t=16, r=0, b=0, l=0))
    fig_capacity.write_html(os.path.join(folder, 'InstalledCapacity.html'),
                            full_html=False, include_plotlyjs='cdn')
    dff.to_csv(os.path.join(folder, 'InstalledCapacity.csv'), index=False)
    
    # Electricity demand graphs
    dff = df_demand.loc[df_demand['Country']==country]
    fig_demand = px.bar(dff, 
                        x="y", y="Total annual demand (GWh)", 
                        color="Demand type", 
                        template='plotly_white',
                        labels={'y': 'Year',
                                'Total annual demand (GWh)': "Annual demand (GWh)"},
                        color_discrete_map={'Current electrified': 'teal',
                                            'New electrified': 'firebrick'},
                        facet_col='Scenario',
                        facet_row='Carbon Price',
                        category_orders={"Carbon Price": ["No", "High"]},
                        width=900, height=600
                       )
    fig_demand.for_each_annotation(lambda a: a.update(text=a.text.replace("=", ': ')))
    fig_demand.update_layout(margin=dict(t=16, r=0, b=0, l=0))
    fig_demand.write_html(os.path.join(folder, 'ElectricityDemand.html'),
                          full_html=False, include_plotlyjs='cdn')
    dff.to_csv(os.path.join(folder, 'ElectricityDemand.csv'), index=False)
    
    # electricity generation graphs
    dff = df_production.loc[df_production['Country']==country]
    fig_production = px.bar(dff, 
                            x="y", y="Total annual production (GWh)", 
                            color="Source", 
                            template='plotly_white',
                            labels={'y': 'Year',
                                    "Total annual production (GWh)": "Annual production (GWh)"},
                            color_discrete_map=color_map,
                            facet_col='Scenario',
                            facet_row='Carbon Price',
                            category_orders={"Carbon Price": ["No", "High"]},
                            width=900, height=600
                           )
    fig_production.for_each_annotation(lambda a: a.update(text=a.text.replace("=", ': ')))
    fig_production.update_layout(margin=dict(t=16, r=0, b=0, l=0))
    fig_production.write_html(os.path.join(folder, 'ElectricityGeneration.html'),
                              full_html=False, include_plotlyjs='cdn')
    dff.to_csv(os.path.join(folder, 'ElectricityGeneration.csv'), index=False)
   
# Sub Sahara Africa Installed capacity graph
countries = pd.read_csv(snakemake.input.ssa_countries)
dff = df_capacity.loc[df_capacity['Country'].isin(countries['Country'])].copy()
dff = dff.groupby(['Scenario', 'Carbon Price', 'Source', 'y'])[['Total installed capacity (GW)']].sum()

fig_capacity = px.bar(dff.reset_index(), 
                      x="y", y="Total installed capacity (GW)", 
                      color="Source",
                      template='plotly_white',
                      labels={'y': 'Year'},
                      color_discrete_map=color_map,
                      facet_col='Scenario',
                      facet_row='Carbon Price',
                      category_orders={"Carbon Price": ["No", "High"]},
                      width=900, height=600
                     )
fig_capacity.for_each_annotation(lambda a: a.update(text=a.text.replace("=", ': ')))
fig_capacity.update_layout(margin=dict(t=16, r=0, b=0, l=0))
fig_capacity.write_image(snakemake.output.capacity_plot,
                         scale=4, width=900, height=600)
                         
# Sub Saharan Africa demand plot
dff = df_demand.loc[df_demand['Country'].isin(countries['Country'])].copy()
dff = dff.groupby(['Scenario', 'Carbon Price', 'Demand type', 'y'])[['Total annual demand (GWh)']].sum()
dff.reset_index(inplace=True)
dff = dff.loc[dff['Demand type'].isin(['Current electrified', 'New electrified'])]
fig_demand = px.bar(dff, 
                    x="y", y="Total annual demand (GWh)", 
                    color="Demand type", 
                    template='plotly_white',
                    labels={'y': 'Year',
                            'Total annual demand (GWh)': "Annual demand (GWh)"},
                    color_discrete_map={'Current electrified': 'teal',
                                        'New electrified': 'firebrick'},
                    facet_col='Scenario',
                    facet_row='Carbon Price',
                    category_orders={"Carbon Price": ["No", "High"]},
                    width=900, height=600
                   )
fig_demand.for_each_annotation(lambda a: a.update(text=a.text.replace("=", ': ')))
fig_demand.update_layout(margin=dict(t=16, r=0, b=0, l=0))
fig_demand.write_image(snakemake.output.demand_plot,
                         scale=4, width=900, height=600)
                         
# Sub Saharan Africa generation plot
dff = df_production.loc[df_production['Country'].isin(countries['Country'])].copy()
dff = dff.groupby(['Scenario', 'Carbon Price', 'Source', 'y'])[['Total annual production (GWh)']].sum()
dff.reset_index(inplace=True)
dff = dff.loc[~dff['Source'].str.contains('Imports|Exports')]
fig_production = px.bar(dff, 
                        x="y", y="Total annual production (GWh)", 
                        color="Source", 
                        template='plotly_white',
                        labels={'y': 'Year',
                                "Total annual production (GWh)": "Annual production (GWh)"},
                        color_discrete_map=color_map,
                        facet_col='Scenario',
                        facet_row='Carbon Price',
                        category_orders={"Carbon Price": ["No", "High"]},
                        width=900, height=600
                       )
fig_production.for_each_annotation(lambda a: a.update(text=a.text.replace("=", ': ')))
fig_production.update_layout(margin=dict(t=16, r=0, b=0, l=0))
fig_production.write_html(os.path.join(folder, 'ElectricityGeneration.html'),
                          full_html=False, include_plotlyjs='cdn')
fig_production.write_image(snakemake.output.generation_plot,
                         scale=4, width=900, height=600)

# Sub Sahara Africa LCOEs plot
dfs = []
scenario_names = ['BU', 'High', 'Low']
for scenario in scenario_names:
    df = pd.read_csv("results/Scenario_CSVs/{}/LCOEs.csv".format(scenario))
    df['Scenario'] = scenarios[scenario]
    dfs.append(df)
    
df = pd.concat(dfs, axis=0)
df = df.loc[df['Country'].isin(countries['Country'])]
df['Country'] = df['Country'] + ' -'
lenght = len(countries['Country'])
for i, country in enumerate(countries['Country']):
    if i < int(lenght / 2):
        df.loc[df['Country']==country + ' -', 'group'] = 1
    else:
        df.loc[df['Country']==country + ' -', 'group'] = 2
        
dfs_ct = []
scenario_names = ['BU_CT_high', 'High_CT_high', 'Low_CT_high']
color_pallete=px.colors.qualitative.Bold
for scenario in scenario_names:
    df_ct = pd.read_csv("results/Scenario_CSVs/{}/LCOEs.csv".format(scenario))
    df_ct['Scenario'] = scenarios[scenario]
    dfs_ct.append(df_ct)
    
df_ct = pd.concat(dfs_ct, axis=0)
df_ct = df_ct.loc[df_ct['Country'].isin(countries['Country'])]
df_ct['Country'] = df_ct['Country'] + ' -'
lenght = len(countries)
for i, country in enumerate(countries['Country']):
    if i < int(lenght / 2):
        df_ct.loc[df_ct['Country']==country + ' -', 'group'] = 1
    else:
        df_ct.loc[df_ct['Country']==country + ' -', 'group'] = 2

color_pallete=px.colors.qualitative.Bold

fig = px.bar(df, 
             y="Country", x="LCOE (USD/kWh)", 
             orientation='h',
             color='Scenario',
             color_discrete_sequence=color_pallete,
             facet_col='group',
             facet_col_spacing=0.1,
             width=700, height=700,
             labels={'Scenario': 'No carbon price:<br>High carbon price:'}
            )
fig.add_traces(px.box(df_ct, y="Country", x="LCOE (USD/kWh)", color='Scenario',
               color_discrete_sequence=color_pallete,
               facet_col='group',).data)


fig.for_each_annotation(lambda a: a.update(text=''))
fig.update_yaxes(matches=None, showticklabels=True, autorange="reversed")
fig.update_layout(barmode='group', 
                  template='plotly_white',
                  legend=dict(orientation="h"),
                  boxmode="group",
                  bargroupgap=0,
                  boxgroupgap=0,
                  bargap=0.2,
                  boxgap=0.2)
fig.update_traces(
    selector=dict(type="box"), # update only boxes
    boxpoints="all", # show points
    pointpos=0, # centered
    jitter=0, # no jitter
    line_color="rgba(255,255,255,0)", # hide box lines
    fillcolor="rgba(255,255,255,0)", # hide box fill
    marker=dict(size=8, opacity=0.8, symbol='circle',
                line=dict(width=1,
                color='black'))
)

fig.write_image(snakemake.output.lcoes_plot, scale=4)

# Sub Sahara Africa Emission Factors plot
dfs = []
scenario_names = ['BU', 'High', 'Low']
for scenario in scenario_names:
    df = pd.read_csv("results/Scenario_CSVs/{}/GridEmissionFactors.csv".format(scenario))
    df['Scenario'] = scenarios[scenario]
    dfs.append(df)
    
df = pd.concat(dfs, axis=0)
df = df.loc[df['Country'].isin(countries['Country'])]
df['Country'] = df['Country'] + ' -'
lenght = len(countries['Country'])
for i, country in enumerate(countries['Country']):
    if i < int(lenght / 2):
        df.loc[df['Country']==country + ' -', 'group'] = 1
    else:
        df.loc[df['Country']==country + ' -', 'group'] = 2
        
dfs_ct = []
scenario_names = ['BU_CT_high', 'High_CT_high', 'Low_CT_high']
color_pallete=px.colors.qualitative.Bold
for scenario in scenario_names:
    df_ct = pd.read_csv("results/Scenario_CSVs/{}/GridEmissionFactors.csv".format(scenario))
    df_ct['Scenario'] = scenarios[scenario]
    dfs_ct.append(df_ct)
    
df_ct = pd.concat(dfs_ct, axis=0)
df_ct = df_ct.loc[df_ct['Country'].isin(countries['Country'])]
df_ct['Country'] = df_ct['Country'] + ' -'
lenght = len(countries['Country'])
for i, country in enumerate(countries['Country']):
    if i < int(lenght / 2):
        df_ct.loc[df_ct['Country']==country + ' -', 'group'] = 1
    else:
        df_ct.loc[df_ct['Country']==country + ' -', 'group'] = 2

color_pallete=px.colors.qualitative.Bold

fig = px.bar(df, 
             y="Country", x="Grid emission factor (gCO2eq/kWh)", 
             orientation='h',
             color='Scenario',
             color_discrete_sequence=color_pallete,
             facet_col='group',
             facet_col_spacing=0.1,
             width=700, height=700,
             labels={'Scenario': 'No carbon price:<br>High carbon price:'}
            )
fig.add_traces(px.box(df_ct, y="Country", x="Grid emission factor (gCO2eq/kWh)", color='Scenario',
               color_discrete_sequence=color_pallete,
               facet_col='group',).data)


fig.for_each_annotation(lambda a: a.update(text=''))
fig.update_yaxes(matches=None, showticklabels=True, autorange="reversed")
fig.update_layout(barmode='group', 
                  template='plotly_white',
                  legend=dict(orientation="h"),
                  boxmode="group",
                  bargroupgap=0,
                  boxgroupgap=0,
                  bargap=0.2,
                  boxgap=0.2)
fig.update_traces(
    selector=dict(type="box"), # update only boxes
    boxpoints="all", # show points
    pointpos=0, # centered
    jitter=0, # no jitter
    line_color="rgba(255,255,255,0)", # hide box lines
    fillcolor="rgba(255,255,255,0)", # hide box fill
    marker=dict(size=8, opacity=0.8, symbol='circle',
                line=dict(width=1,
                color='black'))
)

fig.write_image(snakemake.output.emission_factors_plot, scale=4)