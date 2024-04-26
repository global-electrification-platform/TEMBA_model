import os
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objs as go


def discount_factor(discount_rate, years):
    df = (1 + discount_rate) ** - years
    return df
    
def used_life(years, end_year, tech_life):
    return (end_year - years + 1) % tech_life

def discount_value(value, discount_rate, years, base_year):
    df = discount_factor(discount_rate, years - base_year)
    discounted_value = value * df
    return discounted_value
    
def salvage_value(capital_costs, years, base_year, end_year, tech_life, discount_rate):
    tech_used_life = used_life(years, end_year, tech_life)
    salvage = capital_costs * (1 - tech_used_life / tech_life)
    return salvage
    
def lcoe(data, capital_costs, salvage_value, operational_costs, 
         electricity_generation, discount_rate, years, base_year, groupby):
    data = data.copy()
    data['DiscountedCosts'] = discount_value(data[capital_costs] - data[salvage_value] + data[operational_costs], 
                                             discount_rate, 
                                             years, 
                                             base_year)
    data['DiscountedEnergy'] = discount_value(data[electricity_generation], 
                                              discount_rate, 
                                              years, 
                                              base_year)
    lcoe = data.groupby(groupby)[['DiscountedCosts', 'DiscountedEnergy']].sum()
    lcoe['lcoe'] = lcoe['DiscountedCosts'] /  lcoe['DiscountedEnergy']
    return lcoe
    
def categorize(df, categories, groupby):
    df = df.copy()
    rows = df['t'].isin(categories)
    df = df.loc[rows]
    df['Country'] = df['t'].str[:3].copy()
    return df.groupby(groupby).sum()
    
def get_fuel_costs(prod_by_tech, var_cost, fuel_costs, gen_fuels, fuel):
    if fuel in ['HFO', 'LFO']:
        co1_prod_techs = prod_by_tech.loc[prod_by_tech['t'].str.contains('CO1IMP|CO1PRP'), 't'].unique()
        co1_amount = categorize(prod_by_tech, co1_prod_techs, ['Country', 'y'])

        co1_costs_techs = var_cost.loc[var_cost['t'].str.contains('CO1IMP|CO1PRP'), 't'].unique()
        co1_costs = categorize(var_cost, co1_costs_techs, ['Country', 'y'])
        co1_unit_cost = co1_costs['AnnualVariableOperatingCost'] / co1_amount['ProductionByTechnologyAnnual']
        co1_unit_cost.name = 'ProdCost'

        fuel_techs = prod_by_tech.loc[prod_by_tech['f'].str.contains(fuel), 't'].unique()
        fuel_use = categorize(prod_by_tech, fuel_techs, ['t', 'Country', 'y', 'f']).reset_index()
        fuel_use = fuel_use.loc[fuel_use['f'].str.contains(fuel)]
        fuel_local = fuel_use.loc[fuel_use['t'].str.contains('CO1R1P00X')].groupby(['Country', 'y'])['ProductionByTechnologyAnnual'].sum()
        fuel_local.name = 'ProdAmount'
        fuel_imp = fuel_use.loc[fuel_use['t'].str.contains(f'{fuel}IMP00X')].groupby(['Country', 'y'])['ProductionByTechnologyAnnual'].sum()
        fuel_imp.name = 'ImpAmount'

        fuel_imp_techs = fuel_costs.loc[fuel_costs['t'].str.contains(f'{fuel}IMP'), 't'].unique()
        fuel_imp_cost = categorize(fuel_costs, fuel_imp_techs, ['Country', 'y'])['FuelCost']
        fuel_imp_cost.name = 'ImpCost'

        gen_fuel = gen_fuels.loc[gen_fuels['f'].str.contains(fuel)].set_index(['Country', 'y'])['UseByTechnologyAnnual']
        gen_fuel.name = 'FuelUse'

        fuel_df = pd.concat([fuel_local, co1_unit_cost, fuel_imp, fuel_imp_cost, gen_fuel], axis=1).fillna(0)
        fuel_cost = (fuel_df['ProdAmount'] * fuel_df['ProdCost'] + fuel_df['ImpAmount'] * fuel_df['ImpCost']) / (fuel_df['ProdAmount'] + fuel_df['ImpAmount'])
        return fuel_cost * fuel_df['FuelUse']
    else:
        fuel_techs = prod_by_tech.loc[prod_by_tech['f'].str.fullmatch(f'[A-Z]*{fuel}'), 't'].unique()
        fuel_use = categorize(prod_by_tech, fuel_techs, ['t', 'Country', 'y', 'f']).reset_index()
        fuel_use = fuel_use.loc[fuel_use['f'].str.contains(fuel)]
        fuel_use_all = fuel_use.groupby(['Country', 'y'])['ProductionByTechnologyAnnual'].sum()
        fuel_use_all.name = 'TotalFuel'
        fuel_use = fuel_use.drop('f', axis=1).set_index(['t', 'Country', 'y'])
        fuel_use.rename({'ProductionByTechnologyAnnual': 'FuelProd'}, axis=1, inplace=True)
        
        fuel_costs = categorize(fuel_costs, fuel_techs, ['t', 'Country', 'y'])
        
        fuel_df = pd.concat([fuel_use, fuel_costs], axis=1).fillna(0).reset_index()
        fuel_df['TotalCost'] = fuel_df['FuelProd'] * fuel_df['FuelCost']
        
        unit_cost = fuel_df.groupby(['Country', 'y'])['TotalCost'].sum() / fuel_use_all
        
        gen_fuel = gen_fuels.loc[gen_fuels['f'].str.contains(f'[A-Z]*{fuel}')].set_index(['Country', 'y'])['UseByTechnologyAnnual']
        gen_fuel.name = 'FuelUse'
        
        return (unit_cost * gen_fuel).fillna(0)
        
def get_grid_cost(exp_country, lcoes):
    if exp_country not in lcoes.index:
        cost = lcoes['lcoe'].mean()
    else:
        cost = lcoes.loc[exp_country, 'lcoe']
    return cost
    
def get_emission_factor(country, y, tech, emission_factors):
    return emission_factors.loc[(country, y, tech), 'EmissionFactor']
        
def get_imports_emissions(exp_country, country_emissions):
    if exp_country not in country_emissions.index:
        return country_emissions['GridEmissionFactor'].mean()
    else:
        return country_emissions.loc[exp_country]['GridEmissionFactor']
        
def compare_bar_plot(df, value_vars, labels, color_pallete=px.colors.qualitative.Bold):
    df = df.copy()
    df['Country'] = df['Country'] + ' -'
    dff = df.melt(id_vars=['Country', 'Diff', 'Group'], value_vars=value_vars)
    fig = px.bar(dff, 
                 y="Country", x="value", 
                 orientation='h',
                 color='variable',
                 color_discrete_sequence=color_pallete,
                 labels=labels,
                 facet_col='Group',
                 facet_col_wrap=2,
                 facet_col_spacing=0.1,
                 category_orders={'variable': value_vars})

    fig.update_yaxes(matches=None, showticklabels=True)

    fig.update_layout(barmode='group', 
                      template='plotly_white', 
                      legend=dict(orientation="h", yanchor='top', y=-0.2))

    i = 44
    for axis in fig.layout:
        if type(fig.layout[axis]) == go.layout.XAxis:
            fig.layout[axis].title.text = ''
        if type(fig.layout[axis]) == go.layout.YAxis:
            fig.layout[axis].categoryorder = 'array'
            fig.layout[axis].categoryarray = [x for x in df.iloc[i-22:i]['Country'].iloc[::-1]]
            i -= 22

    fig.for_each_annotation(lambda a: a.update(text=''))

    max_value = dff.groupby(['Country'])['value'].max()
    max_value += max_value.max() * 0.07
    dff.set_index('Country', inplace=True)
    dff['max'] = max_value
    dff.reset_index(inplace=True)
    dff_diff = dff.loc[dff['variable']=='TEMBA'].set_index('Group')
    dff_diff = dff_diff.loc[(dff_diff['Diff']!='nan%') & (dff_diff['Diff']!='inf%')]
   
    fig.add_trace(go.Scatter(y=dff_diff.loc[1]["Country"], x=dff_diff.loc[1]["max"], 
                             text=dff_diff.loc[1]["Diff"], mode='text', textfont=dict(size=9, color='gray'),
                             showlegend=False),
                  row=1, col=1)
    fig.add_trace(go.Scatter(y=dff_diff.loc[2]["Country"], x=dff_diff.loc[2]["max"], 
                             text=dff_diff.loc[2]["Diff"], mode='text', textfont=dict(size=9, color='gray'),
                             showlegend=False),
                  row=1, col=2)

    fig.add_annotation(x=0.5, y=-0.15, text=labels['value'], xanchor='center',
                       font=dict(size=14), xref='paper', yref='paper', showarrow=False)
    return fig