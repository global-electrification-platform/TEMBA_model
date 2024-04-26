import os
import pandas as pd
import numpy as np
from utils import *

    
## Define scenario parameters
discount_rate = 0.08
end_year = 2030
base_year = 2020


## Reading results    
scenario = snakemake.params.scenario
results_folder = snakemake.params.results_folder

sheets = ['NewCapacity', 'TotalCapacityAnnual', 'NewCapacity',
          'ProductionByTechnologyAnnual', 'UseByTechnologyAnnual', 
          'CapitalInvestment', 'AnnualFixedOperatingCost', 'AnnualVariableOperatingCost', 
          'AnnualEmissions', 'AnnualTechnologyEmission']
dfs = {}
for sheet in sheets:
    dfs[sheet] = pd.read_csv(os.path.join(results_folder, sheet + ".csv"))
    dfs[sheet] = dfs[sheet].loc[(dfs[sheet]['y'] >= base_year) & (dfs[sheet]['y'] <= end_year)]

temba_file = r"input_data\Scenarios\{}.xlsx".format(scenario)

# Reading technology life
operational_life = pd.read_excel(temba_file,
                                 sheet_name='OperationalLife')
operational_life.rename({'TECHNOLOGY': 't', 'VALUE': 'OperationalLife'}, axis=1, inplace=True)

# Reading fuel costs
fuel_costs = pd.read_excel(temba_file,
                                 sheet_name='VariableCost')
fuel_costs.rename({'TECHNOLOGY': 't'}, axis=1, inplace=True)
fuel_costs.replace(0.00001, 0, inplace=True)
fuel_costs = fuel_costs.loc[fuel_costs['MODEOFOPERATION']==1]
fuel_costs.drop('MODEOFOPERATION', axis=1, inplace=True)
fuel_costs = fuel_costs.melt(id_vars=['t']).rename({'variable': 'y', 'value': 'FuelCost'}, axis=1)
fuel_costs['Country'] = fuel_costs['t'].str[:3]
fuel_costs = fuel_costs.loc[(fuel_costs['y'] >= base_year) & (fuel_costs['y'] <= end_year)]

# Reading emission factors
emission_factors = pd.read_excel(temba_file,
                                 sheet_name='EmissionActivityRatio')
emission_factors.rename({'TECHNOLOGY': 't', 'EMISSION': 'e'}, axis=1, inplace=True)
emission_factors.replace(0.00001, 0, inplace=True)
emission_factors.drop('MODEOFOPERATION', axis=1, inplace=True)
emission_factors = emission_factors.melt(id_vars=['t', 'e']).rename({'variable': 'y', 'value': 'EmissionFactor'}, axis=1)
emission_factors['Country'] = emission_factors['t'].str[:3]
emission_factors['e'] = emission_factors['e'].str[3:]
emission_factors = emission_factors.loc[(emission_factors['y'] >= base_year) & (emission_factors['y'] <= end_year)]
emission_factors['t'] = emission_factors['t'].str.replace('LNG', 'GAS')

# Reading emission penalty
emission_penalty = pd.read_excel(temba_file,
                                 sheet_name='EmissionsPenalty')
emission_penalty.rename({'EMISSION': 'e'}, axis=1, inplace=True)
emission_penalty = emission_penalty.melt(id_vars=['e']).rename({'variable': 'y', 'value': 'EmissionPenalty'}, axis=1)
emission_penalty['Country'] = emission_penalty['e'].str[:3]
emission_penalty['e'] = emission_penalty['e'].str[3:]
emission_factors = emission_factors.loc[(emission_factors['e']=='CO2')]
emission_penalty = emission_penalty.loc[(emission_penalty['y'] >= base_year) & (emission_penalty['y'] <= end_year)]
emission_penalty.set_index(['Country', 'y'], inplace=True)

# Get used fuels for power generation

all_techs = dfs['ProductionByTechnologyAnnual'].loc[dfs['ProductionByTechnologyAnnual']['f'].str.contains('EL1|EL3E')]['t'].unique()
all_techs = [x for x in all_techs if (x[3:]!='BACKSTOP')]
gen_techs = [x for x in all_techs if (x[3:5]!='EL' and x[-4:]!='BP00')]
# gen_techs.append(['SOLPUP00X', 'SOLPRP00X'])
gen_techs_df = categorize(dfs['UseByTechnologyAnnual'], gen_techs, ['Country', 'y', 't'])
gen_techs_df.reset_index(inplace=True)

gen_fuels = categorize(dfs['UseByTechnologyAnnual'], gen_techs, ['Country', 'y', 'f'])
gen_fuels.reset_index(inplace=True)
gen_fuels = gen_fuels.loc[~gen_fuels['f'].str.contains('SOL|WIN')]
fuels = gen_fuels['f'].str[3:].unique()

# Extracting emissions and emissions costs

emission_factors.set_index(['Country', 'y', 't'], inplace=True)

emissions = gen_techs_df.copy()
emissions['EmissionFactor'] = emissions.apply(lambda row: get_emission_factor(row['Country'], 
                                                                              row['y'], 
                                                                              row['t'],
                                                                              emission_factors),
                                              axis=1)
emissions['Emissions'] = emissions['EmissionFactor'] * emissions['UseByTechnologyAnnual']
emissions['EmissionPenalty'] = emissions.apply(lambda row: emission_penalty.loc[(row['Country'], row['y']), 'EmissionPenalty'],
                                              axis=1)
emissions['EmissionsCost'] = emissions['EmissionPenalty'] * emissions['Emissions']
emissions_cost = emissions.groupby(['Country', 'y'])[['EmissionsCost']].sum()

# Calculate total O&M costs (fuel and fixed)

om_costs = []
for fuel in fuels:
    costs = get_fuel_costs(dfs['ProductionByTechnologyAnnual'],
                           dfs['AnnualVariableOperatingCost'], 
                           fuel_costs,
                           gen_fuels,
                           fuel)
    costs.name = fuel
    om_costs.append(costs)

om_costs = pd.concat(om_costs, axis=1).fillna(0)
om_costs['VariableCosts'] = om_costs[fuels].sum(axis=1)

fix_costs = categorize(dfs['AnnualFixedOperatingCost'], gen_techs, ['Country', 'y'])
om_costs = pd.concat([om_costs, fix_costs, emissions_cost], axis=1).fillna(0)
om_costs['TotalOMCosts'] = om_costs[['VariableCosts', 'AnnualFixedOperatingCost', 'EmissionsCost']].sum(axis=1)

# Get capital investments and salvage value

dfs['CapitalInvestment'] = dfs['CapitalInvestment'].merge(operational_life, on='t', how='left')
capital_investment = categorize(dfs['CapitalInvestment'], gen_techs, ['Country', 'y', 't'])
capital_investment.reset_index(inplace=True)

capital_investment['Salvage'] = salvage_value(capital_investment['CapitalInvestment'],
                                                        capital_investment['y'],
                                                        base_year,
                                                        end_year,
                                                        capital_investment['OperationalLife'],
                                                        discount_rate)

capital_investment = capital_investment.groupby(['Country', 'y'])[['Salvage', 'CapitalInvestment']].sum()

# Get grid electricity production

production_by_tech = categorize(dfs['ProductionByTechnologyAnnual'], gen_techs, ['Country', 'y'])

# Calculate Country LCOE without imports

data = pd.concat([capital_investment, om_costs['TotalOMCosts'], production_by_tech], axis=1).fillna(0)
data.reset_index(inplace=True)
country_lcoe = lcoe(data,
                    capital_costs='CapitalInvestment', 
                    salvage_value='Salvage', 
                    operational_costs='TotalOMCosts', 
                    electricity_generation='ProductionByTechnologyAnnual', 
                    discount_rate=0.08, 
                    years=data['y'],
                    base_year=2020,
                    groupby=['Country'])

# Calculate import costs

import_cost_factor = 1.5

imports = dfs['ProductionByTechnologyAnnual'].loc[dfs['ProductionByTechnologyAnnual']['t'].str.contains('EL[A-Z]*BP00')].copy()
import_techs = imports['t'].unique()
imports['Country'] = imports['f'].str[:3]
imports['ExpCountry'] = [i[0] if i[0]!=j else i[1] for i, j in zip(imports['t'].str[:-4].str.split('EL'), imports['Country'])]
imports['ImportCost'] = imports.apply(lambda row: get_grid_cost(row['ExpCountry'],
                                                                country_lcoe) * 
                                      row['ProductionByTechnologyAnnual'] *
                                      import_cost_factor, axis=1)
imports_all = imports.groupby(['Country', 'y'])[['ProductionByTechnologyAnnual', 'ImportCost']].sum()
imports_all.rename({'ProductionByTechnologyAnnual': 'Imports'}, axis=1, inplace=True)

total_om_costs = pd.concat([om_costs, imports_all['ImportCost']], axis=1).fillna(0)
total_om_costs['TotalOMCosts'] += total_om_costs['ImportCost']

total_production = pd.concat([production_by_tech, imports_all['Imports']], axis=1).fillna(0)
total_production['ProductionByTechnologyAnnual'] += total_production['Imports']

data = pd.concat([capital_investment, total_om_costs['TotalOMCosts'], 
                  total_production['ProductionByTechnologyAnnual'] * 277.777778], axis=1).fillna(0)

data.reset_index(inplace=True)
final_lcoe = lcoe(data,
                  capital_costs='CapitalInvestment', 
                  salvage_value='Salvage', 
                  operational_costs='TotalOMCosts', 
                  electricity_generation='ProductionByTechnologyAnnual', 
                  discount_rate=0.08, 
                  years=data['y'],
                  base_year=2020,
                  groupby=['Country'])

final_lcoe.rename({'DiscountedCosts': 'Total discounted costs (Million USD)',
                   'DiscountedEnergy': 'Total discounted electricity (GWh)',
                   'lcoe': 'LCOE (USD/kWh)'},
                 axis=1, inplace=True)
final_lcoe.to_csv(snakemake.output.lcoes, index=True)

# Saving installed capacity results

source_names = {'BIO': 'Biomass', 
                'GAS': 'Natural Gas', 
                'NGA': 'Natural Gas',
                'HFO': 'Heavy Fuel Oil', 
                'HYD': 'Hydropower', 
                'LFO': 'Light Fuel Oil', 
                'GEO': 'Geothermal', 
                'SOL': 'Solar', 
                'COA': 'Coal', 
                'WIN': 'Wind',
                'URN': 'Nuclear'}

capacity_by_tech = categorize(dfs['TotalCapacityAnnual'], gen_techs, ['Country', 'y', 't'])
capacity_by_tech.reset_index(inplace=True)

used_fuels = categorize(dfs['UseByTechnologyAnnual'], gen_techs, ['Country', 'y', 't', 'f'])
used_fuels.reset_index(inplace=True)
used_fuels['f'] = used_fuels['f'].str[3:]

capacity_by_tech = pd.merge(capacity_by_tech, 
                            used_fuels.groupby('t').agg({'f': 'first'}).reset_index(), 
                            on=['t'], how='left')

capacity_by_tech.loc[capacity_by_tech['t'].str.contains('HYD'), 'f'] = 'HYD'
capacity_by_tech.loc[capacity_by_tech['t'].str.contains('GEO'), 'f'] = 'GEO'
capacity_by_tech = capacity_by_tech.loc[capacity_by_tech['t']!='BACKSTOP']

dff = capacity_by_tech.groupby(['Country', 'y', 'f'])[['TotalCapacityAnnual']].sum().reset_index()

dff['f_names'] = dff['f'].map(source_names)
dff.rename({'TotalCapacityAnnual': 'Total installed capacity (GW)',
            'f_names': 'Source'},
          axis=1, inplace=True)
dff.to_csv(snakemake.output.capacity, index=False)

## Extract new capacity cost

new_capacity = categorize(dfs['NewCapacity'], gen_techs, ['Country', 'y'])
new_capacity.reset_index(inplace=True)
new_capacity = new_capacity.merge(capital_investment.reset_index(), on=['Country', 'y'])

new_capacity_country = new_capacity.groupby('Country')[['NewCapacity', 'CapitalInvestment']].sum().reset_index()
new_capacity_country['AvgCapacityCost'] = new_capacity_country['CapitalInvestment'] / new_capacity_country['NewCapacity']
new_capacity_country.to_csv("results/Scenario_CSVs/{}/AvgCapacityCost.csv".format(scenario), index=False)

# Saving electricity demand results

df = dfs['ProductionByTechnologyAnnual'].loc[dfs['ProductionByTechnologyAnnual']['f'].str.contains('EL3[A-Z]')].copy()
df = df.loc[df['t'].str.contains('EL2DIP|EL3DIP') | df['f'].str.contains('EL3E')]
df['Country'] = df['f'].str[:3]
df['f'] = df['f'].str[3:]
df = df.groupby(['Country', 'y', 'f'])[['ProductionByTechnologyAnnual']].sum().reset_index()
df['f_names'] = df['f'].map({'EL3E': 'Current electrified',
                             'EL3U': 'New electrified'})

exports = imports[['y', 'ProductionByTechnologyAnnual', 'f', 'ExpCountry']].copy().rename({'ExpCountry': 'Country'}, axis=1)
exports_country = exports.groupby(['Country', 'y', 'f'])[['ProductionByTechnologyAnnual']].sum().reset_index()
exports_country['f_names'] = 'Exports (' + exports_country['f'].str[:3] +')'

dff = pd.concat([df, exports_country])
dff['ProductionByTechnologyAnnual'] *= 277.777778
dff.rename({'ProductionByTechnologyAnnual': 'Total annual demand (GWh)',
            'f_names': 'Demand type'},
          axis=1, inplace=True)
dff.to_csv(snakemake.output.demand, index=False)

# Saving electricity production results

production_by_tech = categorize(dfs['ProductionByTechnologyAnnual'], gen_techs, ['Country', 'y', 't']).reset_index()
production_by_tech['f'] = production_by_tech['t'].str[3:6]
import_countries = imports[['Country', 'y', 't', 'ProductionByTechnologyAnnual', 'f', 'ExpCountry']].copy()
import_countries['f'] = 'Imports (' + import_countries['ExpCountry'] + ')'
import_countries.drop('ExpCountry', axis=1, inplace=True)
df = pd.concat([production_by_tech, import_countries], ignore_index=True)
df.loc[df['t'].str.contains('HYD'), 'f'] = 'HYD'
df.loc[df['t'].str.contains('GEO'), 'f'] = 'GEO'
dff = df.groupby(['Country', 'y', 'f'])[['ProductionByTechnologyAnnual']].sum().reset_index()
dff['f_names'] = dff['f'].map(source_names)
dff.loc[dff['f_names'].isna(), 'f_names'] = dff.loc[dff['f_names'].isna(), 'f']
dff['ProductionByTechnologyAnnual'] *= 277.777778

dff.rename({'ProductionByTechnologyAnnual': 'Total annual production (GWh)',
            'f_names': 'Source'},
          axis=1, inplace=True)
dff.to_csv(snakemake.output.generation, index=False)

# Calculating Grid Emissions Factors

country_emissions = emissions.groupby(['Country'])[['Emissions']].sum()
production = production_by_tech.reset_index().groupby(['Country'])['ProductionByTechnologyAnnual'].sum()
country_emissions['GridEmissionFactor'] = country_emissions['Emissions'] / production

imports_country = imports.groupby(['Country', 'ExpCountry'])[['ProductionByTechnologyAnnual']].sum().reset_index()
imports_country['Emissions'] = imports_country.reset_index().apply(lambda row: get_imports_emissions(row['ExpCountry'], country_emissions) *
                                                          row['ProductionByTechnologyAnnual'], axis=1)
imports_country = imports_country.groupby(['Country'])[['ProductionByTechnologyAnnual', 'Emissions']].sum()
imports_country.rename({'ProductionByTechnologyAnnual': 'Imports'}, axis=1, inplace=True)

production = pd.concat([production, imports_country], axis=1).fillna(0)
production['Total'] = production['ProductionByTechnologyAnnual'] + production['Imports']
country_emissions['ProductionByTechnologyAnnual'] = production['Total'] * 277.777778
country_emissions['Emissions'] += production['Emissions']
country_emissions['GridEmissionFactor'] = (country_emissions['Emissions'] * 1000000) / (production['Total'] * 277.777778)

country_emissions.rename({'ProductionByTechnologyAnnual': 'Total annual production (GWh)',
                          'Emissions': 'Emissions (Million tonCO2eq)',
                          'GridEmissionFactor': 'Grid emission factor (gCO2eq/kWh)'},
                        axis=1, inplace=True)
country_emissions[['Total annual production (GWh)',
                   'Emissions (Million tonCO2eq)',
                   'Grid emission factor (gCO2eq/kWh)']].to_csv(snakemake.output.emission_factors, 
                                                                  index=True)



