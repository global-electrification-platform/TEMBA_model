import pandas as pd

scenarios = snakemake.params.scenarios

df_lcoe = pd.DataFrame()
df_emissions = pd.DataFrame()
df_capacity = pd.DataFrame()

for scenario in scenarios:
    lcoe = pd.read_csv(f"results/Scenario_CSVs/{scenario}/LCOEs.csv")
    emission_factors = pd.read_csv(f"results/Scenario_CSVs/{scenario}/GridEmissionFactors.csv")
    capacity_costs = pd.read_csv(f"results/Scenario_CSVs/{scenario}/AvgCapacityCost.csv")
    
    lcoe.set_index('Country', inplace=True)
    emission_factors.set_index('Country', inplace=True)
    capacity_costs.set_index('Country', inplace=True)
    
    lcoe.rename(columns={'LCOE (USD/kWh)': scenario}, inplace=True)
    emission_factors.rename(columns={'Grid emission factor (gCO2eq/kWh)': scenario}, inplace=True)
    capacity_costs.rename(columns={'AvgCapacityCost': scenario}, inplace=True)
    
    df_lcoe = pd.concat([df_lcoe, lcoe[scenario]], axis=1)
    df_emissions = pd.concat([df_emissions, emission_factors[scenario]], axis=1)
    df_capacity = pd.concat([df_capacity, capacity_costs[scenario]], axis=1)

df_lcoe.index.name = 'Country'
df_emissions.index.name = 'Country'
df_capacity.index.name = 'Country'

df_lcoe.to_csv(snakemake.output.lcoes)
df_emissions.to_csv(snakemake.output.emission_factors)
df_capacity.to_csv(snakemake.output.capacity_costs)