import sys, os
import pandas as pd
import numpy as np


def main(demand_file: str, temba_file: str, iso_codes_path: str):

    df = pd.read_csv(demand_file)

    temba_input = pd.read_excel(temba_file, sheet_name='SpecifiedAnnualDemand',
                                engine='openpyxl')

    df.drop('Unnamed: 0', axis=1, inplace=True)
    set_iso_codes(file_path=iso_codes_path, df=df)

    init_year = 2015
    base_year = 2020
    first_step_year = 2025
    second_step_year = 2030
    final_year = 2040
    
    countries = [x[:3] for x in temba_input['FUEL']]
    countries = list(set(countries).intersection(df['ISO3']))
    df.set_index('ISO3', inplace=True)
    df = df.loc[countries].copy()
    temba_input.set_index('FUEL', inplace=True)
    
    df['GridDemandAll'] = extract_demand(df=df,
                                         temba_input=temba_input, 
                                         init_year=init_year, 
                                         final_year=final_year, 
                                         base_year=base_year,
                                         first_step_year=first_step_year,
                                         second_step_year=second_step_year)                                 
                                         
    df['FUEL'] = df.index + 'EL3U'
    for i, row in df.iterrows():
        temba_input.loc[row['FUEL']] = row['GridDemandAll']
       
    with pd.ExcelWriter(temba_file, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        temba_input.to_excel(writer, 'SpecifiedAnnualDemand')
    
 
def extract_demand(df: pd.DataFrame, temba_input: pd.DataFrame, 
                   init_year: int, final_year: int, base_year: int, 
                   first_step_year: int, second_step_year: int):
    demands = []
    for i, row in df.iterrows():
        growth_rate = max(0,
                          carg(temba_input.loc[f"{i}EL3E", final_year], 
                               temba_input.loc[f"{i}EL3E", second_step_year], 
                               final_year - second_step_year))
        zero_period = np.zeros(base_year - init_year)
        first_period = np.linspace(0, 
                                   row['GridDemand2025'], 
                                   first_step_year - base_year + 1, 
                                   endpoint=False)
        second_period = np.linspace(row['GridDemand2025'], 
                                    row['GridDemand2030'], 
                                    second_step_year - first_step_year)
        final_period = np.array([row['GridDemand2030'] * (growth_rate + 1) ** (year) for year in range(1, final_year - second_step_year + 1)])
        complete_period = np.concatenate([zero_period, first_period, second_period, final_period])
        demands.append(complete_period)
        
    return demands
    
def set_iso_codes(file_path: str, df: pd.DataFrame):
    iso_codes = pd.read_csv(file_path)
    iso_codes.set_index('Alpha-2', inplace=True)
    df['CountryCode'] = df['CountryCode'].str.upper()
    df.rename({'CountryCode': 'ISO2'}, inplace=True, axis=1)
    df.set_index('ISO2', inplace=True)
    df['ISO3'] = iso_codes['Alpha-3']
    df.reset_index(inplace=True)
    
    
def carg(v_final, v_begin, years):
    """
    Calculates the Compound Annual Growth Rate (CARG)
    """
    return (v_final / v_begin) ** (1 / years) - 1    


if __name__ == '__main__':
    if len(sys.argv) != 4:
        msg = "Usage: python {} <filepath>"
        print(msg.format(sys.argv[0]))
        sys.exit(1)
    else:
        demand_file = sys.argv[1]
        temba_file = sys.argv[2] 
        iso_codes_path = sys.argv[3]
        try:
            main(demand_file=demand_file, 
                 temba_file=temba_file, 
                 iso_codes_path=iso_codes_path)
        except:
            sys.exit(1)