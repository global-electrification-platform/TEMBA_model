% TEMBA Climate User's Guide
% Camilo Ramirez (expanded from the work from Ioannis Pappis)
% June 22, 2023

This folder contains all the scripts and data necessary to reproduce the
work for the TEMBA Climate model.

## Setup and Installation

To run this analysis you need to install `GLPK` and the `CPLEX` solver

You should also have a `python >=3.6` environment setup with the dependencies listed in the 
`envs/environment.yaml` file

To install the environment you need to install first a distribution of anaconda or miniconda, open
a anaconda promp and run:

```
conda env create -n temba -f envs/environment.yaml
```

Finally, activate the environment with `conda activate temba`.

## Running the TEMBA Climate workflow

The TEMBA Climate model uses two automated workflows to run the analysis and produce all result files.
We use the `snakemake` python package for that. Learn more about snakemake in: [https://snakemake.readthedocs.io/en/stable/index.html](https://snakemake.readthedocs.io/en/stable/index.html)

In the rooth of the `04. TEMBA model` folder, you will find two `snakemake` files:

* `snakefile`: This file contains all the code and rules needed to run the TEMBA Climate model from begining to end. 
It uses `GLPK` to produce the model based on the input excel scenario files and `CPLEX` to solve the model. 
It also runs some postprocessing analysis to clean-up the results and create the results of interest such as 
`Scenario CSVs/scenario/LCOEs.csv` and `Scenario CSVs/scenario/GridEmissionFactors.csv`. 
* `snake_merge_results`: Contains the code and rules to create the `Main results`, divided `Per scenario` and 
`Aggregated` results for each country.

To run any of the workflows, you should first perform a dry run running `snakemake` with the `-n` flag. In a 
conda promp or the command line run `snakemake -s snakefile -n` this will perform the dry run listing all the rules, 
and result files it would create. Inspect this first to be sure it is what you want to run and then run the workflow with 
`snakemake -s snakefile --cores n`. replace `snakefile` for the workflow name you want to run (either `snakefile` or 
`snake_merge_results`) and `n` for the number of cores you wan to use to run the model in parallel. 

## Folder structure

- Input data are stored in `.xlsx` Excel files in the `input_data` folder. The available scenarios are
  * BU: Bottom-up demand
  * Low: Top-down low demand
  * High: Top-down high demand
  * BU_CT_high: Bottom-up demand with high cost of carbon
  * Low_CT_high: Top-down low demand with high cost of carbon
  * High_CT_high: Top-down high demand with high cost of carbon
- A modified OSeMOSYS model file is stored in `model` folder
- Temporary output data is stored in the `output_data` folder
- Final results are stored in the `results` folder (see below for a description of this folder)
- All the scripts are stored in the `scripts` folder
- The conda environment is stored in the `envs` folder

### The Results folder

The results folder contains the results of all TEMBA Climate model runs and is structured in the following way:

- `Main results`: This folder contains the consolidated results of the TEMBA Climate model. Results are presented 
for each country, either aggregated for all scenarios or split per scenario.
  * `Aggreagated`: Contains the results for each country and all scenarios aggregated.
    - *Country*: iso code of the country.
	  * `ElectricityDemand.csv`: lists the electricity demand (GWh) in each year divided per demand type (Current 
	  electrified, New electrified and exports). It also divides the results between scenario (Bottom-up, 
	  Top-down high and Top-down low) and carbon tax/cost (No, High).
	  * `ElectricityDemand.html`: same as above but representing the data in graphical format.
	  * `ElectricityGeneration.csv`: lists the electricity generation (GWh) in each year divided per source type  
	  (fuel types and imports). It also divides the results between scenario (Bottom-up, Top-down high and Top-down low) 
	  and carbon tax/cost (No, High).
	  * `ElectricityGeneration.html`: same as above but representing the data in graphical format.
	  * `InstalledCapacity.csv`: lists the power generation installed capacity (GW) in each year divided per source type  
	  (fuel types). It also divides the results between scenario (Bottom-up, Top-down high and Top-down low) 
	  and carbon tax/cost (No, High).
	  * `InstalledCapacity.html`: same as above but representing the data in graphical format.
  * `Per scenario`: Contains separate results for each scenario and each country.
    - *Scenario*: name of the scenario.
	  * `Capacity`: Contains results related to the yearly power generation installed capacity (GW) of each country. 
	  Results are presented in several text and graphical formats (`CSV`, `HTML`, `PNG`).
	  * `Demand`: Contains results related to the yearly electricity demand (GWh) of each country. 
	  Results are presented in several text and graphical formats (`CSV`, `HTML`, `PNG`).
	  * `Generation`: Contains results related to the yearly power generation (GWh) of each country. 
	  Results are presented in several text and graphical formats (`CSV`, `HTML`, `PNG`).
- `Model files`: This are scenario specific result files from the `CPLEX` solver such as the `scenario.sol` files.
- `Scenario CSVs`: Contains the raw `CSV` results for each scenario. This results are quite desaggregated showing 
results per technology, fuel and year in some cases.
  * *Scenario*: scenario name.
	- `AnnualEmissions.csv`
	- `AnnualFixedOperatingCost.csv`
	- `AnnualTechnologyEmission.csv`
	- `AnnualVariableOperatingCost.csv`
	- `AvgCapacityCost.csv`
	- `CapitalInvestment.csv`
	- `ElectricityDemand.csv`
	- `ElectricityGeneration.csv`
	- `GridEmissionFactors.csv`
	- `InstalledCapacity.csv`
	- `LCOEs.csv`
	- `NewCapacity.csv`
	- `ProductionByTechnologyAnnual.csv`
	- `RateOfActivity.csv`
	- `TotalCapacityAnnual.csv`
	- `UseByTechnologyAnnual.csv`
- `capacity_costs.csv`: Aggregated new installed capacity cost in USD/kW for each scenario.
- `emission_factors.csv`: Aggregated grid emission factors in gCO<sub>2eq</sub>/kWh for each scenario.
- `LCOEs.csv`: Aggregated Levelized Cost of Electricity in USD/kWh for each scenario.
