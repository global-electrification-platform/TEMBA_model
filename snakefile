MODELRUNS = ["BU", "BU_CT_high", "Low", "Low_CT_high", "High", "High_CT_high"]

rule all:
    input: 
        expand("results/Main_results/Per_scenario/{model_run}/Generation/CSV/ZMB.csv", model_run=MODELRUNS)

# rule extract_demand:
    # input:
        # demand = r"input_data/GEP_Demand/{model_run}_Demand.csv",
        # iso_codes = r"input_data/ISO_codes.csv"
    # params:
        # temba = "input_data/{model_run}.xlsx"
    # output: 
        # "input_data/{model_run}.xlsx"
    # shell:
        # "python scripts/extract_demand.py {input.demand} {params.temba} {input.iso_codes}"
        

rule generate_model_file:
    input: 
        "input_data/Scenarios/{model_run}.xlsx"
    params:
        "output_data/{model_run}_CSVs"
    output: 
        "output_data/{model_run}.txt"
    threads: 
        1
    shell:
        "python scripts/excel_to_osemosys.py {input} {output} {params}"

rule modify_model_file:
    input:  
        "output_data/{model_run}.txt"
    output: 
        "output_data/{model_run}_modex.txt"
    threads: 
        1
    shell:
        "python scripts/CBC_results_AS_MODEX.py {input} && cat {input} > {output}"

rule generate_lp_file:
    input: 
        "output_data/{model_run}_modex.txt"
    output: 
        "output_data/{model_run}.lp.gz"
    log: 
        "output_data/glpsol_{model_run}.log"
    threads: 
        1
    shell:
        "glpsol -m model/Temba_0406_modex.txt -d {input} --wlp {output} --check --log {log}"

rule solve_lp:
    input: 
        "output_data/{model_run}.lp.gz"
    output: 
        "output_data/{model_run}.sol"
    log: 
        "output_data/gurobi_{model_run}.log"
    threads: 
        1
    shell:
        'cplex -c "read {input}" "optimize" "write {output}"'

rule remove_zero_values:
    input: "output_data/{model_run}.sol"
    output: "results/Model_files/{model_run}.sol"
    shell:
        "sed '/ * 0$/d' {input} > {output}"
        
rule transform_results:
    input: "results/Model_files/{model_run}.sol"
    output: "results/Model_files/{model_run}_transform.txt"
    shell:
        "Python scripts/transform_31072013.py {input} {output}"
        
rule sort_results:
    input: "results/Model_files/{model_run}_transform.txt"
    output: "results/Model_files/{model_run}_sorted.txt"
    shell:
        "sort < {input} > {output}"
        
rule cplex_to_cbc:
    input: "results/Model_files/{model_run}_sorted.txt"
    output: "results/Model_files/{model_run}_cbc.txt"
    shell:
        "python scripts/cplextocbc.py {input} {output}"

rule generate_results:
    input: 
        results="results/Model_files/{model_run}_cbc.txt",
        datafile="output_data/{model_run}_modex.txt"
    params:
        scenario="{model_run}", folder="results/Scenario_CSVs/{model_run}"
    output: 
        emissions="results/Scenario_CSVs/{model_run}/AnnualEmissions.csv"
    script:
        "scripts/generate_results.py"
        
rule lcoe_and_emissions:
    input:
        results="results/Scenario_CSVs/{model_run}/AnnualEmissions.csv",
        model="input_data/Scenarios/{model_run}.xlsx"
    params:
        scenario="{model_run}",
        results_folder="results/Scenario_CSVs/{model_run}"
    output:
        lcoes="results/Scenario_CSVs/{model_run}/LCOEs.csv",
        emission_factors="results/Scenario_CSVs/{model_run}/GridEmissionFactors.csv",
        capacity="results/Scenario_CSVs/{model_run}/InstalledCapacity.csv",
        demand="results/Scenario_CSVs/{model_run}/ElectricityDemand.csv",
        generation="results/Scenario_CSVs/{model_run}/ElectricityGeneration.csv"
    script:
        "scripts/lcoe_and_emissions.py"
        
rule get_graphs:
    input:
        capacity="results/Scenario_CSVs/{model_run}/InstalledCapacity.csv",
        demand="results/Scenario_CSVs/{model_run}/ElectricityDemand.csv",
        generation="results/Scenario_CSVs/{model_run}/ElectricityGeneration.csv"
    params:
        scenario="{model_run}",
        output_folder="results/Main_results/Per_scenario/{model_run}"
    output:
        "results/Main_results/Per_scenario/{model_run}/Generation/CSV/ZMB.csv"
    script:
        "scripts/graphs.py"  
        