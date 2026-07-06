import biosteam as bst
import biobuilders as biob
from build_system import create_process_tea

process_tea, labor, labor_cost = create_process_tea()
process_sys = process_tea.system

# Simulation results
process_sys.simulate()
process_sys.show(flow='kg/hr',composition=True)

## Process diagram
process_sys.diagram(
    kind = 0,
    number = True,
    file = "results/grape_pomace_valorisation_process",
    format = "svg",
)

## Mass balance
mass_balance = biob.ProcessMassBalance(process_sys.streams)
streams_df = mass_balance.build_stream_table(excel_path="./results/grape_pomace_valorisation_mass_balance.xlsx")

## Unit operation results
unit_operations = biob.DisplayUnitsResults(process_sys.units)
unit_operations.save_report(filename="./results/grape_pomace_valorisation_unit_operations.xlsx")

## Check balance closure
total_in = sum(s.F_mass for s in process_sys.ins)
total_out = sum(s.F_mass for s in process_sys.outs)

error = abs(total_in - total_out)
relative_error = error / total_in * 100

print(f"Total in  : {total_in:.3f} kg/h")
print(f"Total out : {total_out:.3f} kg/h")
print(f"Error     : {error:.3f} kg/h")
print(f"Rel error : {relative_error:.3f} %")

# Techno-economic results
tea = process_tea
sys = process_sys
hours = tea.operating_hours
USD_to_EUR = 0.86

## Feedstock
grape_pomace = sys.units[0].ins[0]

## Products
phenolics = sys.units[14].outs[0]
protein_powder = sys.units[5].outs[0]

## Unit operations
enzyme_unit = sys.units[1]
nades_feed_unit = sys.units[7]
nades_loss_unit = sys.units[9]
flash_nades_cooling = sys.units[8]
flash_antisolvent_cooling = sys.units[20]
flash_steam_unit = sys.units[17]
precipitation_unit = sys.units[11]

## Results object
tea_results = biob.ResultsTEA(tea.get_cashflow_table(), tea)

## Capacity-related CAPEX
capex_by_capacity = {
    "input_1": tea_results.CAPEX_related_to_plant_capacity(input_stream="1"),
    "output_15": tea_results.CAPEX_related_to_plant_capacity(output_stream="15"),
    "output_35": tea_results.CAPEX_related_to_plant_capacity(output_stream="35"),
}

## Production costs
production_costs = tea_results.production_costs(
    streams=[phenolics],
    units=1,
    depreciation=False,
)

# Breakdown of production costs
tea_results.plot_production_costs_scenarios()

## Annual pomace treated
annual_pomace_kg = grape_pomace.F_mass * hours
annual_phenolics_kg = phenolics.F_mass * hours
annual_pomace_tonne = annual_pomace_kg / 1000

## Main economic results
npv_MM_EUR = tea.NPV * USD_to_EUR / 1e6
capex_MM_EUR = (tea.FCI + tea.working_capital) * USD_to_EUR / 1e6
opex_EUR_per_kg_pomace = tea.AOC * USD_to_EUR / annual_pomace_kg
msp_EUR_per_kg = tea.solve_price(phenolics) * USD_to_EUR

## Utility cost in AOC % basis
def utility_cost_fraction(unit, hu_index=0):
    """Return heat utility contribution as % of AOC."""
    return unit.heat_utilities[hu_index].cost * hours / tea.AOC * 100

cooling_nades_pct = utility_cost_fraction(flash_nades_cooling)
cooling_antisolvent_pct = utility_cost_fraction(flash_antisolvent_cooling)
low_pressure_steam_pct = utility_cost_fraction(flash_steam_unit)

# Main reactives cost
enzyme = enzyme_unit.ins[1]
fresh_nades = nades_feed_unit.ins[0]

enzyme_cost_fraction = enzyme.F_mass * hours * enzyme.price / tea.AOC * 100
fresh_nades_cost_fraction = fresh_nades.F_mass * fresh_nades.price * hours / tea.AOC * 100

## Protein sales
protein_sales_EUR = protein_powder.F_mass * hours * protein_powder.price * USD_to_EUR 
protein_sales_pct_AOC = protein_powder.F_mass  * protein_powder.price * hours / tea.AOC * 100

## Enzyme cost covered by protein
enzyme_cost_EUR = enzyme.F_mass * enzyme.price * hours * USD_to_EUR
protein_to_enzyme_revenue_ratio = protein_sales_EUR / enzyme_cost_EUR

## Nades losses
nades_components = ("Choline_Lactate", "Water")

nades_losses = (
    nades_loss_unit.outs[0].imass["Choline_Lactate"]
    + sys.units[9].outs[0].imass["Water"]
)

nades_input = (
    sys.units[9].ins[0].imass["Choline_Lactate"]
    + sys.units[9].ins[0].imass["Water"]
)

nades_loss_fraction = nades_losses / nades_input

## Print results
print("\n=== Techno-economic results ===")
print(f"Net present value: {npv_MM_EUR:.3f} MM EUR")
print(f"CAPEX: {capex_MM_EUR:.3f} MM EUR")
print(f"OPEX: {opex_EUR_per_kg_pomace:.3f} EUR/kg pomace")
print(f"MSP: {msp_EUR_per_kg:.3f} EUR/kg phenolics")
print(f"Plant capacity: {annual_pomace_tonne:.2f} tonne pomace/year")
print(f"Labor: {labor} workers")
print(f"Labor cost: {labor_cost*0.86} EUR/year")
tea_results.TEA_report(
    True,
    "./results/Grape_pomace_valorisation_cashflow.xlsx",
)

print("\n=== Current specifications ===")
print(f"Solvent ratio spec: {nades_feed_unit.specifications[0].args}")
print(f"Antisolvent ratio spec: {sys.units[10].specifications[0].args}")
print(f"Precipitation yield: {precipitation_unit.reaction.X}")

print("\n=== Cost contributions ===")
print(f"Cooling NADES after flash: {cooling_nades_pct:.2f}% AOC")
print(f"Cooling antisolvent after flash: {cooling_antisolvent_pct:.2f}% AOC")
print(f"Low-pressure steam for flash: {low_pressure_steam_pct:.2f}% AOC")
print(f"Fresh NADES: {fresh_nades_cost_fraction:.2f}% AOC")
print(f"Enzyme: {enzyme_cost_fraction:.2f}% AOC")

print("\n=== By-products and losses ===")
print(f"NADES loss fraction: {nades_loss_fraction:.4f}")
print(f"Protein powder sales: {protein_sales_pct_AOC:.2f}% AOC")
print(f"Protein sales: {protein_sales_EUR:.2f} EUR/year")
print(f"Enzyme cost: {enzyme_cost_EUR:.2f} EUR/year")
print(f"Protein sales / enzyme cost: {protein_to_enzyme_revenue_ratio:.3f}")

print("\n=== Prices ===")
print(f"Fresh NADES price: {fresh_nades.price:.4f} USD/kg")
print(f"Fresh enzyme price: {enzyme.price:.4f} USD/kg")
print(f"Peptide powder price: {protein_powder.price:.4f} USD/kg")