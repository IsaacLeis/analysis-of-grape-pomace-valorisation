import biobuilders as biob
from build_model import build_model, run_experiment

model = build_model()
model = run_experiment(model,5000)

uncertainty_plot = biob.UncertaintyPlotter(model.table.copy(), simplified_index={})
uncertainty_plot.plot_correlation_matrix(path = 'results/',filename = 'INRAE_kendall', indicators = ["OPEX [10^6 USD]","MSP [EUR/kg]","CAPEX [10^6 USD]"])
uncertainty_plot.plot_distribution_and_stats(["OPEX [10^6 USD]","MSP [EUR/kg]","CAPEX [10^6 USD]"], save_path = 'results/')

full_problem_space = model.table.copy()

sensitivity = biob.StandRegCoeffs(
    df = full_problem_space,
    parameter_cols = full_problem_space.columns[:-3].tolist(),
    indicator_cols = full_problem_space.columns[-3:].tolist(),
    simplified_labels = {}
)
sensitivity_df = sensitivity.calculate_src()
sensitivity.plot_src(sensitivity_df, path = './results/', filename_tag = 'INRAE')
print(sensitivity_df)

# Rows with at least SRC > +- 0.1
sensitivity_abs_gt_01 = sensitivity_df[sensitivity_df.abs().gt(0.1).any(axis=1)]

print("\nFilas con al menos un valor absoluto > 0.1:")
print(sensitivity_abs_gt_01)

sensitivity.plot_src(
    sensitivity_df,
    path='./results/',
    filename_tag='INRAE'
)

print("\nSRC completo:")
print(sensitivity_df)