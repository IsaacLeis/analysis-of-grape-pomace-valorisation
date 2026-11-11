import biobuilders as biob
from build_model import build_model, run_experiment
import warnings; warnings.filterwarnings('ignore')

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

sensitivity.plot_src(
    sensitivity_df,
    path='./results/',
    filename_tag='INRAE'
)

# Filter SRC values only for MSP
SRC_THRESHOLD = 0.080
MSP_COL = "MSP [EUR/kg]"

# |SRC| > 0.08 for MSP
mask = sensitivity_df[MSP_COL].abs().gt(SRC_THRESHOLD)
sensitivity_msp_abs_gt_008 = sensitivity_df.loc[mask, [MSP_COL]]

sensitivity_msp_abs_gt_008 = (
    sensitivity_msp_abs_gt_008
    .assign(abs_SRC=sensitivity_msp_abs_gt_008[MSP_COL].abs())
    .sort_values(by="abs_SRC", ascending=False)
    .drop(columns="abs_SRC")
)

print("|SRC| > 0.08 for MSP:")
print(sensitivity_msp_abs_gt_008)

sensitivity_msp_abs_gt_008.to_csv(
    "./results/INRAE_SRC_MSP_abs_gt_0p08.csv"
)

sensitivity_msp_abs_gt_008_for_plot = sensitivity_df.loc[
    sensitivity_msp_abs_gt_008.index
]

sensitivity.plot_src(
    sensitivity_msp_abs_gt_008_for_plot,
    path="./results/",
    filename_tag="INRAE_MSP_absSRC_gt_0p08",
)