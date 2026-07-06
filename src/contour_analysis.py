import biosteam as bst
import biobuilders as biob
from build_model import build_model
from biobuilders.diagrams.contourplot import plot_contourf

model = build_model()
model.system.simulate()
model.system.show(composition=True, flow="kg/hr")
contour = biob.ContourAnalysis(model)

results = contour.evaluate_contour(
    parameter_x = "Sle402 solvent ratio",
    parameter_y = "R502 antisolvent ratio",
    nx = 30,
    ny = 30,
    xbounds=(5.,25.),
    ybounds=(10.,32.),
    reset_each_run = True,
    order = "row",
)

figs = plot_contourf(
    results = results,
    indicators = ["MSP"],
    title = "Effect of solvent and antisolvent ratio on MSP",
    xlabel = "Solvent ratio [kg NADES/kg grape pomace]",
    ylabel = "Antisolvent ratio [kg Antisolvent/kg extract]",
    levels = 10,
    path = "./results/",
    desired_ind = {"MSP": {"market price": 65.0}},
    use_results_baseline = True,
)