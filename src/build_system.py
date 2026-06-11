import biosteam as bst
import biobuilders as biob
from process_chemicals import process_chemicals

# Load parameters
parameters = biob.parse_case_spec("./grape_value_chain/data/grape_pomace.json")
units_params = parameters["unit_operations"]
tea_params = parameters["techno_economic"]
prices = parameters["stream_prices"]

def create_system(ID="Grape pomace valorisation"):
    # Initialise thermo
    biob.ChemicalsManager(process_chemicals).set_thermodynamics()

    # Streams
    grape_pomace = bst.Stream(
        ID = '1',
        units = 'kg/hr',
        Grape_Pomace = units_params["GR101"]["Flow"]
    )
    enzyme = bst.Stream(
        ID = '3',
        units = 'kg/hr',
        Protamex = units_params["R201"]["E_S_Ratio"] * units_params["GR101"]["Flow"]
    )
    water = bst.Stream(
        ID = '4',
        units = 'kg/hr',
        H2O = units_params["R201"]["Water_Ratio"] * units_params["GR101"]["Flow"]
    )
    sodium_hydroxide = bst.Stream(
        ID = '5',
        units = 'kg/hr',
        NaOH = units_params["R201"]["NaOH_Ratio"] * units_params["GR101"]["Flow"]
    )
    fresh_nades = bst.Stream(
        ID = '21'
    )
    recycled_nades = bst.Stream(
        ID = '41',
        units = 'kg/hr',
    )
    antisolvent_recycled = bst.Stream(
        ID = '47',
        units = 'kg/hr',
    )
    precipitate_washing_recycled = bst.Stream(
        ID = '48',
        units = 'kg/hr',
        Ethyl_Acetate = ( 0.005 *
            units_params["R502"]["AcetEt"] *
            units_params["R502"]["Antisolvent_Ratio"] *
            units_params["R403"]["Solvent_Ratio"] *
            units_params["GR101"]["Flow"]
        ),
        Ethanol = ( 0.005 *
            units_params["R502"]["EtOH"] *
            units_params["R502"]["Antisolvent_Ratio"] *
            units_params["R403"]["Solvent_Ratio"] *
            units_params["GR101"]["Flow"]
        ),
    )

    # Protein extraction subsystem
    GR101 = biob.AttritionMill(
        ID = 'GR-101',
        ins = (grape_pomace,),
        outs = ("2",)
    )
    GR101.power_consumption = units_params["GR101"]["Specific_Energy"]
    R201 = biob.BatchEnzymaticTreatment(
        ID = 'R-201',
        ins = (GR101-0,enzyme,water,sodium_hydroxide),
        outs = ("6",),
        time = units_params["R201"]["Time"],
        time_loading = units_params["R201"]["Loading"],
        time_CIP = units_params["R201"]["CIP"],
        operating_T = units_params["R201"]["T"],
        reaction = bst.Reaction(
            reaction = '1 Grape_Pomace -> 1 Peptides',
            reactant = 'Grape_Pomace',
            basis = 'wt',
            X = 0.0145
        ),
        N_reactors = 2,
        kW_per_m3 = 0.33
    )
    CE202 = biob.SolidsCentrifuge(
        ID = 'CE-202',
        ins = (R201-0,),
        outs = ("8","7"),
        split = {
            "Grape_Pomace": units_params["CE202"]["Split_Grape_Pomace"],
        },
        solids = set(units_params["CE202"]["Solids"]),
        moisture_content = units_params["CE202"]["Moisture"],
        moisture_ID = ("Water",),
        solute_ID = (),
        kWh_per_kg = units_params["CE202"]["Specific_Energy"]
    )
    MF301 = biob.MembraneFiltration(
        ID = 'MF-301',
        ins = (CE202-1,),
        outs = ("9","10"),
        type = units_params["MF301"]["Type"],
        solids = units_params["MF301"]["Solids"],
        solids_retained = units_params["MF301"]["Solids_Retained"],
        solvent_IDs = ["Water"],
        solute_IDs = ["NaOH"]
    )
    MF301.kWh_per_kg = units_params["MF301"]["Specific_Energy"]
    E302 = biob.MultiEffectEvaporator(
        ID = 'E-302',
        ins = (MF301-0,),
        outs = ("11","12"),
        P = biob.solve_operating_pressures_multieffectevaporator(
            n_effects = units_params["E302"]["Effects"], 
            method_params = units_params["E302"]["Antoine"]
        )[0],
        V = units_params["E302"]["V"],
        V_definition = units_params["E302"]["V_Def"],
        chemical = units_params["E302"]["Chemical"]
    )
    D303 = biob.DrumDryer(
        ID = 'D303', 
        ins = (E302-0,"13","14"), 
        outs = ("15","16","17"), 
        moisture_content = units_params["D303"]["Moisture"],
        split = {
            "Peptides": units_params["D303"]["Split_Peptides"],
            "NaOH": units_params["D303"]["Split_NaOH"]
        },
        RH = units_params["D303"]["RH"]
    )
    
    # Phenolics extraction subsystem
    RVF401 = biob.RotaryVacuumFilter(
        ID = 'RVF-401',
        ins = (CE202-0, '18'),
        outs = ('19', '20'),
        split = units_params["RVF401"]["Split"],
        moisture_content = units_params["RVF401"]["Moisture"],
        moisture_ID = ("Water",),
        solute_ID = (),
        solids = tuple(units_params["RVF401"]["Solids"])
    )
    TK402 = biob.MixTank(
        ID = 'TK-401',
        ins = (fresh_nades,recycled_nades),
        outs = ('22'),
        tau = units_params["TK402"]["Tau"]
    )
    SLE403 = biob.ExtractionReactor(
        ID = 'SLE-402',
        ins = (RVF401-0,TK402-0),
        outs = ('23',),
        extract_reaction = bst.Reaction(
            'Grape_Pomace -> Phenolics_Extract',
            reactant = 'Grape_Pomace',
            X = units_params["R403"]["Yield"],
            basis = 'wt'
        ),
        tau = units_params["R403"]["Tau"],
        operating_T = units_params["R403"]["T"],
        operating_P = units_params["R403"]["P"],
        kW_per_m3 = units_params["R403"]["kW_m3"],
    )
    CE404 = biob.SolidsCentrifuge(
        ID = 'CE-403',
        ins = (SLE403-0,),
        outs = ('24','25'),
        split = units_params["CE404"]["Split"],
        solids = units_params["CE404"]["Solids"],
        moisture_content = units_params["CE404"]["Moisture"],
        moisture_ID = units_params["CE404"]["Moisture_ID"],
        solute_ID = units_params["CE404"]["Solute_ID"],
    )
    CE404.kWh_per_kg = units_params["CE404"]["Specific_Energy"]
    TK501 = biob.MixTank(
        ID = 'TK-501',
        ins = ('26', antisolvent_recycled),
        outs = ('27',),
        tau = units_params["TK501"]["Tau"]
    )
    R502 = biob.ContinuousStirredTankReactor(
        ID = 'R-502',
        ins = (CE404-1, TK501-0),
        outs = ('28',),
        reaction = bst.Reaction(
            reaction = 'Phenolics_Extract -> Phenolics_Precipitate',
            reactant = 'Phenolics_Extract',
            basis = 'wt',
            X = units_params["R502"]["PrecipitationPhenolics"],
        ),
        tau = units_params["R502"]["HRT"],
        operating_T = units_params["R502"]["T"],
        operating_P = units_params["R502"]["P"],
        kW_per_m3 = biob.agitator_volumetric_power_determination(
            method = "Sinnot heuristics",
            agitation_type = "Medium"
        )
    )
    RVF503 = biob.RotaryVacuumFilter(
        ID = 'RVF-503',
        ins = (R502-0,),
        outs = ('29','30'),
        split = units_params["RVF503"]["Split"],
        moisture_content = units_params["RVF503"]["Moisture"],
        moisture_ID = tuple(units_params["RVF503"]["Moisture_ID"]),
        solute_ID = tuple(units_params["RVF503"]["Solute_ID"]),
        solids = tuple(units_params["RVF503"]["Solids"]),
    )
    RVF503.kWh_per_kg = units_params["RVF503"]["Specific_Energy"]
    RVF504 = biob.RotaryVacuumFilter(
        ID = 'RVF-504',
        ins = (RVF503-0,precipitate_washing_recycled),
        outs = ('31','32'),
        split = units_params["RVF504"]["Split"],
        moisture_content = units_params["RVF504"]["Moisture"],
        moisture_ID = tuple(units_params["RVF504"]["Moisture_ID"]),
        solute_ID = tuple(units_params["RVF504"]["Solute_ID"]),
        solids = tuple(units_params["RVF504"]["Solids"]),
    )
    RVF504.kWh_per_kg = units_params["RVF504"]["Specific_Energy"]
    D505 = biob.DrumDryer(
        ID = 'D-505',
        ins = (RVF504-0,'33','34'),
        outs = ('35','36','37'),
        split = units_params["D505"]["Split"],
        RH = units_params["D505"]["RH"],
        moisture_content = units_params["D505"]["Moisture"],
        T = units_params["D505"]["T"],
        moisture_ID = "Ethyl_Acetate"
    )
    MX601 = bst.Mixer(
        ID = 'MX-601',
        ins = (RVF503-1,RVF504-1),
        outs = ('38',)
    )
    P602 = bst.Pump(
        ID = 'P-602',
        ins = (MX601-0,),
        outs = ('39'),
        P = units_params["P602"]["Pout"],
    )
    F603 = biob.Flash(
        ID = 'F-603',
        ins = (P602-0,),
        outs = ('40',recycled_nades),
        P = units_params["F603"]["P"],
        T = units_params["F603"]["T"]
    )
    C604 = biob.GasAdsorptionColumn(
        ID = 'C-604',
        ins = (F603-0,'42'),
        outs = ('43','44'),
        adsorbed_fraction = units_params['C604']['Adsorbed_Fraction'],
        adsorbate = units_params['C604']['Adsorbate'],
        adsorbent = units_params['C604']['Adsorbent'],
        t_ads = units_params['C604']['Ads_Time'],
    #    P_ads = units_params['C403']['P_Ads'],
        T_ads = units_params['C604']['T_Ads'],
        isotherm_args = tuple(units_params['C604']['Isotherm_Args']),
        isotherm_model = units_params['C604']['Isotherm_Model'],
        t_regen = units_params['C604']['Regen_Time'],
    #    P_regen = units_params['C403']['P_Regen'],
        T_regen = units_params['C604']['T_Regen'],
        regeneration_isotherm_model = units_params['C604']['Regen_Isotherm_Model'],
        regeneration_isotherm_args = units_params['C604']['Regen_Isotherm_Args'],
        void_fraction = 0.55,
    )
    C604.L_D_ratio = 1.5
    VLV605 = bst.IsenthalpicValve(
        ID = 'VLV-605',
        ins = (C604-0,),
        outs = ('45',),
        P = units_params["VLV605"]["P"],
    )
    VLV605._graphics = bst.Unit._graphics
    E606 = bst.HXutility(
        ID = 'E-606',
        ins = (VLV605-0,),
        outs = ('46',),
        rigorous = False,
        T = units_params["E606"]["T"],
        V = 0.0
    )
    S607 = bst.Splitter(
        ID = 'S-607',
        ins = (E606-0,),
        outs = (precipitate_washing_recycled,antisolvent_recycled),
        split = units_params["S-607"]["Antisolvent_Washing_Ratio"]/(1+0.02)
    )

    @RVF401.add_specification(run = True, args = (units_params["RVF401"]["Water_Ratio"],))
    def set_washing_water(wash_ratio):
        wet_solids = RVF401.ins[0]
        wash_water = RVF401.ins[1]
        wash_water.imass["Water"] = (wet_solids.F_mass - wet_solids.imass["Water"]) * wash_ratio

    def update_nades(solvent_ratio):

        dry_pomace = CE202.outs[0].imass["Grape_Pomace"]
        water_from_pomace = CE202.outs[0].imass["Water"]

        target_total = solvent_ratio * dry_pomace

        water_fraction = units_params["R403"]["Water"]
        nades_fraction = (
            units_params["R403"]["Choline"]
            + units_params["R403"]["Lactate"]
        )

        target_water = target_total * water_fraction
        target_choline_lactate = target_total * nades_fraction

        recycled_water = recycled_nades.imass["Water"]
        recycled_choline_lactate = recycled_nades.imass["Choline_Lactate"]

        fresh_nades.imass["Water"] = max(
            target_water - recycled_water - water_from_pomace,
            0.
        )

        fresh_nades.imass["Choline_Lactate"] = max(
            target_choline_lactate - recycled_choline_lactate,
            0.
        )  
    TK402.add_specification(
        f=update_nades,
        args=(units_params["R403"]["Solvent_Ratio"],),
        run=True
    )

    def update_antisolvent(antisolvent_ratio):

        solids = (
            CE404.outs[1].F_mass
        )

        target_solvent = (
            antisolvent_ratio
            * solids
        )

        recycled_antisolvent = (
            antisolvent_recycled.imass["Ethanol"]
            +
            antisolvent_recycled.imass["Ethyl_Acetate"]
        )

        fresh_required = max(
            target_solvent
            - recycled_antisolvent,
            0.
        )

        TK501.ins[0].imass["Ethanol"] = (
            fresh_required
            * 0.27
        )

        TK501.ins[0].imass["Ethyl_Acetate"] = (
            fresh_required
            * 0.73
        )   
    TK501.add_specification(
        f=update_antisolvent,
        args=(units_params["R502"]["Antisolvent_Ratio"],),
        run=True
    )

    @C604.add_specification(run=True)
    def calculate_regeneration_fluid():
        feed = C604.ins[0]
        regen_fluid = C604.ins[1]

        water_removed = feed.imass["Water"] * units_params["C604"]["Adsorbed_Fraction"]
        regen_gas = water_removed * (1.61) / 0.8

        regen_fluid.imass["N2"] = regen_gas * 0.79
        regen_fluid.imass["O2"] = regen_gas * 0.21

    # Global system
    grape_pomace_valorisation = bst.System(
        ID,
        path=(
            GR101, R201, CE202, MF301, E302, D303,
            RVF401, TK402, SLE403, CE404, TK501, R502, RVF503, RVF504, D505,
            MX601, P602, F603, C604, VLV605, E606, S607,
        ),
        recycle=[recycled_nades, antisolvent_recycled]
    )

    return grape_pomace_valorisation

def build_TEA(system: bst.System):
    process_settings = biob.ProcessSettingsManager(
        CEPCI = tea_params["CEPCI"],
        electricity_price = prices["Electricity"],
        streams_price = {
            system.units[1].ins[1]: biob.calculate_mean_median_price(prices["Protamex"]),
            system.units[1].ins[2]: prices["Water"]["value"],
            system.units[1].ins[3]: prices["NaOH"]["value"],
            system.units[5].ins[2]: prices["Natural_Gas"]["value"],
            system.units[5].outs[0]: biob.calculate_mean_median_price(prices["Peptides_Extract"]), 
            system.units[6].ins[1]: prices["Water"]["value"],
            system.units[7].ins[0]: biob.calculate_stream_price(
                {
                    0.219 * units_params["R403"]["Choline"]/(units_params["R403"]["Choline"]+units_params["R403"]["Lactate"]): biob.calculate_mean_median_price(prices["Choline_Chloride"]),
                    0.219 * units_params["R403"]["Lactate"]/(units_params["R403"]["Choline"]+units_params["R403"]["Lactate"]): biob.calculate_mean_median_price(prices["Lactic_Acid"]),
                    0.781 : prices["Water"]["value"],
                }
            ),
            system.units[10].ins[0]: biob.calculate_stream_price(
                {
                    units_params["R502"]["EtOH"]: prices["Ethanol"]["value"],
                    units_params["R502"]["AcetEt"]: prices["Ethyl_Acetate"]["value"]
                }
            ),
            system.units[14].ins[2]: prices["Natural_Gas"]["value"],
            system.units[14].outs[0]: prices["Phenolics_Extract"]
        }
    )
    process_settings.load_settings()

    labor, labor_equip = biob.calculate_labor_requirements(
        {
            "Centrifugal separators": 2,
            "Rotatory & belt filter": 2,
            "Process vessel": 2,
            "Continuous reactor": 2,
            "Heat exchanger": 1,
            "Evaporator": 2,
            "Batch reactor": 2,
            "Rotatory dryer": 2,
        }
    )
    labor_cost = labor * tea_params["Wage"]

    techno_economic = biob.TEA(
        system = system,
        IRR =               tea_params["IRR"],
        duration = (        tea_params["Duration"]["Start"],
                            tea_params["Duration"]["End"]
        ), 
        labor_cost =        labor_cost,
        depreciation =      tea_params["Depreciation"],
        supplies =          tea_params["Supplies"],             
        startup_months =    tea_params["Startup_Months"],                       
        startup_salesfrac = tea_params["Startup_Sales"],     
        finance_fraction =  tea_params["Finance_Frac"], 
        finance_interest =  tea_params["Finance_Interest"],
        finance_years =     tea_params["Finance_Years"],
        inflation_rate =    tea_params["Inflation"],      
    )

    return techno_economic