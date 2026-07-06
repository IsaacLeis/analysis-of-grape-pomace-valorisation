import biosteam as bst
import biobuilders as biob
import numpy as np
from build_system import create_process_tea, units_params, prices

techno_economic = create_process_tea()

def build_model():
    """Build and return BioSTEAM model for UNIBO system"""
    model = bst.Model(system=techno_economic.system)
    product_stream = techno_economic.system.outs[1]

    # Indicators
    @model.indicator(units = 'EUR/kg')
    def MSP():
        return techno_economic.solve_price(product_stream)*0.86
    
    @model.indicator(units = '10^6 USD')
    def OPEX():
        return techno_economic.production_costs(product_stream).sum()
    
    @model.indicator(units = '10^6 USD')
    def CAPEX():
        return techno_economic.TCI

    # Parameters
    ## Prices
    enzyme = techno_economic.system.units[1].ins[1]
    enzyme_price = biob.calculate_mean_median_price(prices["Protamex"])
    @model.parameter(
        element = enzyme,
        units = "USD/kg",
        bounds = (6.75,21.46),
        baseline = enzyme_price,
        distribution = 'uniform',
        coupled = True,
    )
    def set_enzyme_price(enzyme_price):
        enzyme.price = enzyme_price

    water_r201 = techno_economic.system.units[1].ins[2]
    water_rvf401 = techno_economic.system.units[6].ins[0]
    water_price = prices["Water"]["value"]
    @model.parameter(
        element = None,
        units = "USD/kg",
        bounds = (prices["Water"]["lb"],prices["Water"]["ub"]),
        baseline = water_price,
        distribution = 'uniform',
        coupled = True,
    )
    def set_water_price(water_price):
        water_r201.price = water_price
        water_rvf401.price = water_price
    
    naoh = techno_economic.system.units[1].ins[3]
    naoh_price = prices["NaOH"]["value"]
    @model.parameter(
        element = naoh,
        units = "USD/kg",
        bounds = (prices["NaOH"]["lb"],prices["NaOH"]["ub"]),
        baseline = naoh_price,
        distribution = 'uniform',
        coupled = True,
    )
    def set_NaOH_price(NaOH_price):
        naoh.price = NaOH_price

    natural_gas_d303 = techno_economic.system.units[5].ins[2]
    natural_gas_d505 = techno_economic.system.units[14].ins[2]
    natural_gas_price = 1.58
    @model.parameter(
        element = None,
        units = "USD/kg",
        bounds = (1.75,9.39),
        baseline = natural_gas_price,
        distribution = 'uniform',
        coupled = True,
    )
    def set_natural_gas_price(natural_gas_price):
        natural_gas_d303.price = natural_gas_price
        natural_gas_d505.price = natural_gas_price
    
    nades = techno_economic.system.units[7].ins[0]
    nades_price = 2.4149
    @model.parameter(
        element = None,
        units = "USD/kg",
        bounds = (1.5591,3.4310),
        baseline = nades_price,
        distribution = 'uniform',
        coupled = True,
    )
    def set_nades_price(nades_price):
        nades.price = nades_price
    
    antisolvent = techno_economic.system.units[10].ins[0]
    antisolvent_price = 1.2445
    @model.parameter(
        element = None,
        units = "USD/kg",
        bounds = (0.9199,1.6538),
        baseline = antisolvent_price,
        distribution = 'uniform',
        coupled = True,
    )
    def set_antisolvent_price(antisolvent_price):
        antisolvent.price = antisolvent_price
    
    ## Yield
    r201 = techno_economic.system.units[1]
    r201_yield = r201.reaction.X
    @model.parameter(
         element = r201,
         units = "kg peptides/kg pomace",
         bounds = (0.01337,0.01553),
         baseline = r201_yield,
         distribution = 'uniform',
         coupled = True
     )
    def set_r201_extraction_yield(r201_extraction_yield):
        r201.reaction.X = r201_extraction_yield
    
    sle402 = techno_economic.system.units[8]
    sle402_yield = units_params["R403"]["Yield"]
    @model.parameter(
        element = sle402,
        units = "kg phenolics extract/kg pomace",
        bounds = (0.115,0.121),
        baseline = sle402_yield,
        distribution = 'uniform',
        coupled = True
    )
    def set_sle402_extraction_yield(sle402_extraction_yield):
        sle402.extract_react.X = sle402_extraction_yield

    r502 = techno_economic.system.units[11]
    r502_yield = r502.reaction.X
    @model.parameter(
        element = r502,
        units = "kg phenolics precipitated/kg phenolics extract",
        bounds = (0.8415, 1.),
        baseline = r502_yield,
        distribution = 'uniform',
        coupled = True,
    )
    def set_r502_precipitation_yield(r502_precipitation_yield):
        r502.reaction.X = r502_precipitation_yield

    ## Dosage
    gr101 = techno_economic.system.units[0]
    r201_dosage = units_params["R201"]["E_S_Ratio"]
    @model.parameter(
        element = r201,
        units = "kg enzyme/kg pomace",
        bounds = (0.04,0.06),
        baseline = r201_dosage,
        distribution = 'uniform',
        coupled = True
    )
    def set_r201_enzyme_dosage(r201_enzyme_dosage):
        r201.ins[1].F_mass = r201_enzyme_dosage * gr101.ins[0].F_mass

    ## time
    r201_time = r201.time
    @model.parameter(
        element = r201,
        units = "h",
        bounds = (1.70,2.30),
        baseline = r201_time,
        distribution = 'uniform',
        coupled = True
    )
    def set_r201_batch_time(r201_batch_time):
        r201.time = r201_batch_time
    
    ## HRT
    sle201_hrt = sle402.tau
    @model.parameter(
        element = sle402,
        units = "h",
        bounds = (1.50,1.92),
        baseline = sle201_hrt,
        distribution = 'uniform',
        coupled = True
    )
    def set_sle402_hrt(sle402_hrt):
        sle402.tau = sle402_hrt

    r502_hrt = r502.tau
    @model.parameter(
        element = r502,
        units = "h",
        bounds = (0.042,0.125),
        baseline = r502_hrt,
        distribution = 'uniform',
        coupled = True
    )
    def set_r502_hrt(r502_hrt):
        r502.tau = r502_hrt

    ## Temperature
    r201_temperature = r201.operating_T
    @model.parameter(
        element = r201,
        units = "h",
        bounds = (r201_temperature*0.85,r201_temperature*1.15),
        baseline = r201_temperature,
        distribution = 'uniform',
        coupled = True
    )
    def set_r201_temperature(r201_temperature):
        r201.operating_T = r201_temperature
    
    sle402_temperature = sle402.operating_T
    @model.parameter(
        element = sle402,
        units = "K",
        bounds = (273.15+40, 273.15+60),
        baseline = sle402_temperature,
        distribution = 'uniform',
        coupled = True
    )
    def set_sle402_temperature(sle402_temperature):
        sle402.operating_T = sle402_temperature

    r502_temperature = r502.operating_T
    @model.parameter(
        element = r502,
        units = "K",
        bounds = (273.15,273.15+10),
        baseline = r502_temperature,
        distribution = "uniform",
        coupled = True
    )
    def set_r502_temperature(r502_temperature):
        r502.operating_T = r502_temperature
    
    ## Moistures
    ce403 = techno_economic.system.units[9]
    ce403_moisture = ce403.moisture_content
    @model.parameter(
        element = ce403,
        units = "kg liquid/kg dry weight",
        bounds = (0.28,0.35),
        baseline = ce403_moisture,
        distribution = 'uniform',
        coupled = True
    )
    def set_ce403_moisture(ce403_moisture):
        ce403.moisture_content = ce403_moisture

    rvf503 = techno_economic.system.units[12]
    rvf503_moisture = rvf503.moisture_content
    @model.parameter(
        element = rvf503,
        units = "kg liquid/kg dry weight",
        bounds = (0.25,0.32),
        baseline = rvf503_moisture,
        distribution = 'uniform',
        coupled = True
    )
    def set_rvf503_moisture(rvf503_moisture):
        rvf503.moisture_content = rvf503_moisture

    rvf504 = techno_economic.system.units[13]
    rvf504_moisture = rvf504.moisture_content
    @model.parameter(
        element = rvf504,
        units = "kg liquid/kg dry weight",
        bounds = (0.25,0.32),
        baseline = rvf504_moisture,
        distribution = 'uniform',
        coupled = True
    )
    def set_rvf504_moisture(rvf504_moisture):
        rvf504.moisture_content = rvf504_moisture

    d505 = techno_economic.system.units[14]
    d505_moisture = d505.moisture_content
    @model.parameter(
        element = d505,
        units = "kg liquid/kg dry weight",
        bounds = (0.0125,0.0375),
        baseline = d505_moisture,
        distribution = 'uniform',
        coupled = True
    )
    def set_d505_moisture(d505_moisture):
        d505.moisture_content = d505_moisture

    ## Volume concentration factor
    mf301 = techno_economic.system.units[3]
    mf301_vcf = mf301.VCF
    @model.parameter(
        element = mf301,
        units = 'times',
        bounds = (mf301_vcf*0.50,mf301_vcf*1.50),
        baseline = mf301_vcf,
        distribution = 'uniform',
        coupled = True
    )
    def set_mf301_vcf(mf301_vcf):
        mf301.VCF = mf301_vcf

    ## Relative humidity
    d303 = techno_economic.system.units[5]
    d303_relative_humidity = d303.RH
    @model.parameter(
        element = d303,
        units = "-",
        bounds = (0.20,0.80),
        baseline = d303_relative_humidity,
        distribution = 'uniform',
        coupled = True
    )
    def set_d303_relative_humidity(d303_relative_humidity):
        d303.RH = d303_relative_humidity

    d505_relative_humidity = d505.RH
    @model.parameter(
        element = d505,
        units = "-",
        bounds = (0.20,0.80),
        baseline = d505_relative_humidity,
        distribution = 'uniform',
        coupled = True
    )
    def set_d505_relative_humidity(d505_relative_humidity):
        d505.RH = d505_relative_humidity

    ## Energy
    gr101 = techno_economic.system.units[0]
    gr101_specific_energy = gr101.power_consumption
    @model.parameter(
        element = gr101,
        units = 'kWh/kg',
        bounds = (0.008,0.016),
        baseline = gr101_specific_energy,
        distribution = 'uniform',
        coupled = True
    )
    def set_gr101_specific_energy(gr101_specific_energy):
        gr101.power_consumption = gr101_specific_energy

    r201_volumetric_power = r201.kW_per_m3
    @model.parameter(
        element = r201,
        units = 'kW/m3',
        bounds = (r201_volumetric_power*0.50,r201_volumetric_power*1.50),
        baseline = r201_volumetric_power,
        distribution = 'uniform',
        coupled = True
    )
    def set_r201_volumetric_power(r201_volumetric_power):
        r201.kW_per_m3 = r201_volumetric_power

    ce202 = techno_economic.system.units[2]
    ce202_specific_energy = ce202.kWh_per_kg
    @model.parameter(
        element = ce202,
        units = 'kWh/kg',
        bounds = (0.001,0.010),
        baseline = ce202_specific_energy,
        distribution = 'uniform',
        coupled = True
    )
    def set_ce202_specific_energy(ce202_specific_energy):
        ce202.power_consumption = ce202_specific_energy

    mf301_tmp = mf301.TMP
    @model.parameter(
        element = mf301,
        units = 'Pa',
        bounds = (mf301_tmp*0.50,mf301_tmp*1.50),
        baseline = mf301_tmp,
        distribution = 'uniform',
        coupled = True
    )
    def set_mf301_tmp(mf301_tmp):
        mf301.TMP = mf301_tmp
    
    mf301_dp = mf301.pressure_drop
    @model.parameter(
        element = mf301,
        units = 'Pa',
        bounds = (mf301_dp*0.50,mf301_dp*1.50),
        baseline = mf301_dp,
        distribution = 'uniform',
        coupled = True
    )
    def set_mf301_dp(mf301_dp):
        mf301.pressure_drop = mf301_dp

    sle402_volumetric_power = sle402.kW_per_m3
    @model.parameter(
        element = sle402,
        units = 'kW/m3',
        bounds = (0.015,0.025),
        baseline = sle402_volumetric_power,
        distribution = 'uniform',
        coupled = True
    )
    def set_sle201_volumetric_power(sle402_volumetric_power):
        sle402.kW_per_m3 = sle402_volumetric_power
    
    ce403_specific_energy = ce403.kWh_per_kg
    @model.parameter(
        element = ce403,
        units = 'kWh/kg',
        bounds = (0.001,0.010),
        baseline = ce403_specific_energy,
        distribution = 'uniform',
        coupled = True
    )
    def set_ce403_specific_energy(ce403_specific_energy):
        ce403.power_consumption = ce403_specific_energy
    
    r502_volumetric_power = r502.kW_per_m3
    @model.parameter(
        element = r502,
        units = 'kW/m3',
        bounds = (r502_volumetric_power*0.50,r502_volumetric_power*1.50),
        baseline = r502_volumetric_power,
        distribution = 'uniform',
        coupled = True
    )
    def set_r502_volumetric_power(r502_volumetric_power):
        r502.kW_per_m3 = r502_volumetric_power
    
    ## Solvent Ratio
    tk401 = techno_economic.system.units[7]
    sle402_solvent_ratio = 20.
    @model.parameter(
        element = sle402,
        units = 'kg solvent/kg solids',
        bounds = (17.5,22.5),
        baseline = sle402_solvent_ratio,
        distribution = 'uniform',
        coupled = True,
    )
    def set_sle402_solvent_ratio(sle402_solvent_ratio):
        techno_economic.system.empty_recycles()
        techno_economic.system.reset_cache()
        tk401.specifications[0].args = (sle402_solvent_ratio,)

    ## Antisolvent ratio
    tk501 = techno_economic.system.units[10]
    r502_antisolvent_ratio = units_params["R502"]["Antisolvent_Ratio"]
    @model.parameter(
        element = r502,
        units = 'kg antisolvent/kg NaDES',
        bounds = (22.725,30.3),
        baseline = r502_antisolvent_ratio,
        distribution = 'uniform',
        coupled = True
    )
    def set_r502_antisolvent_ratio(r502_antisolvent_ratio):
        techno_economic.system.empty_recycles()
        techno_economic.system.reset_cache()
        tk501.specifications[0].args = (r502_antisolvent_ratio,)

    ## Water ratio
    r201_water_ratio = units_params["R201"]["Water_Ratio"]
    @model.parameter(
        element = r201,
        units = 'kg Water/kg pomace',
        bounds = (r201_water_ratio*0.85,r201_water_ratio*1.15),
        baseline = r201_water_ratio,
        distribution = 'uniform',
        coupled = True
    )
    def set_r201_water_ratio(r201_water_ratio):
        r201.ins[2].imass["Water"] = r201_water_ratio * units_params["GR101"]["Flow"]
    
    ## NaOH ratio
    r201_naoh_ratio = units_params["R201"]["NaOH_Ratio"]
    @model.parameter(
        element = r201,
        units = 'kg NaOH/kg pomace',
        bounds = (r201_naoh_ratio*0.85,r201_naoh_ratio*1.15),
        baseline = r201_naoh_ratio,
        distribution = 'uniform',
        coupled = True
    )
    def set_r201_naoh_ratio(r201_naoh_ratio):
        r201.ins[2].imass["NaOH"] = r201_naoh_ratio * units_params["GR101"]["Flow"]
    
    ## Langmuir-Langmuir parameters
    c605 = techno_economic.system.units[18]
    c605_ka = c605.isotherm_args[1]
    @model.parameter(
        element = c605,
        units = '1/Pa',
        bounds = (c605_ka*0.50,c605_ka*1.50),
        baseline = c605_ka,
        distribution = 'uniform',
        coupled = True
    )
    def set_c605_ka(c605_ka):
        c605.isotherm_args[1] = c605_ka

    c605_kb = c605.isotherm_args[3]
    @model.parameter(
        element = c605,
        units = '1/Pa',
        bounds = (c605_kb*0.50,c605_kb*1.50),
        baseline = c605_kb,
        distribution = 'uniform',
        coupled = True
    )
    def set_c605_kb(c605_kb):
        c605.isotherm_args[3] = c605_kb
    
    c605_ka_regeneration = c605.regeneration_isotherm_args[1]
    @model.parameter(
        element = c605,
        units = '1/Pa',
        bounds = (c605_ka_regeneration*0.50,c605_ka_regeneration*1.50),
        baseline = c605_ka_regeneration,
        distribution = 'uniform',
        coupled = True
    )
    def set_c605_ka_regeneration(c605_ka_regeneration):
        c605.regeneration_isotherm_args[1] = c605_ka_regeneration

    c605_kb_regeneration = c605.regeneration_isotherm_args[3]
    @model.parameter(
        element = c605,
        units = '1/Pa',
        bounds = (c605_kb_regeneration*0.50,c605_kb_regeneration*1.50),
        baseline = c605_kb_regeneration,
        distribution = 'uniform',
        coupled = True
    )
    def set_c605_kb_regeneration(c605_kb_regeneration):
        c605.regeneration_isotherm_args[3] = c605_kb_regeneration

    c605_qmax = c605.isotherm_args[-1]
    @model.parameter(
        element = c605,
        units = '1/Pa',
        bounds = (c605_qmax*0.50,c605_qmax*1.50),
        baseline = c605_qmax,
        distribution = 'uniform',
        coupled = True
    )
    def set_c605_qmax(c605_qmax):
        c605.isotherm_args[-1] = c605_qmax
    
    c605_qmax_regeneration = c605.regeneration_isotherm_args[-1]
    @model.parameter(
        element = c605,
        units = '1/Pa',
        bounds = (c605_qmax_regeneration*0.50,c605_qmax_regeneration*1.50),
        baseline = c605_qmax_regeneration,
        distribution = 'uniform',
        coupled = True
    )
    def set_c605_qmax_regeneration(c605_qmax_regeneration):
        c605.regeneration_isotherm_args[-1] = c605_qmax_regeneration
    
    c605_time = c605.t_ads
    @model.parameter(
        element = c605,
        units = '1/Pa',
        bounds = (c605_time*0.50,c605_time*1.50),
        baseline = c605_time,
        distribution = 'uniform',
        coupled = True
    )
    def set_c605_time(c605_time):
        c605.t_ads = c605_time
    
    c605_time_regeneration = c605.t_regen
    @model.parameter(
        element = c605,
        units = '1/Pa',
        bounds = (c605_time_regeneration*0.50,c605_time_regeneration*1.50),
        baseline = c605_time_regeneration,
        distribution = 'uniform',
        coupled = True
    )
    def set_c605_time_regeneration(c605_time_regeneration):
        c605.t_regen = c605_time_regeneration

    ## Equipment cost
    gr101_base_cost = gr101.base_cost
    @model.parameter(
        element = gr101,
        units = '$',
        bounds = (gr101_base_cost*0.50,gr101_base_cost*1.50),
        baseline = gr101_base_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_gr101_base_cost(gr101_base_cost):
        gr101.base_cost = gr101_base_cost
    
    gr101_base_n_cost = gr101.base_n_cost
    @model.parameter(
        element = gr101,
        units = '-',
        bounds = (gr101_base_n_cost*0.50,gr101_base_n_cost*1.50),
        baseline = gr101_base_n_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_gr101_base_n_cost(gr101_base_n_cost):
        gr101.base_n_cost = gr101_base_n_cost

    r201_base_cost = r201.base_cost
    @model.parameter(
        element = r201,
        units = '$',
        bounds = (r201_base_cost*0.50,r201_base_cost*1.50),
        baseline = r201_base_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_r201_base_cost(r201_base_cost):
        r201.base_cost = r201_base_cost
    
    r201_base_n_cost = r201.base_n_cost
    @model.parameter(
        element = r201,
        units = '-',
        bounds = (r201_base_n_cost*0.50,r201_base_n_cost*1.50),
        baseline = r201_base_n_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_r201_base_n_cost(r201_base_n_cost):
        r201.base_n_cost = r201_base_n_cost
    
    ce202_base_cost = ce202.base_cost
    @model.parameter(
        element = ce202,
        units = '$',
        bounds = (ce202_base_cost*0.50,ce202_base_cost*1.50),
        baseline = ce202_base_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_ce202_base_cost(ce202_base_cost):
        ce202.base_cost = ce202_base_cost
    
    ce202_base_n_cost = ce202.base_n_cost
    @model.parameter(
        element = ce202,
        units = '-',
        bounds = (ce202_base_n_cost*0.50,ce202_base_n_cost*1.50),
        baseline = ce202_base_n_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_ce202_base_n_cost(ce202_base_n_cost):
        ce202.base_n_cost = ce202_base_n_cost
    
    mf301_base_cost = mf301.base_cost
    @model.parameter(
        element = mf301,
        units = '$',
        bounds = (mf301_base_cost*0.50,mf301_base_cost*1.50),
        baseline = mf301_base_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_mf301_base_cost(mf301_base_cost):
        mf301.base_cost = mf301_base_cost
    
    mf301_base_n_cost = mf301.base_n_cost
    @model.parameter(
        element = mf301,
        units = '-',
        bounds = (mf301_base_n_cost*0.50,mf301_base_n_cost*1.50),
        baseline = mf301_base_n_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_mf301_base_n_cost(mf301_base_n_cost):
        mf301.base_n_cost = mf301_base_n_cost

    e302 = techno_economic.system.units[4]
    e302_base_cost = e302.base_cost
    @model.parameter(
        element = e302,
        units = '$',
        bounds = (e302_base_cost*0.50,e302_base_cost*1.50),
        baseline = e302_base_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_e302_base_cost(e302_base_cost):
        e302.base_cost = e302_base_cost
    
    e302_base_n_cost = e302.base_n_cost
    @model.parameter(
        element = e302,
        units = '-',
        bounds = (e302_base_n_cost*0.50,e302_base_n_cost*1.50),
        baseline = e302_base_n_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_e302_base_n_cost(e302_base_n_cost):
        e302.base_n_cost = e302_base_n_cost

    d303_base_cost = d303.base_cost
    @model.parameter(
        element = d303,
        units = '$',
        bounds = (d303_base_cost*0.50,d303_base_cost*1.50),
        baseline = d303_base_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_d303_base_cost(d303_base_cost):
        d303.base_cost = d303_base_cost
    
    d303_base_n_cost = d303.base_n_cost
    @model.parameter(
        element = d303,
        units = '-',
        bounds = (d303_base_n_cost*0.50,d303_base_n_cost*1.50),
        baseline = d303_base_n_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_d303_base_n_cost(d303_base_n_cost):
        d303.base_n_cost = d303_base_n_cost

    sle402_base_cost = sle402.base_cost
    @model.parameter(
        element = sle402,
        units = '$',
        bounds = (sle402_base_cost*0.50,sle402_base_cost*1.50),
        baseline = sle402_base_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_sle402_base_cost(sle402_base_cost):
        sle402.base_cost = sle402_base_cost
    
    sle402_base_n_cost = sle402.base_n_cost
    @model.parameter(
        element = sle402,
        units = '-',
        bounds = (sle402_base_n_cost*0.50,sle402_base_n_cost*1.50),
        baseline = sle402_base_n_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_sle402_base_n_cost(sle402_base_n_cost):
        sle402.base_n_cost = sle402_base_n_cost

    ce403_base_cost = ce403.base_cost
    @model.parameter(
        element = ce403,
        units = '$',
        bounds = (ce403_base_cost*0.50,ce403_base_cost*1.50),
        baseline = ce403_base_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_ce403_base_cost(ce403_base_cost):
        ce403.base_cost = ce403_base_cost
    
    ce403_base_n_cost = ce403.base_n_cost
    @model.parameter(
        element = ce403,
        units = '-',
        bounds = (ce403_base_n_cost*0.50,ce403_base_n_cost*1.50),
        baseline = ce403_base_n_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_ce403_base_n_cost(ce403_base_n_cost):
        ce403.base_n_cost = ce403_base_n_cost

    r502_base_cost = r502.base_cost
    @model.parameter(
        element = r502,
        units = '$',
        bounds = (r502_base_cost*0.50,r502_base_cost*1.50),
        baseline = r502_base_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_r502_base_cost(r502_base_cost):
        r502.base_cost = r502_base_cost
    
    r502_base_n_cost = r502.base_n_cost
    @model.parameter(
        element = r502,
        units = '-',
        bounds = (r502_base_n_cost*0.50,r502_base_n_cost*1.50),
        baseline = r502_base_n_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_r502_base_n_cost(r502_base_n_cost):
        r502.base_n_cost = r502_base_n_cost

    rvf503_base_cost = rvf503.base_cost
    @model.parameter(
        element = rvf503,
        units = '$',
        bounds = (rvf503_base_cost*0.50,rvf503_base_cost*1.50),
        baseline = rvf503_base_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_rvf503_base_cost(rvf503_base_cost):
        rvf503.base_cost = rvf503_base_cost
    
    rvf503_base_n_cost = rvf503.base_n_cost
    @model.parameter(
        element = rvf503,
        units = '-',
        bounds = (rvf503_base_n_cost*0.50,rvf503_base_n_cost*1.50),
        baseline = rvf503_base_n_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_rvf503_base_n_cost(rvf503_base_n_cost):
        rvf503.base_n_cost = rvf503_base_n_cost

    rvf504_base_cost = rvf504.base_cost
    @model.parameter(
        element = rvf504,
        units = '$',
        bounds = (rvf504_base_cost*0.50,rvf504_base_cost*1.50),
        baseline = rvf504_base_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_rvf504_base_cost(rvf504_base_cost):
        rvf504.base_cost = rvf504_base_cost
    
    rvf504_base_n_cost = rvf504.base_n_cost
    @model.parameter(
        element = rvf504,
        units = '-',
        bounds = (rvf504_base_n_cost*0.50,rvf504_base_n_cost*1.50),
        baseline = rvf504_base_n_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_rvf504_base_n_cost(rvf504_base_n_cost):
        rvf504.base_n_cost = rvf504_base_n_cost

    d505_base_cost = d505.base_cost
    @model.parameter(
        element = d505,
        units = '$',
        bounds = (d505_base_cost*0.50,d505_base_cost*1.50),
        baseline = d505_base_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_d505_base_cost(d505_base_cost):
        d505.base_cost = d505_base_cost
    
    d505_base_n_cost = d505.base_n_cost
    @model.parameter(
        element = d505,
        units = '-',
        bounds = (d505_base_n_cost*0.50,d505_base_n_cost*1.50),
        baseline = d505_base_n_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_d505_base_n_cost(d505_base_n_cost):
        d505.base_n_cost = d505_base_n_cost

    f603 = techno_economic.system.units[17]
    f603_base_cost = f603.base_cost
    @model.parameter(
        element = f603,
        units = '$',
        bounds = (f603_base_cost*0.50,f603_base_cost*1.50),
        baseline = f603_base_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_f603_base_cost(f603_base_cost):
        f603.base_cost = f603_base_cost
    
    f603_base_n_cost = f603.base_n_cost
    @model.parameter(
        element = f603,
        units = '-',
        bounds = (f603_base_n_cost*0.50,f603_base_n_cost*1.50),
        baseline = f603_base_n_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_f603_base_n_cost(f603_base_n_cost):
        f603.base_n_cost = f603_base_n_cost
    
    c605 = techno_economic.system.units[18]
    c605_base_cost = c605.base_cost_vessel
    @model.parameter(
        element = c605,
        units = '$',
        bounds = (c605_base_cost*0.50,c605_base_cost*1.50),
        baseline = c605_base_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_c605_base_cost(c605_base_cost):
        c605.base_cost = c605_base_cost
    
    c605_base_n_cost = c605.base_n_cost_vessel
    @model.parameter(
        element = c605,
        units = '-',
        bounds = (c605_base_n_cost*0.50,c605_base_n_cost*1.50),
        baseline = c605_base_n_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_c605_base_n_cost(c605_base_n_cost):
        c605.base_n_cost = c605_base_n_cost

    tk401_base_cost = tk401.base_cost
    @model.parameter(
        element = tk401,
        units = '$',
        bounds = (tk401_base_cost*0.50,tk401_base_cost*1.50),
        baseline = tk401_base_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_tk401_base_cost(tk401_base_cost):
        tk401.base_cost = tk401_base_cost
    
    tk401_base_n_cost = tk401.base_n_cost
    @model.parameter(
        element = tk401,
        units = '-',
        bounds = (tk401_base_n_cost*0.50,tk401_base_n_cost*1.50),
        baseline = tk401_base_n_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_tk401_base_n_cost(tk401_base_n_cost):
        tk401.base_n_cost = tk401_base_n_cost

    tk501_base_cost = tk501.base_cost
    @model.parameter(
        element = tk501,
        units = '$',
        bounds = (tk501_base_cost*0.50,tk501_base_cost*1.50),
        baseline = tk501_base_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_tk501_base_cost(tk501_base_cost):
        tk501.base_cost = tk501_base_cost
    
    tk501_base_n_cost = tk501.base_n_cost
    @model.parameter(
        element = tk501,
        units = '-',
        bounds = (tk501_base_n_cost*0.50,tk501_base_n_cost*1.50),
        baseline = tk501_base_n_cost,
        distribution = 'uniform',
        coupled = True
    )
    def set_tk501_base_n_cost(tk501_base_n_cost):
        tk501.base_n_cost = tk501_base_n_cost

    return model

def run_experiment(model, n_samples=5000, seed=1234, rule='L'):
    """Build model evaluation."""

    # Sampling and evaluation
    np.random.seed(seed)
    samples = model.sample(n_samples,rule)
    model.load_samples(samples,sort=True)
    model.exception_hook = "raise"
    model.evaluate(notify=100)

    return model