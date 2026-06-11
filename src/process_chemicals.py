import biosteam as bst

nades = bst.Chemical(
    ID = 'Choline_Lactate',
    Cp = 2.2,
    rho = 1175,
    phase = 'l',
    phase_ref = 'l',
    MW = 319.78,
    default = True,
    search_db = False,
    Psat = 1e-9
)
deproteinised_pomace = bst.Chemical(
    ID = 'Deproteinised_Grape_Pomace',
    Cp = 1.364,
    rho = 1100,
    phase = 's',
    phase_ref = 's',
    MW = 1.,
    default = True,
    search_db = False
)
phenolics_extract = bst.Chemical(
    ID = 'Phenolics_Extract',
    Cp = 2.1,
    rho = 1500,
    phase = 's',
    phase_ref = 's',
    MW = 272.25,
    default = True,
    search_db = False
)
phenolics_precipitate = bst.Chemical(
    ID = 'Phenolics_Precipitate',
    Cp = 2.1,
    rho = 1500,
    phase = 's',
    phase_ref = 's',
    MW = 272.25,
    default = True,
    search_db = False
)
pomace = bst.Chemical(
    ID = 'Grape_Pomace',
    Cp = 1.364,
    rho = 1100,
    phase = 's',
    phase_ref = 's',
    MW = 1.,
    default = True,
    search_db = False
)
protamex = bst.Chemical(
    ID = 'Protamex',
    Cp = 1.364,
    rho = 1350,
    phase = 's',
    phase_ref = 's',
    MW = 1.,
    default = True,
    search_db = False,
)