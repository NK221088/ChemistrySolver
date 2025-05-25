"""
Data module containing solubility rules, Ksp values, and reagent information
for qualitative analysis and solubility calculations.
"""

# Solubility data for common cations with different anions
# True = soluble, False = insoluble
SOLUBILITY_DATA = {
    # Format: 'cation': {'anion': solubility_boolean}
    'Na+': {
        'Cl-': True,
        'Br-': True,
        'I-': True,
        'OH-': True,
        'SO4^2-': True,
        'CO3^2-': True,
        'PO4^3-': True,
        'S^2-': True,
        'NO3-': True,
        'CH3COO-': True
    },
    'K+': {
        'Cl-': True,
        'Br-': True,
        'I-': True,
        'OH-': True,
        'SO4^2-': True,
        'CO3^2-': True,
        'PO4^3-': True,
        'S^2-': True,
        'NO3-': True,
        'CH3COO-': True
    },
    'Ag+': {
        'Cl-': False,
        'Br-': False,
        'I-': False,
        'OH-': False,
        'SO4^2-': True,  # Silver sulfate is slightly soluble
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Sr^2+': {
        'Cl-': True,
        'Br-': True,
        'I-': True,
        'OH-': True,  # Slightly soluble
        'SO4^2-': False,  # Strontium sulfate is insoluble
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Ba^2+': {
        'Cl-': True,
        'Br-': True,
        'I-': True,
        'OH-': True,
        'SO4^2-': False,  # Barium sulfate is insoluble
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Ca^2+': {
        'Cl-': True,
        'Br-': True,
        'I-': True,
        'OH-': False,
        'SO4^2-': True,  # Slightly soluble
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Cu^2+': {
        'Cl-': True,
        'Br-': True,
        'I-': False,
        'OH-': False,
        'SO4^2-': True,
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Fe^2+': {
        'Cl-': True,
        'Br-': True,
        'I-': True,
        'OH-': False,
        'SO4^2-': True,
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Fe^3+': {
        'Cl-': True,
        'Br-': True,
        'I-': True,
        'OH-': False,
        'SO4^2-': True,
        'CO3^3-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Pb^2+': {
        'Cl-': False,  # Slightly soluble in cold water
        'Br-': False,
        'I-': False,
        'OH-': False,
        'SO4^2-': False,  # Lead sulfate is insoluble
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Mg^2+': {
        'Cl-': True,
        'Br-': True,
        'I-': True,
        'OH-': False,
        'SO4^2-': True,
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Hg2^2+': {
        'Cl-': False,
        'Br-': False,
        'I-': False,
        'OH-': False,
        'SO4^2-': True,
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Zn^2+': {
        'Cl-': True,
        'Br-': True,
        'I-': True,
        'OH-': False,
        'SO4^2-': True,
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    }
}

# Ksp values for common slightly soluble compounds
KSP_DATA = {
    'AgCl': {'ksp': 1.6e-10, 'molar_mass': 143.32, 'formula': 'AgCl', 'type': 'AB'},
    'AgBr': {'ksp': 5.0e-13, 'molar_mass': 187.77, 'formula': 'AgBr', 'type': 'AB'},
    'AgI': {'ksp': 8.3e-17, 'molar_mass': 234.77, 'formula': 'AgI', 'type': 'AB'},
    'BaSO4': {'ksp': 1.1e-10, 'molar_mass': 233.39, 'formula': 'BaSO4', 'type': 'AB'},
    'SrSO4': {'ksp': 3.4e-7, 'molar_mass': 183.68, 'formula': 'SrSO4', 'type': 'AB'},
    'CaCO3': {'ksp': 3.3e-9, 'molar_mass': 100.09, 'formula': 'CaCO3', 'type': 'AB'},
    'PbCl2': {'ksp': 1.9e-4, 'molar_mass': 278.11, 'formula': 'PbCl2', 'type': 'AB2'},
    'PbSO4': {'ksp': 2.5e-8, 'molar_mass': 303.26, 'formula': 'PbSO4', 'type': 'AB'},
    'Mg(OH)2': {'ksp': 5.6e-12, 'molar_mass': 58.32, 'formula': 'Mg(OH)2', 'type': 'AB2'},
    'Ca(OH)2': {'ksp': 5.5e-6, 'molar_mass': 74.09, 'formula': 'Ca(OH)2', 'type': 'AB2'},
    'CaF2': {'ksp': 3.9e-11, 'molar_mass': 78.07, 'formula': 'CaF2', 'type': 'AB2'},
    'Ag2CrO4': {'ksp': 1.1e-12, 'molar_mass': 331.73, 'formula': 'Ag2CrO4', 'type': 'A2B'},
    'Fe(OH)3': {'ksp': 2.8e-39, 'molar_mass': 106.87, 'formula': 'Fe(OH)3', 'type': 'AB3'},
    'Al(OH)3': {'ksp': 3.0e-34, 'molar_mass': 78.00, 'formula': 'Al(OH)3', 'type': 'AB3'}
}

# Additional information for reagents
REAGENTS = {
    'H2SO4': {'provides': 'SO4^2-'},
    'NaOH': {'provides': 'OH-'},
    'HCl': {'provides': 'Cl-'},
    'NaCl': {'provides': 'Cl-'},
    'KCl': {'provides': 'Cl-'},
    'AgNO3': {'provides': 'Ag+'},
    'BaCl2': {'provides': 'Ba^2+'}
}

# Predefined test scenarios
SCENARIOS = {
    "W20_8": {
        "description": "Scenario from problem W20_8",
        "cation_candidates": ["Ag+", "Ba^2+", "Pb^2+"],
        "precipitates_with": ["H2SO4"],
        "no_precipitate_with": ["NaOH"]
    },
    "S21_11": {
        "description": "Scenario from problem S21_11 - AgCl solubility",
        "compound": "AgCl",
        "ksp": 1.6e-10,
        "problem_type": "solubility"
    },
    "waste_water": {
        "description": "Waste water analysis problem",
        "cation_candidates": ["Na+", "Ag+", "Sr^2+", "K+"],
        "precipitates_with": ["H2SO4"],
        "no_precipitate_with": ["HCl"]
    }
}

def get_available_cations():
    """
    Get a list of all available cations in the database.
    
    Returns:
        list: List of cation names
    """
    return list(SOLUBILITY_DATA.keys())

def get_available_reagents():
    """
    Get a list of all available reagents in the database.
    
    Returns:
        list: List of reagent names
    """
    return list(REAGENTS.keys())

def get_available_compounds():
    """
    Get a list of all compounds with Ksp data.
    
    Returns:
        list: List of compound names
    """
    return list(KSP_DATA.keys())

def get_available_scenarios():
    """
    Get a list of all predefined scenarios.
    
    Returns:
        list: List of scenario IDs
    """
    return list(SCENARIOS.keys())

def get_scenario_data(scenario_id):
    """
    Get data for a specific scenario.
    
    Args:
        scenario_id (str): Identifier for the scenario
        
    Returns:
        dict: Scenario data or None if not found
    """
    return SCENARIOS.get(scenario_id, None)