"""
Redox Reaction Data Module - Stores all reference data for redox reactions

This module centralizes all reference data used for redox calculations:
- Standard reduction potentials
- Common redox reactions
- Metal information
- Oxidizing and reducing agents
"""

# Standard reduction potentials (V vs. SHE) at 25°C
STANDARD_REDUCTION_POTENTIALS = {
    # Non-metals
    "F2(g) + 2e- → 2F-(aq)": 2.87,
    "O3(g) + 2H+(aq) + 2e- → O2(g) + H2O(l)": 2.07,
    "H2O2(aq) + 2H+(aq) + 2e- → 2H2O(l)": 1.78,
    "MnO4-(aq) + 8H+(aq) + 5e- → Mn2+(aq) + 4H2O(l)": 1.51,
    "Cl2(g) + 2e- → 2Cl-(aq)": 1.36,
    "Cr2O7^2-(aq) + 14H+(aq) + 6e- → 2Cr3+(aq) + 7H2O(l)": 1.33,
    "O2(g) + 4H+(aq) + 4e- → 2H2O(l)": 1.23,
    "Br2(l) + 2e- → 2Br-(aq)": 1.09,
    "NO3-(aq) + 4H+(aq) + 3e- → NO(g) + 2H2O(l)": 0.96,
    "2H+(aq) + 2e- → H2(g)": 0.00,  # Reference half-reaction
    "S(s) + 2H+(aq) + 2e- → H2S(g)": 0.14,
    "Sn4+(aq) + 2e- → Sn2+(aq)": 0.15,
    "SO4^2-(aq) + 4H+(aq) + 2e- → H2SO3(aq) + H2O(l)": 0.17,
    "Cu2+(aq) + e- → Cu+(aq)": 0.16,
    "I2(s) + 2e- → 2I-(aq)": 0.54,
    
    # Metals
    "Ag+(aq) + e- → Ag(s)": 0.80,
    "Hg2^2+(aq) + 2e- → 2Hg(l)": 0.79,
    "Fe3+(aq) + e- → Fe2+(aq)": 0.77,
    "O2(g) + 2H+(aq) + 2e- → H2O2(aq)": 0.70,
    "MnO4-(aq) + e- → MnO4^2-(aq)": 0.56,
    "I2(s) + 2e- → 2I-(aq)": 0.54,
    "Cu+(aq) + e- → Cu(s)": 0.52,
    "I3-(aq) + 2e- → 3I-(aq)": 0.53,
    "Cu2+(aq) + 2e- → Cu(s)": 0.34,
    "AgCl(s) + e- → Ag(s) + Cl-(aq)": 0.22,
    "SO4^2-(aq) + 4H+(aq) + 2e- → SO2(g) + 2H2O(l)": 0.20,
    "Cu2+(aq) + e- → Cu+(aq)": 0.16,
    "Sn4+(aq) + 2e- → Sn2+(aq)": 0.15,
    "S(s) + 2H+(aq) + 2e- → H2S(g)": 0.14,
    "2H+(aq) + 2e- → H2(g)": 0.00,  # Reference half-reaction
    "Pb2+(aq) + 2e- → Pb(s)": -0.13,
    "Sn2+(aq) + 2e- → Sn(s)": -0.14,
    "Ni2+(aq) + 2e- → Ni(s)": -0.25,
    "Co2+(aq) + 2e- → Co(s)": -0.28,
    "PbSO4(s) + 2e- → Pb(s) + SO4^2-(aq)": -0.35,
    "Cd2+(aq) + 2e- → Cd(s)": -0.40,
    "Fe2+(aq) + 2e- → Fe(s)": -0.44,
    "Cr3+(aq) + e- → Cr2+(aq)": -0.41,
    "Cr3+(aq) + 3e- → Cr(s)": -0.74,
    "Zn2+(aq) + 2e- → Zn(s)": -0.76,
    "Cr2+(aq) + 2e- → Cr(s)": -0.91,
    "Mn2+(aq) + 2e- → Mn(s)": -1.18,
    "Al3+(aq) + 3e- → Al(s)": -1.66,
    "H2(g) + 2e- → 2H-(aq)": -2.25,
    "Mg2+(aq) + 2e- → Mg(s)": -2.37,
    "Na+(aq) + e- → Na(s)": -2.71,
    "Ca2+(aq) + 2e- → Ca(s)": -2.87,
    "Sr2+(aq) + 2e- → Sr(s)": -2.89,
    "Ba2+(aq) + 2e- → Ba(s)": -2.90,
    "K+(aq) + e- → K(s)": -2.93,
    "Li+(aq) + e- → Li(s)": -3.05,
}

# Common redox reactions with their balanced half-reactions and molar ratios
COMMON_REDOX_REACTIONS = {
    "SO2 + MnO4-": {
        "reaction": "SO2 + MnO4- + H2O → SO4^2- + Mn2+",
        "acidic": True,
        "oxidation_half": "SO2 + 2H2O → SO4^2- + 4H+ + 2e-",
        "reduction_half": "MnO4- + 8H+ + 5e- → Mn2+ + 4H2O",
        "balanced_equation": "5SO2 + 2MnO4- + 2H2O → 5SO4^2- + 2Mn2+ + 4H+",
        "molar_ratio": {"SO2:MnO4-": "5:2", "SO2:H+": "5:4", "SO2:Mn2+": "5:2", "SO2:SO4^2-": "1:1"}
    },
    "Fe2+ + MnO4-": {
        "reaction": "Fe2+ + MnO4- → Fe3+ + Mn2+",
        "acidic": True,
        "oxidation_half": "Fe2+ → Fe3+ + e-",
        "reduction_half": "MnO4- + 8H+ + 5e- → Mn2+ + 4H2O",
        "balanced_equation": "5Fe2+ + MnO4- + 8H+ → 5Fe3+ + Mn2+ + 4H2O",
        "molar_ratio": {"Fe2+:MnO4-": "5:1", "Fe2+:H+": "5:8", "Fe2+:Mn2+": "5:1"}
    },
    "SO2 + Cr2O7^2-": {
        "reaction": "SO2 + Cr2O7^2- → SO4^2- + Cr3+",
        "acidic": True,
        "oxidation_half": "SO2 + 2H2O → SO4^2- + 4H+ + 2e-",
        "reduction_half": "Cr2O7^2- + 14H+ + 6e- → 2Cr3+ + 7H2O",
        "balanced_equation": "3SO2 + Cr2O7^2- + 2H+ → 3SO4^2- + 2Cr3+ + H2O",
        "molar_ratio": {"SO2:Cr2O7^2-": "3:1", "SO2:H+": "3:2", "SO2:Cr3+": "3:2"}
    },
    "I- + MnO4-": {
        "reaction": "I- + MnO4- → I2 + Mn2+",
        "acidic": True,
        "oxidation_half": "2I- → I2 + 2e-",
        "reduction_half": "MnO4- + 8H+ + 5e- → Mn2+ + 4H2O",
        "balanced_equation": "10I- + 2MnO4- + 16H+ → 5I2 + 2Mn2+ + 8H2O",
        "molar_ratio": {"I-:MnO4-": "5:1", "I-:H+": "5:8", "I-:Mn2+": "5:1"}
    },
    "C2O4^2- + MnO4-": {
        "reaction": "C2O4^2- + MnO4- → CO2 + Mn2+",
        "acidic": True,
        "oxidation_half": "C2O4^2- → 2CO2 + 2e-",
        "reduction_half": "MnO4- + 8H+ + 5e- → Mn2+ + 4H2O",
        "balanced_equation": "5C2O4^2- + 2MnO4- + 16H+ → 10CO2 + 2Mn2+ + 8H2O",
        "molar_ratio": {"C2O4^2-:MnO4-": "5:2", "C2O4^2-:H+": "5:16", "C2O4^2-:Mn2+": "5:2"}
    },
    "Cu + NO3-": {
        "reaction": "Cu + NO3- → Cu2+ + NO",
        "acidic": True,
        "oxidation_half": "Cu → Cu2+ + 2e-",
        "reduction_half": "NO3- + 4H+ + 3e- → NO + 2H2O",
        "balanced_equation": "3Cu + 2NO3- + 8H+ → 3Cu2+ + 2NO + 4H2O",
        "molar_ratio": {"Cu:NO3-": "3:2", "Cu:H+": "3:8", "Cu:NO": "3:2"}
    },
    "Zn + H+": {
        "reaction": "Zn + H+ → Zn2+ + H2",
        "acidic": True,
        "oxidation_half": "Zn → Zn2+ + 2e-",
        "reduction_half": "2H+ + 2e- → H2",
        "balanced_equation": "Zn + 2H+ → Zn2+ + H2",
        "molar_ratio": {"Zn:H+": "1:2", "Zn:H2": "1:1"}
    }
}

# Common metal symbols for quick lookup
COMMON_METALS = {
    "Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", 
    "Al", "Ga", "In", "Sn", "Pb", "Bi", "Fe", "Co", "Ni", "Cu", 
    "Ag", "Au", "Zn", "Cd", "Hg", "Mn", "Tc", "Re", "Cr", "Mo", 
    "W", "V", "Nb", "Ta", "Ti", "Zr", "Hf"
}

# Activity series (from most active/reactive to least)
ACTIVITY_SERIES = [
    "Li", "K", "Ba", "Sr", "Ca", "Na", "Mg", "Al", "Mn", "Zn", 
    "Cr", "Fe", "Cd", "Co", "Ni", "Sn", "Pb", "H", "Cu", "Ag", 
    "Hg", "Pt", "Au"
]

# Common oxidizing and reducing agents with their half-reactions
COMMON_OXIDIZING_AGENTS = {
    "MnO4-": {
        "acidic": "MnO4- + 8H+ + 5e- → Mn2+ + 4H2O",
        "neutral": "MnO4- + 2H2O + 3e- → MnO2 + 4OH-",
        "basic": "MnO4- + 2H2O + 3e- → MnO2 + 4OH-"
    },
    "Cr2O7^2-": {
        "acidic": "Cr2O7^2- + 14H+ + 6e- → 2Cr3+ + 7H2O",
        "neutral": "Not commonly used in neutral conditions",
        "basic": "Not commonly used in basic conditions"
    },
    "H2O2": {
        "acidic": "H2O2 + 2H+ + 2e- → 2H2O",
        "neutral": "H2O2 + 2e- → 2OH-",
        "basic": "H2O2 + 2e- → 2OH-"
    },
    "O2": {
        "acidic": "O2 + 4H+ + 4e- → 2H2O",
        "neutral": "O2 + 2H2O + 4e- → 4OH-",
        "basic": "O2 + 2H2O + 4e- → 4OH-"
    },
    "Cl2": {
        "acidic": "Cl2 + 2e- → 2Cl-",
        "neutral": "Cl2 + 2e- → 2Cl-",
        "basic": "Cl2 + 2e- → 2Cl-"
    },
    "Br2": {
        "acidic": "Br2 + 2e- → 2Br-",
        "neutral": "Br2 + 2e- → 2Br-",
        "basic": "Br2 + 2e- → 2Br-"
    },
    "I2": {
        "acidic": "I2 + 2e- → 2I-",
        "neutral": "I2 + 2e- → 2I-",
        "basic": "I2 + 2e- → 2I-"
    },
    "NO3-": {
        "acidic": "NO3- + 4H+ + 3e- → NO + 2H2O",
        "neutral": "Not commonly used in neutral conditions",
        "basic": "Not commonly used in basic conditions"
    },
    "ClO-": {
        "acidic": "ClO- + 2H+ + 2e- → Cl- + H2O",
        "neutral": "ClO- + H2O + 2e- → Cl- + 2OH-",
        "basic": "ClO- + H2O + 2e- → Cl- + 2OH-"
    }
}

COMMON_REDUCING_AGENTS = {
    "Fe2+": {
        "acidic": "Fe2+ → Fe3+ + e-",
        "neutral": "Fe2+ → Fe3+ + e-",
        "basic": "Fe(OH)2 + OH- → Fe(OH)3 + e-"
    },
    "SO2": {
        "acidic": "SO2 + 2H2O → SO4^2- + 4H+ + 2e-",
        "neutral": "SO2 + 2OH- → SO4^2- + H2O + 2e-",
        "basic": "SO3^2- + 2OH- → SO4^2- + H2O + 2e-"
    },
    "H2S": {
        "acidic": "H2S → S + 2H+ + 2e-",
        "neutral": "H2S → S + 2H+ + 2e-",
        "basic": "HS- + OH- → S + H2O + 2e-"
    },
    "Zn": {
        "acidic": "Zn → Zn2+ + 2e-",
        "neutral": "Zn + 2OH- → Zn(OH)2 + 2e-",
        "basic": "Zn + 4OH- → [Zn(OH)4]^2- + 2e-"
    },
    "Fe": {
        "acidic": "Fe → Fe2+ + 2e-",
        "neutral": "Fe + 2OH- → Fe(OH)2 + 2e-",
        "basic": "Fe + 2OH- → Fe(OH)2 + 2e-"
    },
    "Cu": {
        "acidic": "Cu → Cu2+ + 2e-",
        "neutral": "Cu + 2OH- → Cu(OH)2 + 2e-",
        "basic": "Cu + 2OH- → Cu(OH)2 + 2e-"
    },
    "I-": {
        "acidic": "2I- → I2 + 2e-",
        "neutral": "2I- → I2 + 2e-",
        "basic": "2I- → I2 + 2e-"
    },
    "Sn2+": {
        "acidic": "Sn2+ → Sn4+ + 2e-",
        "neutral": "Sn(OH)2 + 2OH- → Sn(OH)4 + 2e-",
        "basic": "Sn(OH)2 + 2OH- → Sn(OH)4 + 2e-"
    }
}

# Common coupled redox half-reactions (oxidation and reduction pairs)
COMMON_REDOX_PAIRS = {
    "Zn/Cu2+": {
        "oxidation": "Zn → Zn2+ + 2e-",
        "reduction": "Cu2+ + 2e- → Cu",
        "net_reaction": "Zn + Cu2+ → Zn2+ + Cu",
        "e_cell": 1.10,  # V
        "favorable": True
    },
    "Fe/Cu2+": {
        "oxidation": "Fe → Fe2+ + 2e-",
        "reduction": "Cu2+ + 2e- → Cu",
        "net_reaction": "Fe + Cu2+ → Fe2+ + Cu",
        "e_cell": 0.78,  # V
        "favorable": True
    },
    "Zn/H+": {
        "oxidation": "Zn → Zn2+ + 2e-",
        "reduction": "2H+ + 2e- → H2",
        "net_reaction": "Zn + 2H+ → Zn2+ + H2",
        "e_cell": 0.76,  # V
        "favorable": True
    },
    "Cu/Ag+": {
        "oxidation": "Cu → Cu2+ + 2e-",
        "reduction": "2Ag+ + 2e- → 2Ag",
        "net_reaction": "Cu + 2Ag+ → Cu2+ + 2Ag",
        "e_cell": 0.46,  # V
        "favorable": True
    },
    "Fe2+/MnO4-": {
        "oxidation": "5Fe2+ → 5Fe3+ + 5e-",
        "reduction": "MnO4- + 8H+ + 5e- → Mn2+ + 4H2O",
        "net_reaction": "5Fe2+ + MnO4- + 8H+ → 5Fe3+ + Mn2+ + 4H2O",
        "e_cell": 0.74,  # V
        "favorable": True
    },
    "Cu/H+": {
        "oxidation": "Cu → Cu2+ + 2e-",
        "reduction": "2H+ + 2e- → H2",
        "net_reaction": "Cu + 2H+ → Cu2+ + H2",
        "e_cell": -0.34,  # V
        "favorable": False  # Copper won't displace hydrogen from acids
    }
}

# Standard electrode potentials in basic solution (pH = 14)
BASIC_REDUCTION_POTENTIALS = {
    "O2 + 2H2O + 4e- → 4OH-": 0.40,
    "MnO4- + 2H2O + 3e- → MnO2 + 4OH-": 0.59,
    "ClO- + H2O + 2e- → Cl- + 2OH-": 0.89,
    "2NO2- + 3H2O + 4e- → N2O + 6OH-": 0.15,
    "NO3- + H2O + 2e- → NO2- + 2OH-": -0.01
    # Add more as needed
}

# Standard Formation Enthalpies (kJ/mol) for calculating thermodynamic favorability
FORMATION_ENTHALPIES = {
    "H2O(l)": -285.8,
    "H2O(g)": -241.8,
    "CO2(g)": -393.5,
    "CO(g)": -110.5,
    "CH4(g)": -74.8,
    "C2H4(g)": 52.3,
    "C2H6(g)": -84.7,
    "C3H8(g)": -103.8,
    "NH3(g)": -46.1,
    "NO(g)": 90.3,
    "NO2(g)": 33.2,
    "HCl(g)": -92.3,
    "SO2(g)": -296.8,
    "SO3(g)": -395.7,
    "H2SO4(aq)": -814.0,
    "Fe2O3(s)": -824.2,
    "CaCO3(s)": -1207.6,
    "NaCl(s)": -411.2
    # Add more as needed
}

# Electrochemical series with reduction half-reactions and potentials
# Sorted from strongest oxidizing agents (top) to strongest reducing agents (bottom)
ELECTROCHEMICAL_SERIES = [
    {"half_reaction": "F2(g) + 2e- → 2F-(aq)", "potential": 2.87},
    {"half_reaction": "O3(g) + 2H+(aq) + 2e- → O2(g) + H2O(l)", "potential": 2.07},
    {"half_reaction": "H2O2(aq) + 2H+(aq) + 2e- → 2H2O(l)", "potential": 1.78},
    {"half_reaction": "MnO4-(aq) + 8H+(aq) + 5e- → Mn2+(aq) + 4H2O(l)", "potential": 1.51},
    {"half_reaction": "Cl2(g) + 2e- → 2Cl-(aq)", "potential": 1.36},
    {"half_reaction": "Cr2O7^2-(aq) + 14H+(aq) + 6e- → 2Cr3+(aq) + 7H2O(l)", "potential": 1.33},
    {"half_reaction": "O2(g) + 4H+(aq) + 4e- → 2H2O(l)", "potential": 1.23},
    {"half_reaction": "Br2(l) + 2e- → 2Br-(aq)", "potential": 1.09},
    {"half_reaction": "NO3-(aq) + 4H+(aq) + 3e- → NO(g) + 2H2O(l)", "potential": 0.96},
    {"half_reaction": "Ag+(aq) + e- → Ag(s)", "potential": 0.80},
    {"half_reaction": "Fe3+(aq) + e- → Fe2+(aq)", "potential": 0.77},
    {"half_reaction": "O2(g) + 2H+(aq) + 2e- → H2O2(aq)", "potential": 0.70},
    {"half_reaction": "I2(s) + 2e- → 2I-(aq)", "potential": 0.54},
    {"half_reaction": "Cu+(aq) + e- → Cu(s)", "potential": 0.52},
    {"half_reaction": "Cu2+(aq) + 2e- → Cu(s)", "potential": 0.34},
    {"half_reaction": "Sn4+(aq) + 2e- → Sn2+(aq)", "potential": 0.15},
    {"half_reaction": "2H+(aq) + 2e- → H2(g)", "potential": 0.00},  # Reference point
    {"half_reaction": "Pb2+(aq) + 2e- → Pb(s)", "potential": -0.13},
    {"half_reaction": "Sn2+(aq) + 2e- → Sn(s)", "potential": -0.14},
    {"half_reaction": "Ni2+(aq) + 2e- → Ni(s)", "potential": -0.25},
    {"half_reaction": "Co2+(aq) + 2e- → Co(s)", "potential": -0.28},
    {"half_reaction": "Cd2+(aq) + 2e- → Cd(s)", "potential": -0.40},
    {"half_reaction": "Fe2+(aq) + 2e- → Fe(s)", "potential": -0.44},
    {"half_reaction": "Cr3+(aq) + 3e- → Cr(s)", "potential": -0.74},
    {"half_reaction": "Zn2+(aq) + 2e- → Zn(s)", "potential": -0.76},
    {"half_reaction": "Mn2+(aq) + 2e- → Mn(s)", "potential": -1.18},
    {"half_reaction": "Al3+(aq) + 3e- → Al(s)", "potential": -1.66},
    {"half_reaction": "Mg2+(aq) + 2e- → Mg(s)", "potential": -2.37},
    {"half_reaction": "Na+(aq) + e- → Na(s)", "potential": -2.71},
    {"half_reaction": "Ca2+(aq) + 2e- → Ca(s)", "potential": -2.87},
    {"half_reaction": "K+(aq) + e- → K(s)", "potential": -2.93},
    {"half_reaction": "Li+(aq) + e- → Li(s)", "potential": -3.05}
]

# Helper functions for looking up data
def get_half_reaction_potential(half_reaction):
    """Look up the standard reduction potential for a half-reaction"""
    return STANDARD_REDUCTION_POTENTIALS.get(half_reaction)

def get_redox_pair(name):
    """Look up a common redox pair by name"""
    return COMMON_REDOX_PAIRS.get(name)

def get_common_redox_reaction(name):
    """Look up a common redox reaction by name"""
    return COMMON_REDOX_REACTIONS.get(name)

def get_oxidizing_agent(name, condition="acidic"):
    """Look up a common oxidizing agent half-reaction"""
    agent = COMMON_OXIDIZING_AGENTS.get(name)
    if agent:
        return agent.get(condition)
    return None

def get_reducing_agent(name, condition="acidic"):
    """Look up a common reducing agent half-reaction"""
    agent = COMMON_REDUCING_AGENTS.get(name)
    if agent:
        return agent.get(condition)
    return None