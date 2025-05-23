"""
Colligative Properties Constants

This module contains solvent constants and other data used for colligative property calculations.
Includes freezing point depression constants (Kf), boiling point elevation constants (Kb),
and normal freezing/boiling points for various solvents.
"""

# Solvent constants for colligative property calculations
# Each entry contains:
# - Kf: Freezing point depression constant (°C·kg/mol)
# - Kb: Boiling point elevation constant (°C·kg/mol)  
# - freezing_point: Normal freezing point (°C)
# - boiling_point: Normal boiling point (°C)
SOLVENT_CONSTANTS = {
    "water": {
        "Kf": 1.86, 
        "Kb": 0.512, 
        "freezing_point": 0.0, 
        "boiling_point": 100.0
    },
    "benzene": {
        "Kf": 5.12, 
        "Kb": 2.53, 
        "freezing_point": 5.5, 
        "boiling_point": 80.1
    },
    "cyclohexane": {
        "Kf": 20.0, 
        "Kb": 2.79, 
        "freezing_point": 6.6, 
        "boiling_point": 80.7
    },
    "camphor": {
        "Kf": 40.0, 
        "Kb": 5.95, 
        "freezing_point": 179.8, 
        "boiling_point": 209.0
    },
    "acetic_acid": {
        "Kf": 3.90, 
        "Kb": 3.07, 
        "freezing_point": 16.6, 
        "boiling_point": 118.0
    },
    "naphthalene": {
        "Kf": 6.94, 
        "Kb": 5.80, 
        "freezing_point": 80.2, 
        "boiling_point": 218.0
    },
    "chloroform": {
        "Kf": 4.68, 
        "Kb": 3.63, 
        "freezing_point": -63.5, 
        "boiling_point": 61.2
    },
    "ethanol": {
        "Kf": 1.99, 
        "Kb": 1.22, 
        "freezing_point": -114.1, 
        "boiling_point": 78.4
    }
}

# Gas constant for osmotic pressure calculations
R_OSMOTIC = 0.08206  # L·atm/(mol·K)

# Conversion factors
GRAMS_TO_KG = 1000
CELSIUS_TO_KELVIN = 273.15