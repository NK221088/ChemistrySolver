"""
Module for solving thermodynamics and heat transfer problems in chemistry.
"""

import math
from typing import Dict, List, Tuple, Any

# Gas constant in different units
R_IDEAL_GAS = {
    'J/(mol·K)': 8.314,
    'kJ/(mol·K)': 0.008314,
    'L·atm/(mol·K)': 0.0821
}

def calculate_heat(mass, specific_heat, delta_t):
    """
    Calculate heat energy transfer using q = m × c × ΔT.
    
    Args:
        mass (float): Mass in grams
        specific_heat (float): Specific heat capacity in J/(g·K) or J/(g·°C)
        delta_t (float): Temperature change in Celsius or Kelvin
        
    Returns:
        float: Heat energy in Joules
    """
    return mass * specific_heat * delta_t

def calculate_temperature_change(heat, mass, specific_heat):
    """
    Calculate temperature change using ΔT = q / (m × c).
    
    Args:
        heat (float): Heat energy in Joules
        mass (float): Mass in grams
        specific_heat (float): Specific heat capacity in J/(g·K) or J/(g·°C)
        
    Returns:
        float: Temperature change in Celsius or Kelvin
    """
    return heat / (mass * specific_heat)

def calculate_final_temperature(initial_temp, delta_t):
    """
    Calculate final temperature using Tf = Ti + ΔT.
    
    Args:
        initial_temp (float): Initial temperature in Celsius
        delta_t (float): Temperature change in Celsius
        
    Returns:
        float: Final temperature in Celsius
    """
    return initial_temp + delta_t

def calculate_molar_heat(mass, molar_mass, molar_heat_capacity, delta_t):
    """
    Calculate heat energy transfer using molar heat capacity.
    
    Args:
        mass (float): Mass in grams
        molar_mass (float): Molar mass in g/mol
        molar_heat_capacity (float): Molar heat capacity in J/(mol·K) or J/(mol·°C)
        delta_t (float): Temperature change in Celsius or Kelvin
        
    Returns:
        float: Heat energy in Joules
    """
    moles = mass / molar_mass
    return moles * molar_heat_capacity * delta_t

def solve_thermal_equilibrium(substances):
    """
    Solve for the final temperature when multiple substances reach thermal equilibrium.
    
    Args:
        substances (list): List of dictionaries with the following keys:
            - mass (float): Mass in grams
            - specific_heat (float): Specific heat capacity in J/(g·K) or J/(g·°C)
            - initial_temp (float): Initial temperature in Celsius
            
    Returns:
        float: Final equilibrium temperature in Celsius
    """
    numerator = 0
    denominator = 0
    
    for substance in substances:
        m = substance['mass']
        c = substance['specific_heat']
        T = substance['initial_temp']
        
        numerator += m * c * T
        denominator += m * c
    
    return numerator / denominator

def solve_thermal_equilibrium_with_molar_heat(substances):
    """
    Solve for the final temperature when multiple substances reach thermal equilibrium,
    using molar heat capacities.
    
    Args:
        substances (list): List of dictionaries with the following keys:
            - mass (float): Mass in grams
            - molar_mass (float): Molar mass in g/mol
            - molar_heat_capacity (float): Molar heat capacity in J/(mol·K) or J/(mol·°C)
            - initial_temp (float): Initial temperature in Celsius
            
    Returns:
        float: Final equilibrium temperature in Celsius
    """
    numerator = 0
    denominator = 0
    
    for substance in substances:
        m = substance['mass']
        mm = substance['molar_mass']
        mhc = substance['molar_heat_capacity']
        T = substance['initial_temp']
        
        moles = m / mm
        numerator += moles * mhc * T
        denominator += moles * mhc
    
    return numerator / denominator

def handle_heat_transfer_problem(mass1, specific_heat1, initial_temp1, 
                                mass2, specific_heat2, initial_temp2):
    """
    Solve a problem where two substances exchange heat until they reach equilibrium.
    
    Args:
        mass1 (float): Mass of substance 1 in grams
        specific_heat1 (float): Specific heat capacity of substance 1 in J/(g·K) or J/(g·°C)
        initial_temp1 (float): Initial temperature of substance 1 in Celsius
        mass2 (float): Mass of substance 2 in grams
        specific_heat2 (float): Specific heat capacity of substance 2 in J/(g·K) or J/(g·°C)
        initial_temp2 (float): Initial temperature of substance 2 in Celsius
        
    Returns:
        dict: Contains final_temp and steps
    """
    # Calculate final temperature
    substances = [
        {'mass': mass1, 'specific_heat': specific_heat1, 'initial_temp': initial_temp1},
        {'mass': mass2, 'specific_heat': specific_heat2, 'initial_temp': initial_temp2}
    ]
    
    final_temp = solve_thermal_equilibrium(substances)
    
    # Calculate heat transferred
    delta_t1 = final_temp - initial_temp1
    delta_t2 = final_temp - initial_temp2
    
    q1 = calculate_heat(mass1, specific_heat1, delta_t1)
    q2 = calculate_heat(mass2, specific_heat2, delta_t2)
    
    # Generate solution steps
    steps = [
        f"1. Calculate the final equilibrium temperature:",
        f"   - Using conservation of energy: q₁ + q₂ = 0",
        f"   - m₁ × c₁ × (Tf - T₁) + m₂ × c₂ × (Tf - T₂) = 0",
        f"   - Tf = (m₁ × c₁ × T₁ + m₂ × c₂ × T₂) / (m₁ × c₁ + m₂ × c₂)",
        f"   - Tf = ({mass1} × {specific_heat1} × {initial_temp1} + {mass2} × {specific_heat2} × {initial_temp2}) / ({mass1} × {specific_heat1} + {mass2} × {specific_heat2})",
        f"   - Tf = {final_temp:.2f} °C",
        f"",
        f"2. Calculate heat transferred by substance 1:",
        f"   - q₁ = m₁ × c₁ × (Tf - T₁)",
        f"   - q₁ = {mass1} × {specific_heat1} × ({final_temp:.2f} - {initial_temp1})",
        f"   - q₁ = {q1:.2f} J",
        f"",
        f"3. Calculate heat transferred by substance 2:",
        f"   - q₂ = m₂ × c₂ × (Tf - T₂)",
        f"   - q₂ = {mass2} × {specific_heat2} × ({final_temp:.2f} - {initial_temp2})",
        f"   - q₂ = {q2:.2f} J",
        f"",
        f"4. Verify conservation of energy:",
        f"   - q₁ + q₂ = {q1:.2f} + {q2:.2f} = {q1 + q2:.2f} J (approximately 0 due to rounding)"
    ]
    
    return {
        "final_temp": final_temp,
        "heat_transferred_1": q1,
        "heat_transferred_2": q2,
        "steps": steps
    }

def handle_heat_transfer_with_molar_heat(mass1, molar_mass1, molar_heat_capacity1, initial_temp1, 
                                        mass2, molar_mass2, molar_heat_capacity2, initial_temp2):
    """
    Solve a problem where two substances exchange heat until they reach equilibrium,
    using molar heat capacities.
    
    Args:
        mass1 (float): Mass of substance 1 in grams
        molar_mass1 (float): Molar mass of substance 1 in g/mol
        molar_heat_capacity1 (float): Molar heat capacity of substance 1 in J/(mol·K) or J/(mol·°C)
        initial_temp1 (float): Initial temperature of substance 1 in Celsius
        mass2 (float): Mass of substance 2 in grams
        molar_mass2 (float): Molar mass of substance 2 in g/mol
        molar_heat_capacity2 (float): Molar heat capacity of substance 2 in J/(mol·K) or J/(mol·°C)
        initial_temp2 (float): Initial temperature of substance 2 in Celsius
        
    Returns:
        dict: Contains final_temp and steps
    """
    # Calculate final temperature
    substances = [
        {'mass': mass1, 'molar_mass': molar_mass1, 'molar_heat_capacity': molar_heat_capacity1, 'initial_temp': initial_temp1},
        {'mass': mass2, 'molar_mass': molar_mass2, 'molar_heat_capacity': molar_heat_capacity2, 'initial_temp': initial_temp2}
    ]
    
    final_temp = solve_thermal_equilibrium_with_molar_heat(substances)
    
    # Calculate moles
    moles1 = mass1 / molar_mass1
    moles2 = mass2 / molar_mass2
    
    # Calculate heat transferred
    delta_t1 = final_temp - initial_temp1
    delta_t2 = final_temp - initial_temp2
    
    q1 = moles1 * molar_heat_capacity1 * delta_t1
    q2 = moles2 * molar_heat_capacity2 * delta_t2
    
    # Generate solution steps
    steps = [
        f"1. Calculate moles for each substance:",
        f"   - Moles of substance 1: n₁ = m₁ / M₁ = {mass1} g / {molar_mass1} g/mol = {moles1:.4f} mol",
        f"   - Moles of substance 2: n₂ = m₂ / M₂ = {mass2} g / {molar_mass2} g/mol = {moles2:.4f} mol",
        f"",
        f"2. Calculate the final equilibrium temperature:",
        f"   - Using conservation of energy: q₁ + q₂ = 0",
        f"   - n₁ × C₁ × (Tf - T₁) + n₂ × C₂ × (Tf - T₂) = 0",
        f"   - Tf = (n₁ × C₁ × T₁ + n₂ × C₂ × T₂) / (n₁ × C₁ + n₂ × C₂)",
        f"   - Tf = ({moles1:.4f} × {molar_heat_capacity1} × {initial_temp1} + {moles2:.4f} × {molar_heat_capacity2} × {initial_temp2}) / ({moles1:.4f} × {molar_heat_capacity1} + {moles2:.4f} × {molar_heat_capacity2})",
        f"   - Tf = {final_temp:.2f} °C",
        f"",
        f"3. Calculate heat transferred by substance 1:",
        f"   - q₁ = n₁ × C₁ × (Tf - T₁)",
        f"   - q₁ = {moles1:.4f} × {molar_heat_capacity1} × ({final_temp:.2f} - {initial_temp1})",
        f"   - q₁ = {q1:.2f} J",
        f"",
        f"4. Calculate heat transferred by substance 2:",
        f"   - q₂ = n₂ × C₂ × (Tf - T₂)",
        f"   - q₂ = {moles2:.4f} × {molar_heat_capacity2} × ({final_temp:.2f} - {initial_temp2})",
        f"   - q₂ = {q2:.2f} J",
        f"",
        f"5. Verify conservation of energy:",
        f"   - q₁ + q₂ = {q1:.2f} + {q2:.2f} = {q1 + q2:.2f} J (approximately 0 due to rounding)"
    ]
    
    return {
        "final_temp": final_temp,
        "heat_transferred_1": q1,
        "heat_transferred_2": q2,
        "steps": steps
    }

def solve_mixture_problem(substances):
    """
    Solve a general thermal equilibrium problem with multiple substances.
    
    Args:
        substances (list): List of dictionaries, each containing:
            - name (str): Substance name
            - mass (float): Mass in g
            - specific_heat (float, optional): Specific heat capacity in J/(g·K)
            - molar_mass (float, optional): Molar mass in g/mol
            - molar_heat_capacity (float, optional): Molar heat capacity in J/(mol·K)
            - initial_temp (float): Initial temperature in °C
            
    Returns:
        dict: Contains final_temp, heat_transfers, and steps
    """
    # Prepare substances with consistent heat capacity parameters
    processed_substances = []
    
    for s in substances:
        substance = s.copy()
        
        # If molar heat capacity is provided, convert to specific heat
        if 'molar_heat_capacity' in substance and 'molar_mass' in substance:
            if 'specific_heat' not in substance:
                substance['specific_heat'] = substance['molar_heat_capacity'] / substance['molar_mass']
        
        # Make sure we have specific heat capacity
        if 'specific_heat' not in substance:
            raise ValueError(f"Missing specific heat for {substance['name']}")
        
        processed_substances.append(substance)
    
    # Calculate final temperature
    final_temp = solve_thermal_equilibrium(processed_substances)
    
    # Calculate heat transferred for each substance
    heat_transfers = []
    
    for substance in processed_substances:
        delta_t = final_temp - substance['initial_temp']
        heat = calculate_heat(substance['mass'], substance['specific_heat'], delta_t)
        heat_transfers.append({
            'name': substance['name'],
            'heat': heat,
            'delta_t': delta_t
        })
    
    # Generate solution steps
    steps = [
        f"1. Calculate the final equilibrium temperature:",
        f"   - Using conservation of energy: Σ q = 0",
        f"   - Σ (m × c × (Tf - Ti)) = 0",
        f"   - Tf = Σ(m × c × Ti) / Σ(m × c)"
    ]
    
    # Add details for each substance
    for i, substance in enumerate(processed_substances):
        steps.append(f"   - For {substance['name']}: m = {substance['mass']} g, c = {substance['specific_heat']} J/(g·K), Ti = {substance['initial_temp']} °C")
    
    steps.append(f"   - Final temperature: Tf = {final_temp:.2f} °C")
    steps.append("")
    
    # Add heat transfer calculations
    steps.append(f"2. Calculate heat transferred for each substance:")
    
    for transfer in heat_transfers:
        steps.append(f"   - {transfer['name']}: q = m × c × (Tf - Ti) = {transfer['heat']:.2f} J (ΔT = {transfer['delta_t']:.2f} °C)")
    
    steps.append("")
    steps.append(f"3. Verify conservation of energy:")
    total_heat = sum(t['heat'] for t in heat_transfers)
    steps.append(f"   - Sum of all heat transfers: {total_heat:.2f} J (approximately 0 due to rounding)")
    
    return {
        "final_temp": final_temp,
        "heat_transfers": heat_transfers,
        "steps": steps
    }

# Clausius-Clapeyron equation functions for pressure-temperature relationships

def calculate_boiling_point_with_pressure(
    normal_boiling_point_c: float, 
    heat_of_vaporization: float,  # in kJ/mol
    initial_pressure: float,      # in atm
    final_pressure: float,        # in atm
) -> Dict[str, Any]:
    """
    Calculate the new boiling point when pressure changes using the Clausius-Clapeyron equation.
    
    Parameters:
        normal_boiling_point_c (float): The normal boiling point in degrees Celsius
        heat_of_vaporization (float): The heat of vaporization in kJ/mol
        initial_pressure (float): The initial pressure in atm
        final_pressure (float): The final pressure in atm
        
    Returns:
        Dict[str, Any]: Dictionary containing results and solution steps
    """
    # Convert boiling point to Kelvin
    normal_boiling_point_k = normal_boiling_point_c + 273.15
    
    # Gas constant in kJ/(mol·K)
    R = R_IDEAL_GAS['kJ/(mol·K)']
    
    # Calculate the new boiling point using the integrated Clausius-Clapeyron equation
    # ln(P2/P1) = -(ΔHvap/R) * (1/T2 - 1/T1)
    # Rearranging: 1/T2 = 1/T1 - (R/ΔHvap) * ln(P2/P1)
    
    # Calculate 1/T2
    inv_t2 = (1/normal_boiling_point_k) - (R/heat_of_vaporization) * math.log(final_pressure/initial_pressure)
    
    # Calculate T2
    new_boiling_point_k = 1/inv_t2
    new_boiling_point_c = new_boiling_point_k - 273.15
    
    # Build solution steps
    steps = [
        f"Step 1: Convert the normal boiling point to Kelvin:",
        f"T₁ = {normal_boiling_point_c} °C + 273.15 = {normal_boiling_point_k:.2f} K",
        f"",
        f"Step 2: Use the Clausius-Clapeyron equation to find the new boiling point:",
        f"ln(P₂/P₁) = -(ΔHvap/R) × (1/T₂ - 1/T₁)",
        f"",
        f"Rearranging for 1/T₂:",
        f"1/T₂ = 1/T₁ - (R/ΔHvap) × ln(P₂/P₁)",
        f"",
        f"Step 3: Substitute the values:",
        f"1/T₂ = 1/{normal_boiling_point_k:.2f} K - ({R:.6f} kJ/(mol·K)/{heat_of_vaporization} kJ/mol) × ln({final_pressure}/{initial_pressure})",
        f"1/T₂ = {1/normal_boiling_point_k:.6f} K⁻¹ - {R/heat_of_vaporization:.6f} × {math.log(final_pressure/initial_pressure):.6f}",
        f"1/T₂ = {1/normal_boiling_point_k:.6f} K⁻¹ - {(R/heat_of_vaporization) * math.log(final_pressure/initial_pressure):.6f} K⁻¹",
        f"1/T₂ = {inv_t2:.6f} K⁻¹",
        f"",
        f"Step 4: Calculate T₂:",
        f"T₂ = 1/({inv_t2:.6f} K⁻¹) = {new_boiling_point_k:.2f} K",
        f"",
        f"Step 5: Convert back to Celsius:",
        f"T₂ = {new_boiling_point_k:.2f} K - 273.15 = {new_boiling_point_c:.2f} °C"
    ]
    
    return {
        "normal_boiling_point_c": normal_boiling_point_c,
        "normal_boiling_point_k": normal_boiling_point_k,
        "heat_of_vaporization": heat_of_vaporization,
        "initial_pressure": initial_pressure,
        "final_pressure": final_pressure,
        "new_boiling_point_k": new_boiling_point_k,
        "new_boiling_point_c": new_boiling_point_c,
        "steps": steps
    }

def calculate_pressure_with_temperature(
    normal_boiling_point_c: float,
    heat_of_vaporization: float,  # in kJ/mol
    initial_pressure: float,      # in atm
    final_temperature_c: float    # in °C
) -> Dict[str, Any]:
    """
    Calculate the vapor pressure at a given temperature using the Clausius-Clapeyron equation.
    
    Parameters:
        normal_boiling_point_c (float): The normal boiling point in degrees Celsius
        heat_of_vaporization (float): The heat of vaporization in kJ/mol
        initial_pressure (float): The initial pressure in atm (typically 1 atm)
        final_temperature_c (float): The temperature at which to calculate the vapor pressure, in °C
        
    Returns:
        Dict[str, Any]: Dictionary containing results and solution steps
    """
    # Convert temperatures to Kelvin
    normal_boiling_point_k = normal_boiling_point_c + 273.15
    final_temperature_k = final_temperature_c + 273.15
    
    # Gas constant in kJ/(mol·K)
    R = R_IDEAL_GAS['kJ/(mol·K)']
    
    # Calculate the new pressure using the integrated Clausius-Clapeyron equation
    # ln(P2/P1) = -(ΔHvap/R) * (1/T2 - 1/T1)
    
    exponent = -(heat_of_vaporization/R) * (1/final_temperature_k - 1/normal_boiling_point_k)
    final_pressure = initial_pressure * math.exp(exponent)
    
    # Build solution steps
    steps = [
        f"Step 1: Convert temperatures to Kelvin:",
        f"T₁ = {normal_boiling_point_c} °C + 273.15 = {normal_boiling_point_k:.2f} K",
        f"T₂ = {final_temperature_c} °C + 273.15 = {final_temperature_k:.2f} K",
        f"",
        f"Step 2: Use the Clausius-Clapeyron equation to find the new pressure:",
        f"ln(P₂/P₁) = -(ΔHvap/R) × (1/T₂ - 1/T₁)",
        f"",
        f"Step 3: Substitute the values:",
        f"ln(P₂/{initial_pressure}) = -({heat_of_vaporization} kJ/mol/{R:.6f} kJ/(mol·K)) × (1/{final_temperature_k:.2f} K - 1/{normal_boiling_point_k:.2f} K)",
        f"ln(P₂/{initial_pressure}) = -{heat_of_vaporization/R:.2f} × ({1/final_temperature_k:.6f} K⁻¹ - {1/normal_boiling_point_k:.6f} K⁻¹)",
        f"ln(P₂/{initial_pressure}) = -{heat_of_vaporization/R:.2f} × {1/final_temperature_k - 1/normal_boiling_point_k:.6f} K⁻¹",
        f"ln(P₂/{initial_pressure}) = {exponent:.6f}",
        f"",
        f"Step 4: Calculate P₂:",
        f"P₂ = {initial_pressure} × e^({exponent:.6f})",
        f"P₂ = {initial_pressure} × {math.exp(exponent):.6f}",
        f"P₂ = {final_pressure:.6f} atm"
    ]
    
    return {
        "normal_boiling_point_c": normal_boiling_point_c,
        "normal_boiling_point_k": normal_boiling_point_k,
        "heat_of_vaporization": heat_of_vaporization,
        "initial_pressure": initial_pressure,
        "final_temperature_c": final_temperature_c,
        "final_temperature_k": final_temperature_k,
        "final_pressure": final_pressure,
        "steps": steps
    }

def calculate_heat_of_vaporization(
    temp1_c: float,
    pressure1: float,  # in atm
    temp2_c: float,
    pressure2: float   # in atm
) -> Dict[str, Any]:
    """
    Calculate the heat of vaporization using two pressure-temperature data points.
    
    Parameters:
        temp1_c (float): First temperature in degrees Celsius
        pressure1 (float): First pressure in atm
        temp2_c (float): Second temperature in degrees Celsius
        pressure2 (float): Second pressure in atm
        
    Returns:
        Dict[str, Any]: Dictionary containing results and solution steps
    """
    # Convert temperatures to Kelvin
    temp1_k = temp1_c + 273.15
    temp2_k = temp2_c + 273.15
    
    # Gas constant in kJ/(mol·K)
    R = R_IDEAL_GAS['kJ/(mol·K)']
    
    # Calculate heat of vaporization using the Clausius-Clapeyron equation
    # ln(P2/P1) = -(ΔHvap/R) * (1/T2 - 1/T1)
    # Rearranging: ΔHvap = -R * ln(P2/P1) / (1/T2 - 1/T1)
    
    heat_of_vaporization = -R * math.log(pressure2/pressure1) / (1/temp2_k - 1/temp1_k)
    
    # Build solution steps
    steps = [
        f"Step 1: Convert temperatures to Kelvin:",
        f"T₁ = {temp1_c} °C + 273.15 = {temp1_k:.2f} K",
        f"T₂ = {temp2_c} °C + 273.15 = {temp2_k:.2f} K",
        f"",
        f"Step 2: Rearrange the Clausius-Clapeyron equation to solve for the heat of vaporization:",
        f"ln(P₂/P₁) = -(ΔHvap/R) × (1/T₂ - 1/T₁)",
        f"ΔHvap = -R × ln(P₂/P₁) / (1/T₂ - 1/T₁)",
        f"",
        f"Step 3: Substitute the values:",
        f"ΔHvap = -{R:.6f} kJ/(mol·K) × ln({pressure2}/{pressure1}) / (1/{temp2_k:.2f} K - 1/{temp1_k:.2f} K)",
        f"ΔHvap = -{R:.6f} kJ/(mol·K) × {math.log(pressure2/pressure1):.6f} / ({1/temp2_k:.6f} K⁻¹ - {1/temp1_k:.6f} K⁻¹)",
        f"ΔHvap = -{R:.6f} kJ/(mol·K) × {math.log(pressure2/pressure1):.6f} / {1/temp2_k - 1/temp1_k:.6f} K⁻¹",
        f"ΔHvap = {heat_of_vaporization:.2f} kJ/mol"
    ]
    
    return {
        "temp1_c": temp1_c,
        "temp1_k": temp1_k,
        "pressure1": pressure1,
        "temp2_c": temp2_c,
        "temp2_k": temp2_k,
        "pressure2": pressure2,
        "heat_of_vaporization": heat_of_vaporization,
        "steps": steps
    }