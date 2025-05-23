"""
Module for gas concentration calculations and conversions.
Extended with gas density calculations.
"""
from chemistry_solver.molar_mass import calculate_molar_mass

# Gas constant R (L·atm/(mol·K))
R = 0.08206

def ppm_to_mass(formula, ppm, volume_m3, temperature_celsius=25, pressure_atm=1):
    """
    Convert concentration in ppm to mass in grams.
    
    Args:
        formula (str): Chemical formula of the gas
        ppm (float): Concentration in parts per million
        volume_m3 (float): Volume in cubic meters
        temperature_celsius (float): Temperature in Celsius (default 25°C)
        pressure_atm (float): Pressure in atmospheres (default 1 atm)
    
    Returns:
        dict: Result dictionary with mass in grams and calculation details
    """
    try:
        # Get molar mass of the gas
        molar_mass_result = calculate_molar_mass(formula)
        if not molar_mass_result['success']:
            return {'success': False, 'error': f"Molar mass calculation failed: {molar_mass_result['error']}"}
        
        molar_mass = molar_mass_result['molar_mass']  # g/mol
        
        # Convert temperature to Kelvin
        temperature_kelvin = temperature_celsius + 273.15
        
        # Convert volume from m³ to liters
        volume_liters = volume_m3 * 1000
        
        # Calculate moles of air in the volume (using ideal gas law PV = nRT)
        moles_air = (pressure_atm * volume_liters) / (R * temperature_kelvin)
        
        # Calculate moles of gas at the given ppm
        # ppm is parts per million, so divide by 10^6
        moles_gas = moles_air * (ppm / 1_000_000)
        
        # Calculate mass of gas in grams
        mass_grams = moles_gas * molar_mass
        
        return {
            'success': True,
            'mass_grams': mass_grams,
            'formula': formula,
            'molar_mass': molar_mass,
            'volume_m3': volume_m3,
            'volume_liters': volume_liters,
            'temperature_celsius': temperature_celsius,
            'temperature_kelvin': temperature_kelvin,
            'pressure_atm': pressure_atm,
            'ppm': ppm,
            'moles_air': moles_air,
            'moles_gas': moles_gas
        }
    except Exception as e:
        return {'success': False, 'error': str(e)}

def mass_to_ppm(formula, mass_grams, volume_m3, temperature_celsius=25, pressure_atm=1):
    """
    Convert mass in grams to concentration in ppm.
    
    Args:
        formula (str): Chemical formula of the gas
        mass_grams (float): Mass in grams
        volume_m3 (float): Volume in cubic meters
        temperature_celsius (float): Temperature in Celsius (default 25°C)
        pressure_atm (float): Pressure in atmospheres (default 1 atm)
    
    Returns:
        dict: Result dictionary with concentration in ppm and calculation details
    """
    try:
        # Get molar mass of the gas
        molar_mass_result = calculate_molar_mass(formula)
        if not molar_mass_result['success']:
            return {'success': False, 'error': f"Molar mass calculation failed: {molar_mass_result['error']}"}
        
        molar_mass = molar_mass_result['molar_mass']  # g/mol
        
        # Convert temperature to Kelvin
        temperature_kelvin = temperature_celsius + 273.15
        
        # Convert volume from m³ to liters
        volume_liters = volume_m3 * 1000
        
        # Calculate moles of air in the volume (using ideal gas law PV = nRT)
        moles_air = (pressure_atm * volume_liters) / (R * temperature_kelvin)
        
        # Calculate moles of gas from mass
        moles_gas = mass_grams / molar_mass
        
        # Calculate ppm (parts per million)
        ppm = (moles_gas / moles_air) * 1_000_000
        
        return {
            'success': True,
            'ppm': ppm,
            'formula': formula,
            'molar_mass': molar_mass,
            'volume_m3': volume_m3,
            'volume_liters': volume_liters,
            'temperature_celsius': temperature_celsius,
            'temperature_kelvin': temperature_kelvin,
            'pressure_atm': pressure_atm,
            'mass_grams': mass_grams,
            'moles_air': moles_air,
            'moles_gas': moles_gas
        }
    except Exception as e:
        return {'success': False, 'error': str(e)}

def calculate_gas_density(formula, temperature_celsius=25, pressure_atm=1):
    """
    Calculate the density of a gas at given conditions using the ideal gas law.
    
    Args:
        formula (str): Chemical formula of the gas
        temperature_celsius (float): Temperature in Celsius (default 25°C)
        pressure_atm (float): Pressure in atmospheres (default 1 atm)
    
    Returns:
        dict: Result dictionary with density in g/L and calculation details
    """
    try:
        # Get molar mass of the gas
        molar_mass_result = calculate_molar_mass(formula)
        if not molar_mass_result['success']:
            return {'success': False, 'error': f"Molar mass calculation failed: {molar_mass_result['error']}"}
        
        molar_mass = molar_mass_result['molar_mass']  # g/mol
        
        # Convert temperature to Kelvin
        temperature_kelvin = temperature_celsius + 273.15
        
        # Calculate density using: density = (P × M) / (R × T)
        # Where P = pressure (atm), M = molar mass (g/mol), R = gas constant, T = temperature (K)
        density_g_per_L = (pressure_atm * molar_mass) / (R * temperature_kelvin)
        
        return {
            'success': True,
            'density_g_per_L': density_g_per_L,
            'formula': formula,
            'molar_mass': molar_mass,
            'temperature_celsius': temperature_celsius,
            'temperature_kelvin': temperature_kelvin,
            'pressure_atm': pressure_atm
        }
    except Exception as e:
        return {'success': False, 'error': str(e)}

def identify_gas_by_density(candidate_formulas, target_density, temperature_celsius=25, pressure_atm=1, tolerance=0.01):
    """
    Identify which gas formula matches a target density within a given tolerance.
    
    Args:
        candidate_formulas (list): List of chemical formulas to test
        target_density (float): Target density in g/L
        temperature_celsius (float): Temperature in Celsius (default 25°C)
        pressure_atm (float): Pressure in atmospheres (default 1 atm)
        tolerance (float): Acceptable difference from target density (default 0.01 g/L)
    
    Returns:
        dict: Result dictionary with matching formula(s) and all calculated densities
    """
    try:
        results = []
        matches = []
        
        for formula in candidate_formulas:
            density_result = calculate_gas_density(formula, temperature_celsius, pressure_atm)
            
            if density_result['success']:
                calculated_density = density_result['density_g_per_L']
                difference = abs(calculated_density - target_density)
                
                result_entry = {
                    'formula': formula,
                    'calculated_density': calculated_density,
                    'difference': difference,
                    'molar_mass': density_result['molar_mass']
                }
                results.append(result_entry)
                
                if difference <= tolerance:
                    matches.append(result_entry)
            else:
                results.append({
                    'formula': formula,
                    'error': density_result['error']
                })
        
        # Sort results by difference (closest match first)
        valid_results = [r for r in results if 'calculated_density' in r]
        valid_results.sort(key=lambda x: x['difference'])
        
        return {
            'success': True,
            'target_density': target_density,
            'temperature_celsius': temperature_celsius,
            'temperature_kelvin': temperature_celsius + 273.15,
            'pressure_atm': pressure_atm,
            'tolerance': tolerance,
            'matches': matches,
            'all_results': results,
            'best_match': valid_results[0] if valid_results else None
        }
    except Exception as e:
        return {'success': False, 'error': str(e)}

def solve_density_problem(candidate_formulas, target_density, temperature_celsius, pressure_atm):
    """
    Solve a gas density identification problem like the one in your chemistry question.
    
    Args:
        candidate_formulas (list): List of chemical formulas to test
        target_density (float): Target density in g/L
        temperature_celsius (float): Temperature in Celsius
        pressure_atm (float): Pressure in atmospheres
    
    Returns:
        dict: Detailed solution with step-by-step calculations
    """
    print(f"\n=== SOLVING GAS DENSITY IDENTIFICATION PROBLEM ===")
    print(f"Target density: {target_density} g/L at {temperature_celsius}°C and {pressure_atm} atm")
    print(f"Candidate compounds: {', '.join(candidate_formulas)}")
    print(f"\nUsing ideal gas law: density = (P × M) / (R × T)")
    print(f"Where P = {pressure_atm} atm, R = {R} L·atm/(mol·K), T = {temperature_celsius + 273.15} K")
    
    result = identify_gas_by_density(candidate_formulas, target_density, temperature_celsius, pressure_atm, tolerance=0.1)
    
    if result['success']:
        print(f"\n=== CALCULATIONS FOR EACH COMPOUND ===")
        for entry in result['all_results']:
            if 'calculated_density' in entry:
                print(f"\n{entry['formula']}:")
                print(f"  Molar mass = {entry['molar_mass']:.4f} g/mol")
                print(f"  Calculated density = ({pressure_atm} × {entry['molar_mass']:.4f}) / ({R} × {temperature_celsius + 273.15}) = {entry['calculated_density']:.4f} g/L")
                print(f"  Difference from target = |{entry['calculated_density']:.4f} - {target_density}| = {entry['difference']:.4f} g/L")
            else:
                print(f"\n{entry['formula']}: Error - {entry['error']}")
        
        print(f"\n=== RESULTS ===")
        if result['best_match']:
            best = result['best_match']
            print(f"Best match: {best['formula']} (difference: {best['difference']:.4f} g/L)")
            
            if result['matches']:
                print(f"Exact matches within tolerance: {[m['formula'] for m in result['matches']]}")
            else:
                print("No exact matches within default tolerance, but closest match shown above.")
        else:
            print("No valid calculations could be performed.")
    
    return result