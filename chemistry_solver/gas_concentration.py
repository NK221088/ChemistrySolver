"""
Module for gas concentration calculations and conversions.
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