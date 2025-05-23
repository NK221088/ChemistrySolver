"""
Chemical Kinetics Problem Solver Module

This module provides functions to solve various chemical kinetics problems,
including first-order kinetics and Arrhenius equation problems.
"""
import math
from fractions import Fraction

# ===== FIRST-ORDER KINETICS FUNCTIONS =====

def calculate_rate_constant_from_half_life(half_life):
    """
    Calculate rate constant (k) from half-life for first-order reactions.
    
    Args:
        half_life (float): Half-life value
    
    Returns:
        float: Rate constant k
    """
    if half_life <= 0:
        raise ValueError("Half-life must be positive")
    
    # k = ln(2) / t1/2
    return math.log(2) / half_life

def calculate_half_life_from_rate_constant(rate_constant):
    """
    Calculate half-life from rate constant (k) for first-order reactions.
    
    Args:
        rate_constant (float): Rate constant value
    
    Returns:
        float: Half-life
    """
    if rate_constant <= 0:
        raise ValueError("Rate constant must be positive")
    
    # t1/2 = ln(2) / k
    return math.log(2) / rate_constant

def calculate_concentration_after_time(initial_concentration, rate_constant, time):
    """
    Calculate concentration after a specified time for first-order reactions.
    
    Args:
        initial_concentration (float): Initial concentration
        rate_constant (float): Rate constant k
        time (float): Time elapsed
    
    Returns:
        float: Final concentration
    """
    if initial_concentration < 0:
        raise ValueError("Initial concentration cannot be negative")
    if rate_constant <= 0:
        raise ValueError("Rate constant must be positive")
    if time < 0:
        raise ValueError("Time cannot be negative")
    
    # [A]t = [A]0 * e^(-kt)
    return initial_concentration * math.exp(-rate_constant * time)

def calculate_fraction_remaining(rate_constant, time):
    """
    Calculate the fraction of initial concentration remaining after time.
    
    Args:
        rate_constant (float): Rate constant k
        time (float): Time elapsed
    
    Returns:
        float: Fraction remaining (between 0 and 1)
    """
    if rate_constant <= 0:
        raise ValueError("Rate constant must be positive")
    if time < 0:
        raise ValueError("Time cannot be negative")
    
    # Fraction = e^(-kt)
    return math.exp(-rate_constant * time)

def as_simplified_fraction(decimal, max_denominator=100):
    """
    Convert a decimal to a simplified fraction if possible.
    
    Args:
        decimal (float): Decimal value to convert
        max_denominator (int): Maximum denominator to consider
    
    Returns:
        tuple or None: (numerator, denominator) if a "nice" fraction exists, None otherwise
    """
    if decimal <= 0 or decimal >= 1:
        return None
    
    try:
        frac = Fraction(decimal).limit_denominator(max_denominator)
        return (frac.numerator, frac.denominator)
    except (ValueError, ZeroDivisionError):
        return None

# ===== ARRHENIUS EQUATION FUNCTIONS =====

def celsius_to_kelvin(celsius):
    """Convert Celsius to Kelvin."""
    return celsius + 273.15

def kelvin_to_celsius(kelvin):
    """Convert Kelvin to Celsius."""
    return kelvin - 273.15

def calculate_rate_constant_arrhenius(activation_energy, temperature_k, pre_exponential_factor):
    """
    Calculate rate constant using Arrhenius equation.
    
    Args:
        activation_energy (float): Activation energy in J/mol
        temperature_k (float): Temperature in Kelvin
        pre_exponential_factor (float): Pre-exponential factor A
    
    Returns:
        float: Rate constant k
    """
    R = 8.314  # Gas constant in J/(mol·K)
    
    if activation_energy < 0:
        raise ValueError("Activation energy cannot be negative")
    if temperature_k <= 0:
        raise ValueError("Temperature must be positive (in Kelvin)")
    if pre_exponential_factor <= 0:
        raise ValueError("Pre-exponential factor must be positive")
    
    # k = A * e^(-Ea/RT)
    return pre_exponential_factor * math.exp(-activation_energy / (R * temperature_k))

def calculate_temperature_for_rate_constant(activation_energy, k1, T1_k, k2):
    """
    Calculate temperature needed to achieve a specific rate constant.
    Uses the two-point form of Arrhenius equation.
    
    Args:
        activation_energy (float): Activation energy in J/mol
        k1 (float): Initial rate constant
        T1_k (float): Initial temperature in Kelvin
        k2 (float): Target rate constant
    
    Returns:
        float: Temperature in Kelvin needed to achieve k2
    """
    R = 8.314  # Gas constant in J/(mol·K)
    
    if activation_energy <= 0:
        raise ValueError("Activation energy must be positive")
    if k1 <= 0 or k2 <= 0:
        raise ValueError("Rate constants must be positive")
    if T1_k <= 0:
        raise ValueError("Temperature must be positive (in Kelvin)")
    
    # ln(k2/k1) = -(Ea/R) * (1/T2 - 1/T1)
    # Solving for T2: 1/T2 = 1/T1 - (R/Ea) * ln(k2/k1)
    ln_ratio = math.log(k2 / k1)
    inverse_T2 = (1 / T1_k) - (R / activation_energy) * ln_ratio
    
    if inverse_T2 <= 0:
        raise ValueError("Calculated temperature is not physically meaningful")
    
    return 1 / inverse_T2

def calculate_activation_energy(k1, T1_k, k2, T2_k):
    """
    Calculate activation energy from two rate constants at different temperatures.
    
    Args:
        k1 (float): Rate constant at temperature T1
        T1_k (float): First temperature in Kelvin
        k2 (float): Rate constant at temperature T2
        T2_k (float): Second temperature in Kelvin
    
    Returns:
        float: Activation energy in J/mol
    """
    R = 8.314  # Gas constant in J/(mol·K)
    
    if k1 <= 0 or k2 <= 0:
        raise ValueError("Rate constants must be positive")
    if T1_k <= 0 or T2_k <= 0:
        raise ValueError("Temperatures must be positive (in Kelvin)")
    if T1_k == T2_k:
        raise ValueError("Temperatures must be different")
    
    # Ea = -R * ln(k2/k1) / (1/T2 - 1/T1)
    ln_ratio = math.log(k2 / k1)
    temp_term = (1 / T2_k) - (1 / T1_k)
    
    if temp_term == 0:
        raise ValueError("Temperature difference is too small")
    
    return -R * ln_ratio / temp_term

def solve_arrhenius_problem(activation_energy=None, k1=None, T1_celsius=None, T1_kelvin=None,
                           k2=None, T2_celsius=None, T2_kelvin=None, 
                           pre_exponential_factor=None):
    """
    Solve Arrhenius equation problems given various parameters.
    
    Args:
        activation_energy (float, optional): Activation energy in J/mol
        k1 (float, optional): First rate constant
        T1_celsius (float, optional): First temperature in Celsius
        T1_kelvin (float, optional): First temperature in Kelvin
        k2 (float, optional): Second rate constant
        T2_celsius (float, optional): Second temperature in Celsius
        T2_kelvin (float, optional): Second temperature in Kelvin
        pre_exponential_factor (float, optional): Pre-exponential factor A
    
    Returns:
        dict: Solution with all calculated parameters and steps
    """
    # Convert temperatures to Kelvin if given in Celsius
    if T1_celsius is not None and T1_kelvin is None:
        T1_kelvin = celsius_to_kelvin(T1_celsius)
    elif T1_kelvin is not None and T1_celsius is None:
        T1_celsius = kelvin_to_celsius(T1_kelvin)
    
    if T2_celsius is not None and T2_kelvin is None:
        T2_kelvin = celsius_to_kelvin(T2_celsius)
    elif T2_kelvin is not None and T2_celsius is None:
        T2_celsius = kelvin_to_celsius(T2_kelvin)
    
    # Initialize result dictionary
    result = {
        "activation_energy": activation_energy,
        "k1": k1,
        "T1_celsius": T1_celsius,
        "T1_kelvin": T1_kelvin,
        "k2": k2,
        "T2_celsius": T2_celsius,
        "T2_kelvin": T2_kelvin,
        "pre_exponential_factor": pre_exponential_factor,
        "steps": []
    }
    
    steps = result["steps"]
    steps.append("Arrhenius equation forms:")
    steps.append("1. k = A × e^(-Ea/RT)")
    steps.append("2. ln(k2/k1) = -(Ea/R) × (1/T2 - 1/T1)")
    steps.append("3. R = 8.314 J/(mol·K)")
    steps.append("")
    
    # Case 1: Calculate T2 from Ea, k1, T1, k2
    if (activation_energy is not None and k1 is not None and 
        T1_kelvin is not None and k2 is not None and T2_kelvin is None):
        
        T2_kelvin = calculate_temperature_for_rate_constant(activation_energy, k1, T1_kelvin, k2)
        T2_celsius = kelvin_to_celsius(T2_kelvin)
        
        result["T2_kelvin"] = T2_kelvin
        result["T2_celsius"] = T2_celsius
        
        steps.append("Calculating T2 from given Ea, k1, T1, and k2:")
        steps.append(f"Given: Ea = {activation_energy} J/mol, k1 = {k1}, T1 = {T1_kelvin} K, k2 = {k2}")
        steps.append("Using: ln(k2/k1) = -(Ea/R) × (1/T2 - 1/T1)")
        steps.append(f"ln({k2}/{k1}) = -({activation_energy}/8.314) × (1/T2 - 1/{T1_kelvin})")
        
        ln_ratio = math.log(k2/k1)
        steps.append(f"ln({k2/k1:.6f}) = {ln_ratio:.6f}")
        steps.append(f"Solving for 1/T2:")
        steps.append(f"1/T2 = 1/{T1_kelvin} - (8.314/{activation_energy}) × {ln_ratio:.6f}")
        steps.append(f"1/T2 = {1/T1_kelvin:.8f} - {8.314/activation_energy * ln_ratio:.8f}")
        steps.append(f"1/T2 = {1/T2_kelvin:.8f}")
        steps.append(f"T2 = {T2_kelvin:.2f} K = {T2_celsius:.2f} °C")
    
    # Case 2: Calculate Ea from k1, T1, k2, T2
    elif (k1 is not None and T1_kelvin is not None and 
          k2 is not None and T2_kelvin is not None and activation_energy is None):
        
        activation_energy = calculate_activation_energy(k1, T1_kelvin, k2, T2_kelvin)
        result["activation_energy"] = activation_energy
        
        steps.append("Calculating activation energy from given rate constants and temperatures:")
        steps.append(f"Given: k1 = {k1}, T1 = {T1_kelvin} K, k2 = {k2}, T2 = {T2_kelvin} K")
        steps.append("Using: Ea = -R × ln(k2/k1) / (1/T2 - 1/T1)")
        
        ln_ratio = math.log(k2/k1)
        temp_term = (1/T2_kelvin) - (1/T1_kelvin)
        steps.append(f"ln(k2/k1) = ln({k2}/{k1}) = {ln_ratio:.6f}")
        steps.append(f"1/T2 - 1/T1 = 1/{T2_kelvin} - 1/{T1_kelvin} = {temp_term:.8f}")
        steps.append(f"Ea = -8.314 × {ln_ratio:.6f} / {temp_term:.8f}")
        steps.append(f"Ea = {activation_energy:.0f} J/mol")
    
    # Case 3: Calculate k2 from Ea, k1, T1, T2
    elif (activation_energy is not None and k1 is not None and 
          T1_kelvin is not None and T2_kelvin is not None and k2 is None):
        
        # k2 = k1 × e^[-(Ea/R) × (1/T2 - 1/T1)]
        R = 8.314
        temp_term = (1/T2_kelvin) - (1/T1_kelvin)
        exponent = -(activation_energy/R) * temp_term
        k2 = k1 * math.exp(exponent)
        
        result["k2"] = k2
        
        steps.append("Calculating k2 from given Ea, k1, T1, and T2:")
        steps.append(f"Given: Ea = {activation_energy} J/mol, k1 = {k1}, T1 = {T1_kelvin} K, T2 = {T2_kelvin} K")
        steps.append("Using: k2 = k1 × e^[-(Ea/R) × (1/T2 - 1/T1)]")
        steps.append(f"1/T2 - 1/T1 = 1/{T2_kelvin} - 1/{T1_kelvin} = {temp_term:.8f}")
        steps.append(f"Exponent = -({activation_energy}/8.314) × {temp_term:.8f} = {exponent:.6f}")
        steps.append(f"k2 = {k1} × e^({exponent:.6f}) = {k1} × {math.exp(exponent):.6f} = {k2:.6f}")
    
    return result

# ===== COMBINED SOLVER FUNCTION =====

def solve_first_order_kinetics(half_life=None, rate_constant=None, time=None, 
                              fraction_remaining=None, initial_concentration=None, 
                              final_concentration=None):
    """
    Solve first-order kinetics problems given any two parameters.
    
    Args:
        half_life (float, optional): Half-life
        rate_constant (float, optional): Rate constant k
        time (float, optional): Time elapsed
        fraction_remaining (float, optional): Fraction of initial concentration remaining
        initial_concentration (float, optional): Initial concentration
        final_concentration (float, optional): Final concentration
    
    Returns:
        dict: Solution with all calculated parameters
    """
    # Count how many parameters are provided
    params = [half_life, rate_constant, time, fraction_remaining, 
              initial_concentration, final_concentration]
    provided = sum(p is not None for p in params)
    
    if provided < 2:
        raise ValueError("At least two parameters must be provided to solve the problem")
    
    # Solution steps
    steps = ["First-order kinetics equations:"]
    steps.append("1. k = ln(2) / t₁/₂")
    steps.append("2. t₁/₂ = ln(2) / k")
    steps.append("3. [A]t = [A]₀ × e^(-kt)")
    steps.append("4. Fraction remaining = [A]t / [A]₀ = e^(-kt)")
    steps.append("5. t = ln([A]₀/[A]t) / k = -ln(fraction) / k")
    
    # Create a result dictionary
    result = {
        "half_life": half_life,
        "rate_constant": rate_constant,
        "time": time,
        "fraction_remaining": fraction_remaining,
        "initial_concentration": initial_concentration,
        "final_concentration": final_concentration,
        "fraction_as_tuple": None,
        "steps": steps
    }
    
    # Step 1: If we have half-life but not rate constant, calculate rate constant
    if half_life is not None and rate_constant is None:
        rate_constant = calculate_rate_constant_from_half_life(half_life)
        result["rate_constant"] = rate_constant
        steps.append(f"Calculating rate constant from half-life:")
        steps.append(f"k = ln(2) / t₁/₂ = ln(2) / {half_life} = {rate_constant:.6f}")
    
    # Step 2: If we have rate constant but not half-life, calculate half-life
    elif rate_constant is not None and half_life is None:
        half_life = calculate_half_life_from_rate_constant(rate_constant)
        result["half_life"] = half_life
        steps.append(f"Calculating half-life from rate constant:")
        steps.append(f"t₁/₂ = ln(2) / k = ln(2) / {rate_constant} = {half_life:.6f}")
    
    # Step 3: Calculate fraction remaining if we have rate constant and time
    if rate_constant is not None and time is not None and fraction_remaining is None:
        fraction_remaining = calculate_fraction_remaining(rate_constant, time)
        result["fraction_remaining"] = fraction_remaining
        steps.append(f"Calculating fraction remaining from rate constant and time:")
        steps.append(f"Fraction = e^(-kt) = e^(-{rate_constant} × {time}) = {fraction_remaining:.6f}")
        
        # Try to represent as simplified fraction
        fraction_as_tuple = as_simplified_fraction(fraction_remaining)
        result["fraction_as_tuple"] = fraction_as_tuple
    
    # Step 4: Calculate time if we have rate constant and fraction remaining
    elif rate_constant is not None and fraction_remaining is not None and time is None:
        if fraction_remaining <= 0 or fraction_remaining >= 1:
            raise ValueError("Fraction remaining must be between 0 and 1 exclusively")
        
        time = -math.log(fraction_remaining) / rate_constant
        result["time"] = time
        steps.append(f"Calculating time from rate constant and fraction remaining:")
        steps.append(f"t = -ln(fraction) / k = -ln({fraction_remaining}) / {rate_constant} = {time:.6f}")
    
    # Step 5: If we have initial and final concentrations but not fraction remaining
    if initial_concentration is not None and final_concentration is not None and fraction_remaining is None:
        if initial_concentration <= 0:
            raise ValueError("Initial concentration must be positive")
        if final_concentration <= 0:
            raise ValueError("Final concentration must be positive")
        if final_concentration > initial_concentration:
            raise ValueError("Final concentration cannot be greater than initial concentration")
        
        fraction_remaining = final_concentration / initial_concentration
        result["fraction_remaining"] = fraction_remaining
        steps.append(f"Calculating fraction remaining from concentrations:")
        steps.append(f"Fraction = [A]t / [A]₀ = {final_concentration} / {initial_concentration} = {fraction_remaining:.6f}")
        
        # Try to represent as simplified fraction
        fraction_as_tuple = as_simplified_fraction(fraction_remaining)
        result["fraction_as_tuple"] = fraction_as_tuple
    
    # Step 6: Calculate final concentration if we have initial concentration and fraction remaining
    elif initial_concentration is not None and fraction_remaining is not None and final_concentration is None:
        final_concentration = initial_concentration * fraction_remaining
        result["final_concentration"] = final_concentration
        steps.append(f"Calculating final concentration:")
        steps.append(f"[A]t = [A]₀ × fraction = {initial_concentration} × {fraction_remaining} = {final_concentration:.6f}")
    
    # Step 7: Calculate initial concentration if we have final concentration and fraction remaining
    elif final_concentration is not None and fraction_remaining is not None and initial_concentration is None:
        if fraction_remaining <= 0:
            raise ValueError("Fraction remaining must be positive")
        
        initial_concentration = final_concentration / fraction_remaining
        result["initial_concentration"] = initial_concentration
        steps.append(f"Calculating initial concentration:")
        steps.append(f"[A]₀ = [A]t / fraction = {final_concentration} / {fraction_remaining} = {initial_concentration:.6f}")
    
    # Calculate any missing parameters if possible
    
    # If we now have rate constant and time but not fraction remaining
    if rate_constant is not None and time is not None and result["fraction_remaining"] is None:
        fraction_remaining = calculate_fraction_remaining(rate_constant, time)
        result["fraction_remaining"] = fraction_remaining
        steps.append(f"Calculating fraction remaining from rate constant and time:")
        steps.append(f"Fraction = e^(-kt) = e^(-{rate_constant} × {time}) = {fraction_remaining:.6f}")
        
        # Try to represent as simplified fraction
        fraction_as_tuple = as_simplified_fraction(fraction_remaining)
        result["fraction_as_tuple"] = fraction_as_tuple
    
    # If we now have initial concentration, fraction, but not final concentration
    if result["initial_concentration"] is not None and result["fraction_remaining"] is not None and result["final_concentration"] is None:
        final_concentration = result["initial_concentration"] * result["fraction_remaining"]
        result["final_concentration"] = final_concentration
        steps.append(f"Calculating final concentration:")
        steps.append(f"[A]t = [A]₀ × fraction = {result['initial_concentration']} × {result['fraction_remaining']} = {final_concentration:.6f}")
    
    return result

# ===== EXAMPLE USAGE =====

if __name__ == "__main__":
    # Example 1: Solve the textbook problem
    print("=== Example 1: Textbook Arrhenius Problem ===")
    result = solve_arrhenius_problem(
        activation_energy=100000,  # 100 kJ/mol = 100,000 J/mol
        k1=50,                     # s^-1 at 25°C
        T1_celsius=25,
        k2=100                     # s^-1 (double the rate)
    )
    
    print(f"To double the rate constant from {result['k1']} to {result['k2']} s⁻¹:")
    print(f"Temperature must increase from {result['T1_celsius']:.1f}°C to {result['T2_celsius']:.1f}°C")
    print()
    for step in result['steps']:
        print(step)
    print()
    
    # Example 2: First-order kinetics problem
    print("=== Example 2: First-Order Kinetics ===")
    result2 = solve_first_order_kinetics(half_life=10, time=20)
    print(f"After {result2['time']} time units (2 half-lives):")
    print(f"Fraction remaining: {result2['fraction_remaining']:.4f}")
    print(f"Rate constant: {result2['rate_constant']:.6f}")