from chemistry_solver.molar_mass import calculate_molar_mass
from chemistry_solver.balancer import parse_chemical_equation, balance_equation, format_balanced_equation


def solve_stoichiometry_problem(equation, given_compound, given_mass, target_compound):
    reactants, products = parse_chemical_equation(equation)
    balanced_reactants, balanced_products = balance_equation(reactants, products)
    balanced_equation = format_balanced_equation(balanced_reactants, balanced_products)

    all_compounds = balanced_reactants + balanced_products
    all_dict = {f: c for c, f in all_compounds}

    if given_compound not in all_dict or target_compound not in all_dict:
        raise ValueError("Given or target compound not in balanced equation.")

    given_info = all_dict[given_compound]
    target_info = all_dict[target_compound]

    given_result = calculate_molar_mass(given_compound)
    target_result = calculate_molar_mass(target_compound)

    if not given_result['success']:
        raise ValueError(f"Error with given compound: {given_result['error']}")
    if not target_result['success']:
        raise ValueError(f"Error with target compound: {target_result['error']}")

    given_molar_mass = given_result['molar_mass']
    target_molar_mass = target_result['molar_mass']

    given_moles = given_mass / given_molar_mass
    target_moles = given_moles * (target_info / given_info)
    target_mass = target_moles * target_molar_mass

    steps = [
        f"1. Balanced equation: {balanced_equation}",
        f"2. Molar mass of {given_compound}: {given_molar_mass:.4f} g/mol",
        f"3. Molar mass of {target_compound}: {target_molar_mass:.4f} g/mol",
        f"4. Moles of {given_compound}: {given_mass} g / {given_molar_mass:.4f} g/mol = {given_moles:.6f} mol",
        f"5. Stoichiometric ratio: {target_info}:{given_info}",
        f"6. Moles of {target_compound}: {given_moles:.6f} × ({target_info}/{given_info}) = {target_moles:.6f} mol",
        f"7. Mass of {target_compound}: {target_moles:.6f} × {target_molar_mass:.4f} = {target_mass:.4f} g"
    ]

    return {
        "balanced_equation": balanced_equation,
        "target_mass": target_mass,
        "target_moles": target_moles,
        "steps": steps
    }


def solve_multireactant_problem(equation, reactant_data, target_compound):
    reactants, products = parse_chemical_equation(equation)
    br, bp = balance_equation(reactants, products)
    beq = format_balanced_equation(br, bp)

    all = {f: c for c, f in br + bp}
    if target_compound not in all:
        raise ValueError("Target not found.")

    target_coef = all[target_compound]
    target_result = calculate_molar_mass(target_compound)
    if not target_result['success']:
        raise ValueError(target_result['error'])
    target_molar_mass = target_result['molar_mass']

    yields = {}
    limiting = None
    min_mass = float('inf')

    steps = [f"Balanced: {beq}"]

    for reac, mass in reactant_data.items():
        if reac not in all:
            raise ValueError(f"{reac} not in equation")
        reac_coef = all[reac]
        reac_result = calculate_molar_mass(reac)
        if not reac_result['success']:
            raise ValueError(reac_result['error'])

        reac_molar_mass = reac_result['molar_mass']
        moles = mass / reac_molar_mass
        t_moles = moles * (target_coef / reac_coef)
        t_mass = t_moles * target_molar_mass

        yields[reac] = {
            "moles": moles,
            "mass": mass,
            "target_mass": t_mass
        }

        steps.extend([
            f"\nReactant: {reac}",
            f"  - Moles: {mass:.2f} / {reac_molar_mass:.4f} = {moles:.6f}",
            f"  - Target mass: {t_moles:.6f} × {target_molar_mass:.4f} = {t_mass:.4f} g"
        ])

        if t_mass < min_mass:
            min_mass = t_mass
            limiting = reac

    steps.append(f"\nLimiting reactant: {limiting}, yields {min_mass:.4f} g of {target_compound}")

    return {
        "balanced_equation": beq,
        "limiting_reactant": limiting,
        "target_mass": min_mass,
        "steps": steps
    }

def calculate_gas_volume(moles, temperature_c=0, pressure_atm=1.0):
    """
    Calculate the volume of gas in liters using the ideal gas law.
    
    Args:
        moles (float): Number of moles of gas
        temperature_c (float): Temperature in degrees Celsius
        pressure_atm (float): Pressure in atmospheres
    
    Returns:
        float: Volume in liters
    """
    # Gas constant R = 0.08206 L·atm/(mol·K)
    R = 0.08206
    
    # Convert temperature from Celsius to Kelvin
    temperature_k = temperature_c + 273.15
    
    # Calculate volume using PV = nRT
    volume = (moles * R * temperature_k) / pressure_atm
    
    return volume

def solve_gas_stoichiometry_problem(equation, given_compound, given_mass, target_gas, temperature_c=0, pressure_atm=1.0):
    """
    Solve a stoichiometry problem involving gas products under specific temperature and pressure.
    
    Args:
        equation (str): Chemical equation
        given_compound (str): Formula of the given compound
        given_mass (float): Mass of the given compound in grams
        target_gas (str): Formula of the target gas compound
        temperature_c (float): Temperature in degrees Celsius
        pressure_atm (float): Pressure in atmospheres
    
    Returns:
        dict: Results including steps, gas volume, etc.
    """
    # First solve the regular stoichiometry problem to get moles
    result = solve_stoichiometry_problem(equation, given_compound, given_mass, target_gas)
    
    # Calculate gas volume using ideal gas law
    target_moles = result["target_moles"]
    gas_volume = calculate_gas_volume(target_moles, temperature_c, pressure_atm)
    
    # Add gas volume calculation steps
    gas_steps = [
        f"\n8. Calculate volume of {target_gas} gas at {temperature_c}°C and {pressure_atm} atm:",
        f"   - Using ideal gas law: PV = nRT",
        f"   - V = (n × R × T) / P",
        f"   - V = ({target_moles:.6f} mol × 0.08206 L·atm/(mol·K) × {temperature_c + 273.15} K) / {pressure_atm} atm",
        f"   - V = {gas_volume:.6f} L"
    ]
    
    result["steps"].extend(gas_steps)
    result["gas_volume"] = gas_volume
    
    return result