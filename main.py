from molar_mass import calculate_molar_mass
from balancer import parse_chemical_equation, balance_equation, format_balanced_equation
from stoichiometry import (solve_stoichiometry_problem, solve_multireactant_problem, 
                        solve_gas_stoichiometry_problem)
from acid_base import identify_acid_base, analyze_compound_list
from oxidation_state import calculate_oxidation_number, display_oxidation_result
from redox_reactions import (balance_redox_reaction, identify_oxidation_changes, 
                            find_molar_ratio,
                            determine_reaction_favorability)
from functional_groups import (identify_functional_groups, check_functional_groups,
                            find_missing_functional_groups, explain_functional_groups_in_compound,
                            solve_functional_group_problem)

# Import thermodynamics module
from thermodynamics import (calculate_heat, calculate_temperature_change, calculate_final_temperature,
                           calculate_molar_heat, solve_thermal_equilibrium, solve_thermal_equilibrium_with_molar_heat,
                           handle_heat_transfer_problem, handle_heat_transfer_with_molar_heat, solve_mixture_problem,
                           calculate_boiling_point_with_pressure, calculate_pressure_with_temperature, 
                           calculate_heat_of_vaporization)

from colligative_properties import (
    solve_freezing_point_problem, calculate_boiling_point_elevation,
    calculate_osmotic_pressure, calculate_vapor_pressure_lowering,
    print_solution, FREEZING_POINT_CONSTANTS, BOILING_POINT_CONSTANTS
)

from insoluble_salts import handle_qualitative_analysis, analyze_specific_scenario

# Import the chemical_name_to_formula module
from chemical_name_to_formula import get_formula_from_name

from kinetics import (solve_first_order_kinetics, calculate_rate_constant_from_half_life,
                     calculate_half_life_from_rate_constant, calculate_concentration_after_time,
                     calculate_fraction_remaining, as_simplified_fraction)

def handle_heat_calculation():
    """
    Handler function for calculating heat using q = m × c × ΔT.
    """
    print("\n=== Heat Transfer Calculation ===")
    mass = float(input("Enter mass (g): "))
    specific_heat = float(input("Enter specific heat capacity (J/(g·K)): "))
    initial_temp = float(input("Enter initial temperature (°C): "))
    final_temp = float(input("Enter final temperature (°C): "))
    
    delta_t = final_temp - initial_temp
    heat = calculate_heat(mass, specific_heat, delta_t)
    
    print(f"\nResults:")
    print(f"Temperature change (ΔT): {delta_t:.2f} °C")
    print(f"Heat transferred (q): {heat:.2f} J")
    print(f"Heat transferred: {heat/1000:.2f} kJ")

def handle_temperature_change_calculation():
    """
    Handler function for calculating temperature change using ΔT = q / (m × c).
    """
    print("\n=== Temperature Change Calculation ===")
    heat = float(input("Enter heat energy (J): "))
    mass = float(input("Enter mass (g): "))
    specific_heat = float(input("Enter specific heat capacity (J/(g·K)): "))
    
    delta_t = calculate_temperature_change(heat, mass, specific_heat)
    
    print(f"\nResults:")
    print(f"Temperature change (ΔT): {delta_t:.2f} °C")

def handle_final_temperature_calculation():
    """
    Handler function for calculating final temperature.
    """
    print("\n=== Final Temperature Calculation ===")
    initial_temp = float(input("Enter initial temperature (°C): "))
    delta_t = float(input("Enter temperature change (°C): "))
    
    final_temp = calculate_final_temperature(initial_temp, delta_t)
    
    print(f"\nResults:")
    print(f"Final temperature: {final_temp:.2f} °C")

def handle_thermal_equilibrium():
    """
    Handler function for solving thermal equilibrium problems.
    """
    print("\n=== Thermal Equilibrium Problem ===")
    print("This will calculate the final temperature when two substances reach thermal equilibrium.")
    
    # Substance 1
    print("\nSubstance 1:")
    mass1 = float(input("Enter mass (g): "))
    specific_heat1 = float(input("Enter specific heat capacity (J/(g·K)): "))
    initial_temp1 = float(input("Enter initial temperature (°C): "))
    
    # Substance 2
    print("\nSubstance 2:")
    mass2 = float(input("Enter mass (g): "))
    specific_heat2 = float(input("Enter specific heat capacity (J/(g·K)): "))
    initial_temp2 = float(input("Enter initial temperature (°C): "))
    
    # Solve problem
    result = handle_heat_transfer_problem(
        mass1, specific_heat1, initial_temp1,
        mass2, specific_heat2, initial_temp2
    )
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal equilibrium temperature: {result['final_temp']:.2f} °C")
    print(f"Heat transferred by substance 1: {result['heat_transferred_1']:.2f} J")
    print(f"Heat transferred by substance 2: {result['heat_transferred_2']:.2f} J")

def handle_thermal_equilibrium_molar():
    """
    Handler function for solving thermal equilibrium problems using molar heat capacities.
    """
    print("\n=== Thermal Equilibrium Problem (Molar Heat Capacities) ===")
    print("This will calculate the final temperature when two substances reach thermal equilibrium.")
    
    # Substance 1
    print("\nSubstance 1:")
    name1 = input("Enter substance name or chemical formula: ")
    mass1 = float(input("Enter mass (g): "))
    
    # Try to calculate molar mass if a chemical formula is given
    try:
        molar_mass_result = calculate_molar_mass(name1)
        if molar_mass_result['success']:
            molar_mass1 = molar_mass_result['molar_mass']
            print(f"Calculated molar mass: {molar_mass1:.4f} g/mol")
        else:
            print("Could not automatically calculate molar mass.")
            molar_mass1 = float(input("Enter molar mass (g/mol): "))
    except:
        print("Could not automatically calculate molar mass.")
        molar_mass1 = float(input("Enter molar mass (g/mol): "))
    
    molar_heat_capacity1 = float(input("Enter molar heat capacity (J/(mol·K)): "))
    initial_temp1 = float(input("Enter initial temperature (°C): "))
    
    # Substance 2
    print("\nSubstance 2:")
    name2 = input("Enter substance name or chemical formula: ")
    mass2 = float(input("Enter mass (g): "))
    
    # Try to calculate molar mass if a chemical formula is given
    try:
        molar_mass_result = calculate_molar_mass(name2)
        if molar_mass_result['success']:
            molar_mass2 = molar_mass_result['molar_mass']
            print(f"Calculated molar mass: {molar_mass2:.4f} g/mol")
        else:
            print("Could not automatically calculate molar mass.")
            molar_mass2 = float(input("Enter molar mass (g/mol): "))
    except:
        print("Could not automatically calculate molar mass.")
        molar_mass2 = float(input("Enter molar mass (g/mol): "))
    
    molar_heat_capacity2 = float(input("Enter molar heat capacity (J/(mol·K)): "))
    initial_temp2 = float(input("Enter initial temperature (°C): "))
    
    # Solve problem
    result = handle_heat_transfer_with_molar_heat(
        mass1, molar_mass1, molar_heat_capacity1, initial_temp1,
        mass2, molar_mass2, molar_heat_capacity2, initial_temp2
    )
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal equilibrium temperature: {result['final_temp']:.2f} °C")
    print(f"Heat transferred by substance 1: {result['heat_transferred_1']:.2f} J")
    print(f"Heat transferred by substance 2: {result['heat_transferred_2']:.2f} J")

def handle_mixture_thermal_equilibrium():
    """
    Handler function for solving thermal equilibrium problems with multiple substances.
    """
    print("\n=== Multiple Substance Thermal Equilibrium Problem ===")
    print("This will calculate the final temperature when multiple substances reach thermal equilibrium.")
    
    substance_count = int(input("Enter number of substances: "))
    substances = []
    
    for i in range(substance_count):
        print(f"\nSubstance {i+1}:")
        name = input("Enter substance name or chemical formula: ")
        mass = float(input("Enter mass (g): "))
        
        # Try to automatically calculate molar mass
        use_molar_heat = input("Use molar heat capacity? (y/n): ").lower() == 'y'
        
        if use_molar_heat:
            # Try to calculate molar mass if a chemical formula is given
            try:
                molar_mass_result = calculate_molar_mass(name)
                if molar_mass_result['success']:
                    molar_mass = molar_mass_result['molar_mass']
                    print(f"Calculated molar mass: {molar_mass:.4f} g/mol")
                else:
                    print("Could not automatically calculate molar mass.")
                    molar_mass = float(input("Enter molar mass (g/mol): "))
            except:
                print("Could not automatically calculate molar mass.")
                molar_mass = float(input("Enter molar mass (g/mol): "))
                
            molar_heat_capacity = float(input("Enter molar heat capacity (J/(mol·K)): "))
            # Calculate specific heat
            specific_heat = molar_heat_capacity / molar_mass
            print(f"Calculated specific heat: {specific_heat:.4f} J/(g·K)")
        else:
            specific_heat = float(input("Enter specific heat capacity (J/(g·K)): "))
        
        initial_temp = float(input("Enter initial temperature (°C): "))
        
        substance = {
            'name': name,
            'mass': mass,
            'specific_heat': specific_heat,
            'initial_temp': initial_temp
        }
        
        if use_molar_heat:
            substance['molar_mass'] = molar_mass
            substance['molar_heat_capacity'] = molar_heat_capacity
        
        substances.append(substance)
    
    # Solve problem
    result = solve_mixture_problem(substances)
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal equilibrium temperature: {result['final_temp']:.2f} °C")
    print("\nHeat transferred by each substance:")
    for transfer in result["heat_transfers"]:
        print(f"  {transfer['name']}: {transfer['heat']:.2f} J (ΔT = {transfer['delta_t']:.2f} °C)")
        
def handle_predefined_scenarios():
    """
    Handler function for solving predefined qualitative analysis scenarios.
    """
    print("\n=== Predefined Qualitative Analysis Scenarios ===")
    print("1. W20_8: Unknown solution with possible Ag+, Ba2+, or Pb2+ cations")
    # Add more scenarios here as they're defined
    
    choice = input("\nEnter scenario number: ")
    
    if choice == "1":
        scenario_id = "W20_8"
    else:
        print("Invalid choice.")
        return
    
    result = analyze_specific_scenario(scenario_id)
    
    if "error" in result:
        print(f"Error: {result['error']}")
        return
    
    print("\n=== Analysis Results ===")
    print("\nAnalysis steps:")
    for step in result["steps"]:
        print(step)
    
    print("\nConclusion:")
    print(result["conclusion"])
    
    if result["identified_cations"]:
        print("\nIdentified cation(s):", ", ".join(result["identified_cations"]))
    else:
        print("\nNo cation could be identified with the given constraints.")

def handle_chemical_name_to_formula():
    """
    Handler function for converting chemical names to formulas.
    """
    print("\n=== Chemical Name to Formula Converter ===")
    name = input("Enter chemical name: ")
    result = get_formula_from_name(name)
    
    if result['success']:
        print(f"\nName: {result['name']}")
        print(f"Formula: {result['formula']}")
        print(f"IUPAC Name: {result['iupac_name']}")
        print(f"Molecular Weight: {result['weight']} g/mol")
        return result['formula']
    else:
        print(f"\nError: {result['error']}")
        return None

def handle_freezing_point_depression():
    """
    Handler function for solving molecular weight problems using freezing point depression.
    """
    print("\n=== Freezing Point Depression Calculator ===")
    
    # Ask if the user wants to enter a chemical name or formula
    input_type = input("Do you want to enter a chemical name or formula for the solvent? (name/formula): ").lower()
    
    solvent = None
    if input_type == 'name':
        print("\nConverting chemical name to formula...")
        solvent = handle_chemical_name_to_formula()
        if not solvent:
            solvent = input("\nFallback: Please enter solvent formula directly (e.g., 'H2O', 'C6H6' for benzene): ")
    else:
        solvent = input("Enter solvent formula (e.g., 'H2O', 'C6H6' for benzene): ")
    
    # Try to lookup Kf in constants
    k_f = None
    solvent_lower = solvent.lower()
    for key, value in FREEZING_POINT_CONSTANTS.items():
        if key in solvent_lower or solvent_lower in key:
            k_f = value
            print(f"Found freezing point constant for {key}: {value} °C/m")
            use_constant = input("Use this constant? (y/n): ").lower() == 'y'
            if use_constant:
                break
            else:
                k_f = None
    
    if k_f is None:
        k_f = float(input("Enter freezing point depression constant (Kf in °C/m): "))
    
    T_pure = float(input("Enter freezing point of pure solvent (°C): "))
    T_solution = float(input("Enter freezing point of solution (°C): "))
    solute_mass = float(input("Enter mass of solute (g): "))
    solvent_mass = float(input("Enter mass of solvent (g): "))
    
    ion_input = input("Does the solute ionize? (y/n): ").lower()
    ionization_factor = 1
    if ion_input == 'y':
        ionization_factor = float(input("Enter van 't Hoff factor (i): "))
    
    result = solve_freezing_point_problem(
        T_pure=T_pure,
        T_solution=T_solution,
        solvent=solvent,
        K_f=k_f,
        solute_mass=solute_mass,
        solvent_mass=solvent_mass,
        ionization_factor=ionization_factor
    )
    
    print_solution(result, "Freezing Point Depression")

def handle_boiling_point_elevation():
    """
    Handler function for solving molecular weight problems using boiling point elevation.
    """
    print("\n=== Boiling Point Elevation Calculator ===")
    
    # Ask if the user wants to enter a chemical name or formula
    input_type = input("Do you want to enter a chemical name or formula for the solvent? (name/formula): ").lower()
    
    solvent = None
    if input_type == 'name':
        print("\nConverting chemical name to formula...")
        solvent = handle_chemical_name_to_formula()
        if not solvent:
            solvent = input("\nFallback: Please enter solvent formula directly (e.g., 'H2O', 'C6H6' for benzene): ")
    else:
        solvent = input("Enter solvent formula (e.g., 'H2O', 'C6H6' for benzene): ")
    
    # Try to lookup Kb in constants
    k_b = None
    solvent_lower = solvent.lower()
    for key, value in BOILING_POINT_CONSTANTS.items():
        if key in solvent_lower or solvent_lower in key:
            k_b = value
            print(f"Found boiling point constant for {key}: {value} °C/m")
            use_constant = input("Use this constant? (y/n): ").lower() == 'y'
            if use_constant:
                break
            else:
                k_b = None
    
    if k_b is None:
        k_b = float(input("Enter boiling point elevation constant (Kb in °C/m): "))
    
    T_pure = float(input("Enter boiling point of pure solvent (°C): "))
    T_solution = float(input("Enter boiling point of solution (°C): "))
    solute_mass = float(input("Enter mass of solute (g): "))
    solvent_mass = float(input("Enter mass of solvent (g): "))
    
    ion_input = input("Does the solute ionize? (y/n): ").lower()
    ionization_factor = 1
    if ion_input == 'y':
        ionization_factor = float(input("Enter van 't Hoff factor (i): "))
    
    result = calculate_boiling_point_elevation(
        T_pure=T_pure,
        T_solution=T_solution,
        K_b=k_b,
        solute_mass=solute_mass,
        solvent=solvent,
        solvent_mass=solvent_mass,
        ionization_factor=ionization_factor
    )
    
    print_solution(result, "Boiling Point Elevation")

def handle_osmotic_pressure():
    """
    Handler function for solving molecular weight problems using osmotic pressure.
    """
    print("\n=== Osmotic Pressure Calculator ===")
    
    osmotic_pressure_atm = float(input("Enter osmotic pressure (atm): "))
    temperature_c = float(input("Enter temperature (°C): "))
    solution_volume_L = float(input("Enter solution volume (L): "))
    solute_mass = float(input("Enter mass of solute (g): "))
    
    ion_input = input("Does the solute ionize? (y/n): ").lower()
    ionization_factor = 1
    if ion_input == 'y':
        ionization_factor = float(input("Enter van 't Hoff factor (i): "))
    
    result = calculate_osmotic_pressure(
        osmotic_pressure_atm=osmotic_pressure_atm,
        temperature_c=temperature_c,
        solution_volume_L=solution_volume_L,
        solute_mass=solute_mass,
        ionization_factor=ionization_factor
    )
    
    print_solution(result, "Osmotic Pressure")

def handle_vapor_pressure_lowering():
    """
    Handler function for solving molecular weight problems using vapor pressure lowering.
    """
    print("\n=== Vapor Pressure Lowering Calculator ===")
    
    # Ask if the user wants to enter a chemical name or formula
    input_type = input("Do you want to enter a chemical name or formula for the solvent? (name/formula): ").lower()
    
    solvent = None
    if input_type == 'name':
        print("\nConverting chemical name to formula...")
        solvent = handle_chemical_name_to_formula()
        if not solvent:
            solvent = input("\nFallback: Please enter solvent formula directly (e.g., 'H2O', 'C6H6' for benzene): ")
    else:
        solvent = input("Enter solvent formula (e.g., 'H2O', 'C6H6' for benzene): ")
    
    P_pure = float(input("Enter vapor pressure of pure solvent: "))
    P_solution = float(input("Enter vapor pressure of solution: "))
    solute_mass = float(input("Enter mass of solute (g): "))
    solvent_mass = float(input("Enter mass of solvent (g): "))
    
    result = calculate_vapor_pressure_lowering(
        P_pure=P_pure,
        P_solution=P_solution,
        solute_mass=solute_mass,
        solvent=solvent,
        solvent_mass=solvent_mass
    )
    
    print_solution(result, "Vapor Pressure Lowering")
    
def handle_molar_mass():
    """
    Handler function for calculating molar mass of a compound.
    """
    # Ask if the user wants to enter a chemical name or formula
    input_type = input("Do you want to enter a chemical name or formula? (name/formula): ").lower()
    
    formula = None
    if input_type == 'name':
        print("\nConverting chemical name to formula...")
        formula = handle_chemical_name_to_formula()
        if not formula:
            formula = input("\nFallback: Please enter chemical formula directly: ")
    else:
        formula = input("Enter chemical formula: ")
    
    result = calculate_molar_mass(formula)
    if result['success']:
        print(f"Molar mass: {result['molar_mass']:.4f} g/mol")
        for e in result['composition']:
            print(f"  {e['element']}: {e['count']} × {e['atomic_mass']:.4f} u = {e['contribution']:.4f}")
    else:
        print(f"Error: {result['error']}")

def handle_balance():
    eq = input("Enter equation (e.g., H2 + O2 -> H2O): ")
    try:
        r, p = parse_chemical_equation(eq)
        br, bp = balance_equation(r, p)
        print("Balanced equation:", format_balanced_equation(br, bp))
    except Exception as e:
        print("Error:", str(e))

def handle_stoichiometry():
    eq = input("Equation: ")
    given = input("Known compound: ")
    mass = float(input("Mass in g: "))
    target = input("Target compound: ")
    try:
        result = solve_stoichiometry_problem(eq, given, mass, target)
        print("\n".join(result["steps"]))
        print(f"\nMass of {target}: {result['target_mass']:.4f} g")
    except Exception as e:
        print("Error:", str(e))

def handle_gas_stoichiometry():
    eq = input("Equation: ")
    given = input("Known compound: ")
    mass = float(input("Mass in g: "))
    target = input("Target gas compound: ")
    temp = float(input("Temperature (°C): "))
    pressure = float(input("Pressure (atm): "))
    
    try:
        result = solve_gas_stoichiometry_problem(eq, given, mass, target, temp, pressure)
        print("\n".join(result["steps"]))
        print(f"\nVolume of {target} gas: {result['gas_volume']:.4f} L")
    except Exception as e:
        print("Error:", str(e))

def handle_multireactant():
    eq = input("Equation: ")
    reactant_count = int(input("Number of reactants: "))
    reactant_data = {}
    
    for i in range(reactant_count):
        reactant = input(f"Reactant {i+1} formula: ")
        mass = float(input(f"Mass of {reactant} in g: "))
        reactant_data[reactant] = mass
    
    target = input("Target compound: ")
    
    try:
        result = solve_multireactant_problem(eq, reactant_data, target)
        print("\n".join(result["steps"]))
        print(f"\nLimiting reactant: {result['limiting_reactant']}")
        print(f"Mass of {target}: {result['target_mass']:.4f} g")
    except Exception as e:
        print("Error:", str(e))

def handle_acid_base_single():
    compound = input("Enter a chemical formula (e.g., HCl, NaOH): ")
    result = identify_acid_base(compound)
    print(f"\nCompound: {compound}")
    print(f"Classification: {result['classification']}")
    print(f"Explanation: {result['explanation']}")

def handle_acid_base_list():
    compound_input = input("Enter compounds separated by commas (e.g., HCl, H2SO4, CH3COOH): ")
    compounds = [comp.strip() for comp in compound_input.split(",")]
    analyze_compound_list(compounds)
    
def handle_oxidation_number():
    compound = input("Enter a chemical compound (e.g., CrO2Cl2, Fe2O3): ")
    element = input("Enter the element to find oxidation number for: ")
    try:
        result = calculate_oxidation_number(compound, element)
        display_oxidation_result(result)
    except Exception as e:
        print(f"Error: {str(e)}")

def handle_redox_balance():
    equation = input("Enter the redox reaction equation: ")
    environment = input("Environment (acidic/basic): ").lower()
    
    if environment not in ["acidic", "basic"]:
        environment = "acidic"  # default to acidic
    
    try:
        result = balance_redox_reaction(equation, environment)
        
        print("\n=== Balanced Redox Reaction ===")
        print(f"Balanced equation: {result['balanced_equation']}")
        print(f"Environment: {result['environment']}")
        
        print("\nElements undergoing redox:")
        for element in result['redox_elements']:
            change = element['change']
            direction = "oxidized" if change > 0 else "reduced"
            print(f"  {element['element']}: {element['reactant_oxidation']} → {element['product_oxidation']} ({direction})")
        
        print("\nOxidizing agents:", ", ".join(result['oxidizing_agents']))
        print("Reducing agents:", ", ".join(result['reducing_agents']))
    except Exception as e:
        print(f"Error: {str(e)}")

def handle_reaction_favorability():
    equation = input("Enter the redox reaction equation: ")
    
    try:
        result = determine_reaction_favorability(equation)
        
        print("\n=== Reaction Favorability Analysis ===")
        print(f"Balanced equation: {result['balanced_equation']}")
        
        if result['cell_potential'] is not None:
            print(f"\nStandard Cell Potential (E°cell): {result['cell_potential']:.2f} V")
            
            if result['favorable']:
                print("\nResult: FAVORABLE at standard conditions")
            else:
                print("\nResult: NOT FAVORABLE at standard conditions")
            
            print("\nExplanation:")
            print(result['explanation'])
        else:
            print("\nCould not determine favorability. Missing standard potentials for one or more half-reactions.")
    except Exception as e:
        print(f"Error: {str(e)}")

def handle_molar_ratio():
    equation = input("Enter the balanced equation: ")
    compound1 = input("First compound: ")
    compound2 = input("Second compound: ")
    
    try:
        ratio = find_molar_ratio(equation, compound1, compound2)
        
        if ratio is not None:
            print(f"Molar ratio of {compound1} to {compound2}: {ratio:.2f}:{1}")
            print(f"For every {ratio:.2f} moles of {compound1}, you need 1 mole of {compound2}.")
        else:
            print("Could not determine molar ratio. Make sure both compounds are in the equation.")
    except Exception as e:
        print(f"Error: {str(e)}")

# New handlers for functional groups
def handle_identify_functional_groups():
    """
    Handler function for identifying functional groups in a compound.
    """
    print("\n=== Identify Functional Groups ===")
    compound = input("Enter compound name: ")
    result = explain_functional_groups_in_compound(compound_name=compound)
    
    print(f"\nCompound: {result['compound']}")
    if result['functional_groups']:
        print("Identified functional groups:")
        for i, (group, explanation) in enumerate(zip(result['functional_groups'], result['explanations'])):
            print(f"  {i+1}. {group.capitalize()}: {explanation}")
    else:
        print("No functional groups identified.")

def handle_functional_group_problem():
    """
    Handler function for solving which functional group is NOT present in a compound.
    """
    print("\n=== Functional Group Analysis ===")
    compound = input("Enter compound name: ")
    print("Enter functional groups to check (comma-separated):")
    options_input = input("e.g., methyl, carboxyl, hydroxyl, amine, halogen: ")
    options = [opt.strip().lower() for opt in options_input.split(",")]
    
    result = solve_functional_group_problem(compound, options)
    
    print(f"\n=== Functional Group Analysis for {result['compound']} ===")
    print("Present functional groups:")
    for group in result['present_groups']:
        print(f"  - {group}")
    
    print("\nMissing functional groups:")
    for group in result['missing_groups']:
        print(f"  - {group}")
    
    if len(result['missing_groups']) == 1:
        print(f"\nAnswer: {result['missing_groups'][0]} is NOT present in {result['compound']}.")
    elif len(result['missing_groups']) > 1:
        print(f"\nAnswer: Multiple functional groups are not present in {result['compound']}.")
    else:
        print(f"\nAnswer: All specified functional groups are present in {result['compound']}.")

def handle_check_specific_functional_groups():
    """
    Handler function for checking specific functional groups in a compound.
    """
    print("\n=== Check Specific Functional Groups ===")
    compound = input("Enter compound name: ")
    print("Enter specific functional groups to check (comma-separated):")
    groups_input = input("e.g., methyl, carboxyl, hydroxyl, amine, halogen: ")
    groups_to_check = [group.strip().lower() for group in groups_input.split(",")]
    
    result = check_functional_groups(compound_name=compound, groups_to_check=groups_to_check)
    
    print(f"\nResults for {compound}:")
    for group, present in result.items():
        status = "Present" if present else "Not present"
        print(f"  {group.capitalize()}: {status}")

# New handlers for thermodynamics
def handle_heat_calculation():
    """
    Handler function for calculating heat using q = m × c × ΔT.
    """
    print("\n=== Heat Transfer Calculation ===")
    mass = float(input("Enter mass (g): "))
    specific_heat = float(input("Enter specific heat capacity (J/(g·K)): "))
    initial_temp = float(input("Enter initial temperature (°C): "))
    final_temp = float(input("Enter final temperature (°C): "))
    
    delta_t = final_temp - initial_temp
    heat = calculate_heat(mass, specific_heat, delta_t)
    
    print(f"\nResults:")
    print(f"Temperature change (ΔT): {delta_t:.2f} °C")
    print(f"Heat transferred (q): {heat:.2f} J")
    print(f"Heat transferred: {heat/1000:.2f} kJ")

def handle_temperature_change_calculation():
    """
    Handler function for calculating temperature change using ΔT = q / (m × c).
    """
    print("\n=== Temperature Change Calculation ===")
    heat = float(input("Enter heat energy (J): "))
    mass = float(input("Enter mass (g): "))
    specific_heat = float(input("Enter specific heat capacity (J/(g·K)): "))
    
    delta_t = calculate_temperature_change(heat, mass, specific_heat)
    
    print(f"\nResults:")
    print(f"Temperature change (ΔT): {delta_t:.2f} °C")

def handle_thermal_equilibrium():
    """
    Handler function for solving thermal equilibrium problems.
    """
    print("\n=== Thermal Equilibrium Problem ===")
    print("This will calculate the final temperature when two substances reach thermal equilibrium.")
    
    # Substance 1
    print("\nSubstance 1:")
    mass1 = float(input("Enter mass (g): "))
    specific_heat1 = float(input("Enter specific heat capacity (J/(g·K)): "))
    initial_temp1 = float(input("Enter initial temperature (°C): "))
    
    # Substance 2
    print("\nSubstance 2:")
    mass2 = float(input("Enter mass (g): "))
    specific_heat2 = float(input("Enter specific heat capacity (J/(g·K)): "))
    initial_temp2 = float(input("Enter initial temperature (°C): "))
    
    # Solve problem
    result = handle_heat_transfer_problem(
        mass1, specific_heat1, initial_temp1,
        mass2, specific_heat2, initial_temp2
    )
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal equilibrium temperature: {result['final_temp']:.2f} °C")
    print(f"Heat transferred by substance 1: {result['heat_transferred_1']:.2f} J")
    print(f"Heat transferred by substance 2: {result['heat_transferred_2']:.2f} J")

def handle_thermal_equilibrium_molar():
    """
    Handler function for solving thermal equilibrium problems using molar heat capacities.
    """
    print("\n=== Thermal Equilibrium Problem (Molar Heat Capacities) ===")
    print("This will calculate the final temperature when two substances reach thermal equilibrium.")
    
    # Substance 1
    print("\nSubstance 1:")
    name1 = input("Enter substance name or chemical formula: ")
    mass1 = float(input("Enter mass (g): "))
    
    # Try to calculate molar mass if a chemical formula is given
    try:
        molar_mass_result = calculate_molar_mass(name1)
        if molar_mass_result['success']:
            molar_mass1 = molar_mass_result['molar_mass']
            print(f"Calculated molar mass: {molar_mass1:.4f} g/mol")
        else:
            print("Could not automatically calculate molar mass.")
            molar_mass1 = float(input("Enter molar mass (g/mol): "))
    except:
        print("Could not automatically calculate molar mass.")
        molar_mass1 = float(input("Enter molar mass (g/mol): "))
    
    molar_heat_capacity1 = float(input("Enter molar heat capacity (J/(mol·K)): "))
    initial_temp1 = float(input("Enter initial temperature (°C): "))
    
    # Substance 2
    print("\nSubstance 2:")
    name2 = input("Enter substance name or chemical formula: ")
    mass2 = float(input("Enter mass (g): "))
    
    # Try to calculate molar mass if a chemical formula is given
    try:
        molar_mass_result = calculate_molar_mass(name2)
        if molar_mass_result['success']:
            molar_mass2 = molar_mass_result['molar_mass']
            print(f"Calculated molar mass: {molar_mass2:.4f} g/mol")
        else:
            print("Could not automatically calculate molar mass.")
            molar_mass2 = float(input("Enter molar mass (g/mol): "))
    except:
        print("Could not automatically calculate molar mass.")
        molar_mass2 = float(input("Enter molar mass (g/mol): "))
    
    molar_heat_capacity2 = float(input("Enter molar heat capacity (J/(mol·K)): "))
    initial_temp2 = float(input("Enter initial temperature (°C): "))
    
    # Solve problem
    result = handle_heat_transfer_with_molar_heat(
        mass1, molar_mass1, molar_heat_capacity1, initial_temp1,
        mass2, molar_mass2, molar_heat_capacity2, initial_temp2
    )
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal equilibrium temperature: {result['final_temp']:.2f} °C")
    print(f"Heat transferred by substance 1: {result['heat_transferred_1']:.2f} J")
    print(f"Heat transferred by substance 2: {result['heat_transferred_2']:.2f} J")

def handle_mixture_thermal_equilibrium():
    """
    Handler function for solving thermal equilibrium problems with multiple substances.
    """
    print("\n=== Multiple Substance Thermal Equilibrium Problem ===")
    print("This will calculate the final temperature when multiple substances reach thermal equilibrium.")
    
    substance_count = int(input("Enter number of substances: "))
    substances = []
    
    for i in range(substance_count):
        print(f"\nSubstance {i+1}:")
        name = input("Enter substance name or chemical formula: ")
        mass = float(input("Enter mass (g): "))
        
        # Try to automatically calculate molar mass
        use_molar_heat = input("Use molar heat capacity? (y/n): ").lower() == 'y'
        
        if use_molar_heat:
            # Try to calculate molar mass if a chemical formula is given
            try:
                molar_mass_result = calculate_molar_mass(name)
                if molar_mass_result['success']:
                    molar_mass = molar_mass_result['molar_mass']
                    print(f"Calculated molar mass: {molar_mass:.4f} g/mol")
                else:
                    print("Could not automatically calculate molar mass.")
                    molar_mass = float(input("Enter molar mass (g/mol): "))
            except:
                print("Could not automatically calculate molar mass.")
                molar_mass = float(input("Enter molar mass (g/mol): "))
                
            molar_heat_capacity = float(input("Enter molar heat capacity (J/(mol·K)): "))
            # Calculate specific heat
            specific_heat = molar_heat_capacity / molar_mass
            print(f"Calculated specific heat: {specific_heat:.4f} J/(g·K)")
        else:
            specific_heat = float(input("Enter specific heat capacity (J/(g·K)): "))
        
        initial_temp = float(input("Enter initial temperature (°C): "))
        
        substance = {
            'name': name,
            'mass': mass,
            'specific_heat': specific_heat,
            'initial_temp': initial_temp
        }
        
        if use_molar_heat:
            substance['molar_mass'] = molar_mass
            substance['molar_heat_capacity'] = molar_heat_capacity
        
        substances.append(substance)
    
    # Solve problem
    result = solve_mixture_problem(substances)
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal equilibrium temperature: {result['final_temp']:.2f} °C")
    print("\nHeat transferred by each substance:")
    for transfer in result["heat_transfers"]:
        print(f"  {transfer['name']}: {transfer['heat']:.2f} J (ΔT = {transfer['delta_t']:.2f} °C)")

def handle_pressure_boiling_point_calculation():
    """
    Handler function for calculating how pressure affects boiling point using the Clausius-Clapeyron equation.
    """
    print("\n=== Pressure Effect on Boiling Point Calculator ===")
    print("This calculates how pressure changes affect boiling point using the Clausius-Clapeyron equation.")
    
    normal_boiling_point_c = float(input("Enter the normal boiling point (°C): "))
    heat_of_vaporization = float(input("Enter the heat of vaporization (kJ/mol): "))
    initial_pressure = float(input("Enter the initial pressure (atm): "))
    final_pressure = float(input("Enter the final pressure (atm): "))
    
    result = calculate_boiling_point_with_pressure(
        normal_boiling_point_c=normal_boiling_point_c,
        heat_of_vaporization=heat_of_vaporization,
        initial_pressure=initial_pressure,
        final_pressure=final_pressure
    )
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal Result:")
    print(f"At {result['final_pressure']} atm, the boiling point of the substance is {result['new_boiling_point_c']:.2f} °C")

def handle_temperature_pressure_calculation():
    """
    Handler function for calculating vapor pressure at a given temperature using the Clausius-Clapeyron equation.
    """
    print("\n=== Vapor Pressure Calculator ===")
    print("This calculates the vapor pressure at a given temperature using the Clausius-Clapeyron equation.")
    
    normal_boiling_point_c = float(input("Enter the normal boiling point (°C): "))
    heat_of_vaporization = float(input("Enter the heat of vaporization (kJ/mol): "))
    initial_pressure = float(input("Enter the reference pressure (usually 1 atm): "))
    final_temperature_c = float(input("Enter the temperature for vapor pressure calculation (°C): "))
    
    result = calculate_pressure_with_temperature(
        normal_boiling_point_c=normal_boiling_point_c,
        heat_of_vaporization=heat_of_vaporization,
        initial_pressure=initial_pressure,
        final_temperature_c=final_temperature_c
    )
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal Result:")
    print(f"At {result['final_temperature_c']:.2f} °C, the vapor pressure of the substance is {result['final_pressure']:.6f} atm")

def handle_heat_of_vaporization_calculation():
    """
    Handler function for calculating heat of vaporization using two P-T data points.
    """
    print("\n=== Heat of Vaporization Calculator ===")
    print("This calculates the heat of vaporization using two pressure-temperature data points.")
    
    temp1_c = float(input("Enter the first temperature (°C): "))
    pressure1 = float(input("Enter the first pressure (atm): "))
    temp2_c = float(input("Enter the second temperature (°C): "))
    pressure2 = float(input("Enter the second pressure (atm): "))
    
    result = calculate_heat_of_vaporization(
        temp1_c=temp1_c,
        pressure1=pressure1,
        temp2_c=temp2_c,
        pressure2=pressure2
    )
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal Result:")
    print(f"The heat of vaporization is {result['heat_of_vaporization']:.2f} kJ/mol")


from kinetics import (solve_first_order_kinetics, calculate_rate_constant_from_half_life,
                     calculate_half_life_from_rate_constant, calculate_concentration_after_time,
                     calculate_fraction_remaining, as_simplified_fraction)

def handle_heat_of_vaporization_calculation():
    """
    Handler function for calculating heat of vaporization using two P-T data points.
    """
    print("\n=== Heat of Vaporization Calculator ===")
    print("This calculates the heat of vaporization using two pressure-temperature data points.")
    
    temp1_c = float(input("Enter the first temperature (°C): "))
    pressure1 = float(input("Enter the first pressure (atm): "))
    temp2_c = float(input("Enter the second temperature (°C): "))
    pressure2 = float(input("Enter the second pressure (atm): "))
    
    result = calculate_heat_of_vaporization(
        temp1_c=temp1_c,
        pressure1=pressure1,
        temp2_c=temp2_c,
        pressure2=pressure2
    )
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal Result:")
    print(f"The heat of vaporization is {result['heat_of_vaporization']:.2f} kJ/mol")

# Kinetics handler functions
def handle_first_order_rate_constant():
    """
    Handler function for calculating rate constant from half-life.
    """
    print("\n=== Rate Constant Calculator ===")
    print("This calculates the rate constant (k) from the half-life for a first-order reaction.")
    
    half_life = float(input("Enter the half-life: "))
    unit = input("Enter the unit of time (e.g., s, min, h): ")
    
    result = solve_first_order_kinetics(half_life=half_life)
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal Result:")
    print(f"The rate constant (k) is {result['rate_constant']:.6f} {unit}⁻¹")

def handle_first_order_half_life():
    """
    Handler function for calculating half-life from rate constant.
    """
    print("\n=== Half-Life Calculator ===")
    print("This calculates the half-life from the rate constant for a first-order reaction.")
    
    rate_constant = float(input("Enter the rate constant (k): "))
    unit = input("Enter the unit of time (e.g., s, min, h): ")
    
    result = solve_first_order_kinetics(rate_constant=rate_constant)
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal Result:")
    print(f"The half-life is {result['half_life']:.6f} {unit}")

def handle_first_order_concentration():
    """
    Handler function for calculating concentration after a given time.
    """
    print("\n=== Concentration Calculator ===")
    print("This calculates the concentration after a given time using the first-order rate law.")
    
    initial_concentration = float(input("Enter the initial concentration [A]₀: "))
    conc_unit = input("Enter the concentration unit (e.g., M, mol/L): ")
    rate_constant = float(input("Enter the rate constant (k): "))
    time = float(input("Enter the time elapsed: "))
    time_unit = input("Enter the time unit (e.g., s, min, h): ")
    
    result = solve_first_order_kinetics(
        initial_concentration=initial_concentration,
        rate_constant=rate_constant,
        time=time
    )
    
    final_concentration = result["final_concentration"]
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal Result:")
    print(f"The concentration after {time} {time_unit} is {final_concentration:.6f} {conc_unit}")
    print(f"The fraction remaining is {result['fraction_remaining']:.6f}")
    
    if result["fraction_as_tuple"]:
        num, den = result["fraction_as_tuple"]
        if den <= 100:  # Only show if denominator is reasonable
            print(f"This fraction can be expressed as {num}/{den}")

def handle_first_order_time():
    """
    Handler function for calculating time needed to reach a specified concentration.
    """
    print("\n=== Time Calculator ===")
    print("This calculates the time needed to reach a specified concentration.")
    
    initial_concentration = float(input("Enter the initial concentration [A]₀: "))
    conc_unit = input("Enter the concentration unit (e.g., M, mol/L): ")
    
    calculation_type = input("Calculate time using (1) final concentration or (2) fraction remaining? (1/2): ")
    
    if calculation_type == "1":
        final_concentration = float(input(f"Enter the final concentration [A]t in {conc_unit}: "))
        fraction_remaining = final_concentration / initial_concentration
    else:
        fraction_remaining = float(input("Enter the fraction remaining (e.g., 0.5 for half): "))
        final_concentration = initial_concentration * fraction_remaining
        
    rate_constant = float(input("Enter the rate constant (k): "))
    time_unit = input("Enter the time unit (e.g., s, min, h): ")
    
    result = solve_first_order_kinetics(
        initial_concentration=initial_concentration,
        final_concentration=final_concentration,
        rate_constant=rate_constant,
        fraction_remaining=fraction_remaining
    )
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal Result:")
    print(f"The time needed to reach {final_concentration:.6f} {conc_unit} is {result['time']:.6f} {time_unit}")
    print(f"This corresponds to {result['time'] / result['half_life']:.4f} half-lives.")

def handle_first_order_complete_analysis():
    """
    Handler function for complete first-order kinetics analysis.
    """
    print("\n=== Complete First-Order Kinetics Analysis ===")
    print("This performs a complete analysis given any two parameters.")
    
    print("\nWhat information do you have? Enter values for any two of the following:")
    print("(Leave blank if unknown)")
    
    half_life_input = input("Half-life (leave blank if unknown): ")
    half_life = float(half_life_input) if half_life_input else None
    
    rate_constant_input = input("Rate constant k (leave blank if unknown): ")
    rate_constant = float(rate_constant_input) if rate_constant_input else None
    
    time_input = input("Time elapsed (leave blank if unknown): ")
    time = float(time_input) if time_input else None
    
    fraction_input = input("Fraction remaining (e.g., 0.5 for half) (leave blank if unknown): ")
    fraction_remaining = float(fraction_input) if fraction_input else None
    
    initial_conc_input = input("Initial concentration [A]₀ (leave blank if unknown): ")
    initial_concentration = float(initial_conc_input) if initial_conc_input else None
    
    final_conc_input = input("Final concentration [A]t (leave blank if unknown): ")
    final_concentration = float(final_conc_input) if final_conc_input else None
    
    unit = input("Enter the unit of time (e.g., s, min, h): ")
    conc_unit = input("Enter the concentration unit (if applicable, e.g., M, mol/L): ") if initial_concentration or final_concentration else None
    
    try:
        result = solve_first_order_kinetics(
            half_life=half_life,
            rate_constant=rate_constant,
            initial_concentration=initial_concentration,
            final_concentration=final_concentration,
            time=time,
            fraction_remaining=fraction_remaining
        )
        
        print("\n=== Solution ===")
        for step in result["steps"]:
            print(step)
        
        print(f"\nFinal Results:")
        if result["half_life"] is not None:
            print(f"Half-life: {result['half_life']:.6f} {unit}")
        if result["rate_constant"] is not None:
            print(f"Rate constant (k): {result['rate_constant']:.6f} {unit}⁻¹")
        if result["time"] is not None:
            print(f"Time elapsed: {result['time']:.6f} {unit}")
        if result["fraction_remaining"] is not None:
            print(f"Fraction remaining: {result['fraction_remaining']:.6f}")
            if result["fraction_as_tuple"]:
                num, den = result["fraction_as_tuple"]
                if den <= 100:  # Only show if denominator is reasonable
                    print(f"This fraction can be expressed as {num}/{den}")
        if result["initial_concentration"] is not None:
            print(f"Initial concentration [A]₀: {result['initial_concentration']:.6f} {conc_unit}")
        if result["final_concentration"] is not None:
            print(f"Final concentration [A]t: {result['final_concentration']:.6f} {conc_unit}")
            
    except Exception as e:
        print(f"\nError: {e}")
        print("Please provide sufficient information to solve the problem.")
    
def main():
    actions = {
        '1': handle_molar_mass,
        '2': handle_balance,
        '3': handle_stoichiometry,
        '4': handle_multireactant,
        '5': handle_acid_base_single,
        '6': handle_acid_base_list,
        '7': handle_oxidation_number,
        '8': handle_redox_balance,
        '9': handle_molar_ratio,
        '10': handle_gas_stoichiometry,
        '11': handle_reaction_favorability,
        '12': handle_identify_functional_groups,
        '13': handle_check_specific_functional_groups,
        '14': handle_functional_group_problem,
        '15': handle_freezing_point_depression,  
        '16': handle_boiling_point_elevation,    
        '17': handle_osmotic_pressure,           
        '18': handle_vapor_pressure_lowering,
        '19': handle_chemical_name_to_formula,
        '20': handle_qualitative_analysis,
        '21': handle_predefined_scenarios,
        # Thermodynamics handlers
        '22': handle_heat_calculation,
        '23': handle_temperature_change_calculation,
        '24': handle_final_temperature_calculation,
        '25': handle_thermal_equilibrium,
        '26': handle_thermal_equilibrium_molar,
        '27': handle_mixture_thermal_equilibrium,
        # Clausius-Clapeyron handlers
        '28': handle_pressure_boiling_point_calculation,
        '29': handle_temperature_pressure_calculation,
        '30': handle_heat_of_vaporization_calculation,
        # First-order kinetics handlers
        '31': handle_first_order_rate_constant,
        '32': handle_first_order_half_life,
        '33': handle_first_order_concentration,
        '34': handle_first_order_time,
        '35': handle_first_order_complete_analysis,
        '0': exit
    }
    while True:
        print("\n===== Chemical Equation Solver =====")
        print("1. Calculate molar mass")
        print("2. Balance a chemical equation")
        print("3. Solve stoichiometry problem")
        print("4. Solve multireactant problem (limiting reactant)")
        print("5. Identify acid/base compound")
        print("6. Analyze list of compounds")
        print("7. Calculate oxidation number")
        print("8. Balance a redox reaction")
        print("9. Find molar ratio between compounds")
        print("10. Solve gas stoichiometry problem")
        print("11. Determine redox reaction favorability")
        print("12. Identify functional groups in a compound")
        print("13. Check specific functional groups")
        print("14. Find missing functional groups in a compound")
        print("---- Colligative Properties ----")
        print("15. Freezing point depression (molecular weight)")
        print("16. Boiling point elevation (molecular weight)")
        print("17. Osmotic pressure (molecular weight)")
        print("18. Vapor pressure lowering (molecular weight)")
        print("19. Convert chemical name to formula")
        print("---- Qualitative Analysis ----")
        print("20. Solve general qualitative analysis problem")
        print("21. Solve predefined qualitative analysis scenario")
        print("---- Thermodynamics ----")
        print("22. Calculate heat energy (q = m × c × ΔT)")
        print("23. Calculate temperature change (ΔT = q / (m × c))")
        print("24. Calculate final temperature (Tf = Ti + ΔT)")
        print("25. Solve thermal equilibrium problem (two substances)")
        print("26. Solve thermal equilibrium with molar heat capacities")
        print("27. Solve thermal equilibrium with multiple substances")
        print("---- Clausius-Clapeyron Equation ----")
        print("28. Calculate boiling point at different pressure")
        print("29. Calculate vapor pressure at different temperature")
        print("30. Calculate heat of vaporization from P-T data points")
        print("---- First-Order Kinetics ----")
        print("31. Calculate rate constant from half-life")
        print("32. Calculate half-life from rate constant")
        print("33. Calculate concentration after time")
        print("34. Calculate time to reach concentration")
        print("35. Complete first-order kinetics analysis")
        print("0. Exit")
        
        choice = input("Enter choice: ")
        action = actions.get(choice)
        if action:
            action()
        else:
            print("Invalid choice.")

if __name__ == "__main__":
    main()