"""
Terminal User Interface for Thermodynamics Calculations
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user
from chemistry_solver.thermodynamics import (calculate_heat, calculate_temperature_change, calculate_final_temperature,
                      calculate_molar_heat, solve_thermal_equilibrium, solve_thermal_equilibrium_with_molar_heat,
                      handle_heat_transfer_problem, handle_heat_transfer_with_molar_heat, solve_mixture_problem,
                      calculate_boiling_point_with_pressure, calculate_pressure_with_temperature, 
                      calculate_heat_of_vaporization, parse_reaction_string, solve_enthalpy_problem, 
                      solve_equilibrium_constant_problem, solve_equilibrium_from_formation_data)
# Import the name_to_formula module
try:
    from chemistry_solver.name_to_formula import get_formula_from_name
except ImportError:
    # Define a fallback function if the module is not available
    def get_formula_from_name(compound_name):
        return {'success': False, 'error': 'pubchempy module not installed'}

try:
    from chemistry_solver.molar_mass import calculate_molar_mass
except ImportError:
    # Define a fallback function if the module is not available
    def calculate_molar_mass(formula):
        return {'success': False, 'error': 'molmass module not installed'}

# Common substances with their molar masses and formulas as fallbacks
COMMON_SUBSTANCES = {
    'water': {'formula': 'H2O', 'molar_mass': 18.01528},
    'ethanol': {'formula': 'C2H5OH', 'molar_mass': 46.06844},
    'methanol': {'formula': 'CH3OH', 'molar_mass': 32.04186},
    'carbon dioxide': {'formula': 'CO2', 'molar_mass': 44.0095},
    'oxygen': {'formula': 'O2', 'molar_mass': 31.9988},
    'nitrogen': {'formula': 'N2', 'molar_mass': 28.0134},
    'hydrogen': {'formula': 'H2', 'molar_mass': 2.01588},
    'acetone': {'formula': 'C3H6O', 'molar_mass': 58.08},
    'benzene': {'formula': 'C6H6', 'molar_mass': 78.11},
    'glucose': {'formula': 'C6H12O6', 'molar_mass': 180.156},
    'salt': {'formula': 'NaCl', 'molar_mass': 58.44},
    'sodium chloride': {'formula': 'NaCl', 'molar_mass': 58.44},
    'sulfuric acid': {'formula': 'H2SO4', 'molar_mass': 98.079},
    'hydrochloric acid': {'formula': 'HCl', 'molar_mass': 36.46},
    'ammonia': {'formula': 'NH3', 'molar_mass': 17.03},
    'gold': {'formula': 'Au', 'molar_mass': 196.97},
    'silver': {'formula': 'Ag', 'molar_mass': 107.87},
    'copper': {'formula': 'Cu', 'molar_mass': 63.55},
    'iron': {'formula': 'Fe', 'molar_mass': 55.85},
    'aluminum': {'formula': 'Al', 'molar_mass': 26.98},
    'zinc': {'formula': 'Zn', 'molar_mass': 65.38},
    'lead': {'formula': 'Pb', 'molar_mass': 207.2},
    'mercury': {'formula': 'Hg', 'molar_mass': 200.59},
    'calcium carbonate': {'formula': 'CaCO3', 'molar_mass': 100.09},
    'sodium hydroxide': {'formula': 'NaOH', 'molar_mass': 40.00},
    'potassium hydroxide': {'formula': 'KOH', 'molar_mass': 56.11},
    'acetic acid': {'formula': 'CH3COOH', 'molar_mass': 60.05},
}

class ThermodynamicsUI:
    """UI class for thermodynamics calculations."""
    
    def __init__(self):
        self.title = "THERMODYNAMICS CALCULATOR"
    
    def run(self):
        """Run the thermodynamics UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-9): ").strip()  # Changed from 0-8 to 0-9
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_heat_calculation()
            elif choice == "2":
                self._handle_temperature_change_calculation()
            elif choice == "3":
                self._handle_final_temperature_calculation()
            elif choice == "4":
                self._handle_thermal_equilibrium()
            elif choice == "5":
                self._handle_thermal_equilibrium_molar()
            elif choice == "6":
                self._handle_mixture_thermal_equilibrium()
            elif choice == "7":
                self._handle_clausius_clapeyron_calculations()
            elif choice == "8":
                self._handle_hess_law_calculations()
            elif choice == "9":  # NEW OPTION
                self._handle_equilibrium_constant_calculations()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the thermodynamics module menu."""
        menu = """
        [1] Calculate heat transfer
        [2] Calculate temperature change
        [3] Calculate final temperature
        [4] Solve thermal equilibrium problem
        [5] Solve thermal equilibrium with molar heat capacities
        [6] Solve multi-substance thermal equilibrium problem
        [7] Clausius-Clapeyron equation calculations
        [8] Hess's Law calculations
        [9] Equilibrium constant calculations
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_heat_calculation(self):
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
        
        display_results_header()
        print(f"Temperature change (ΔT): {delta_t:.2f} °C")
        print(f"Heat transferred (q): {heat:.2f} J")
        print(f"Heat transferred: {heat/1000:.2f} kJ")

    def _handle_temperature_change_calculation(self):
        """
        Handler function for calculating temperature change using ΔT = q / (m × c).
        """
        print("\n=== Temperature Change Calculation ===")
        heat = float(input("Enter heat energy (J): "))
        mass = float(input("Enter mass (g): "))
        specific_heat = float(input("Enter specific heat capacity (J/(g·K)): "))
        
        delta_t = calculate_temperature_change(heat, mass, specific_heat)
        
        display_results_header()
        print(f"Temperature change (ΔT): {delta_t:.2f} °C")

    def _handle_final_temperature_calculation(self):
        """
        Handler function for calculating final temperature.
        """
        print("\n=== Final Temperature Calculation ===")
        initial_temp = float(input("Enter initial temperature (°C): "))
        delta_t = float(input("Enter temperature change (°C): "))
        
        final_temp = calculate_final_temperature(initial_temp, delta_t)
        
        display_results_header()
        print(f"Final temperature: {final_temp:.2f} °C")

    def _handle_thermal_equilibrium(self):
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
        
        display_results_header()
        for step in result["steps"]:
            print(step)
        
        print(f"\nFinal equilibrium temperature: {result['final_temp']:.2f} °C")
        print(f"Heat transferred by substance 1: {result['heat_transferred_1']:.2f} J")
        print(f"Heat transferred by substance 2: {result['heat_transferred_2']:.2f} J")

    def _get_molar_mass_from_input(self, name_or_formula):
        """
        Helper function to determine molar mass from a substance name or formula.
        Tries multiple methods:
        1. Direct calculation if it's a formula
        2. Looking up in common substances dictionary
        3. Getting formula from name via PubChem and then calculate
        4. Asking the user if all else fails
        """
        name_or_formula = name_or_formula.strip().lower()
        
        # Check if it's in our common substances dictionary
        if name_or_formula in COMMON_SUBSTANCES:
            formula = COMMON_SUBSTANCES[name_or_formula]['formula']
            molar_mass = COMMON_SUBSTANCES[name_or_formula]['molar_mass']
            print(f"Found in common substances: {formula}, molar mass = {molar_mass:.4f} g/mol")
            return formula, molar_mass
        
        # Try direct calculation if it might be a formula
        try:
            molar_mass_result = calculate_molar_mass(name_or_formula)
            if molar_mass_result['success']:
                formula = name_or_formula
                molar_mass = molar_mass_result['molar_mass']
                print(f"Calculated molar mass from formula: {molar_mass:.4f} g/mol")
                return formula, molar_mass
        except Exception:
            pass
        
        # Try to get formula from name using PubChem
        try:
            formula_result = get_formula_from_name(name_or_formula)
            if formula_result['success']:
                formula = formula_result['formula']
                molar_mass = formula_result['weight']
                print(f"Found formula via PubChem: {formula}, molar mass = {molar_mass:.4f} g/mol")
                return formula, molar_mass
        except Exception:
            pass
        
        # If we couldn't determine it automatically, ask the user
        print("Could not automatically determine molar mass.")
        formula = input(f"Enter chemical formula for {name_or_formula}: ")
        
        # Try again with the provided formula
        try:
            molar_mass_result = calculate_molar_mass(formula)
            if molar_mass_result['success']:
                molar_mass = molar_mass_result['molar_mass']
                print(f"Calculated molar mass: {molar_mass:.4f} g/mol")
                return formula, molar_mass
        except Exception:
            pass
        
        # As a last resort, ask for the molar mass directly
        molar_mass = float(input("Enter molar mass (g/mol): "))
        return formula, molar_mass

    def _handle_thermal_equilibrium_molar(self):
        """
        Handler function for solving thermal equilibrium problems using molar heat capacities.
        """
        print("\n=== Thermal Equilibrium Problem (Molar Heat Capacities) ===")
        print("This will calculate the final temperature when two substances reach thermal equilibrium.")
        
        # Substance 1
        print("\nSubstance 1:")
        name1 = input("Enter substance name or chemical formula: ")
        mass1 = float(input("Enter mass (g): "))
        
        # Try to calculate molar mass if a chemical name or formula is given
        formula1, molar_mass1 = self._get_molar_mass_from_input(name1)
        
        molar_heat_capacity1 = float(input("Enter molar heat capacity (J/(mol·K)): "))
        initial_temp1 = float(input("Enter initial temperature (°C): "))
        
        # Substance 2
        print("\nSubstance 2:")
        name2 = input("Enter substance name or chemical formula: ")
        mass2 = float(input("Enter mass (g): "))
        
        # Try to calculate molar mass if a chemical name or formula is given
        formula2, molar_mass2 = self._get_molar_mass_from_input(name2)
        
        molar_heat_capacity2 = float(input("Enter molar heat capacity (J/(mol·K)): "))
        initial_temp2 = float(input("Enter initial temperature (°C): "))
        
        # Solve problem
        result = handle_heat_transfer_with_molar_heat(
            mass1, molar_mass1, molar_heat_capacity1, initial_temp1,
            mass2, molar_mass2, molar_heat_capacity2, initial_temp2
        )
        
        display_results_header()
        for step in result["steps"]:
            print(step)
        
        print(f"\nFinal equilibrium temperature: {result['final_temp']:.2f} °C")
        print(f"Heat transferred by substance 1: {result['heat_transferred_1']:.2f} J")
        print(f"Heat transferred by substance 2: {result['heat_transferred_2']:.2f} J")

    def _handle_mixture_thermal_equilibrium(self):
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
                formula, molar_mass = self._get_molar_mass_from_input(name)
                    
                molar_heat_capacity = float(input("Enter molar heat capacity (J/(mol·K)): "))
                # Calculate specific heat
                specific_heat = molar_heat_capacity / molar_mass
                print(f"Calculated specific heat: {specific_heat:.4f} J/(g·K)")
            else:
                specific_heat = float(input("Enter specific heat capacity (J/(g·K)): "))
                formula = name
                molar_mass = None
                molar_heat_capacity = None
            
            initial_temp = float(input("Enter initial temperature (°C): "))
            
            substance = {
                'name': name,
                'formula': formula,
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
        
        display_results_header()
        for step in result["steps"]:
            print(step)
        
        print(f"\nFinal equilibrium temperature: {result['final_temp']:.2f} °C")
        print("\nHeat transferred by each substance:")
        for transfer in result["heat_transfers"]:
            print(f"  {transfer['name']}: {transfer['heat']:.2f} J (ΔT = {transfer['delta_t']:.2f} °C)")

    def _handle_clausius_clapeyron_calculations(self):
        """
        Handler function for Clausius-Clapeyron equation calculations.
        """
        print("\n=== Clausius-Clapeyron Equation Calculations ===")
        print("Select the type of calculation:")
        print("[1] Calculate boiling point at a given pressure")
        print("[2] Calculate vapor pressure at a given temperature")
        print("[3] Calculate heat of vaporization from two data points")
        
        choice = input("\nEnter choice (1-3): ").strip()
        
        if choice == "1":
            print("\n--- Calculate Boiling Point at a Given Pressure ---")
            normal_bp = float(input("Enter normal boiling point (°C): "))
            heat_vap = float(input("Enter heat of vaporization (kJ/mol): "))
            init_pressure = float(input("Enter initial pressure (atm, default=1): ") or "1")
            final_pressure = float(input("Enter final pressure (atm): "))
            
            result = calculate_boiling_point_with_pressure(
                normal_bp, heat_vap, init_pressure, final_pressure
            )
            
            display_results_header()
            for step in result["steps"]:
                print(step)
            
            print(f"\nThe boiling point at {final_pressure} atm is {result['new_boiling_point_c']:.2f} °C")
            
        elif choice == "2":
            print("\n--- Calculate Vapor Pressure at a Given Temperature ---")
            normal_bp = float(input("Enter normal boiling point (°C): "))
            heat_vap = float(input("Enter heat of vaporization (kJ/mol): "))
            init_pressure = float(input("Enter initial pressure (atm, default=1): ") or "1")
            final_temp = float(input("Enter temperature (°C): "))
            
            result = calculate_pressure_with_temperature(
                normal_bp, heat_vap, init_pressure, final_temp
            )
            
            display_results_header()
            for step in result["steps"]:
                print(step)
            
            print(f"\nThe vapor pressure at {final_temp} °C is {result['final_pressure']:.6f} atm")
            
        elif choice == "3":
            print("\n--- Calculate Heat of Vaporization from Two Data Points ---")
            temp1 = float(input("Enter first temperature (°C): "))
            pressure1 = float(input("Enter pressure at first temperature (atm): "))
            temp2 = float(input("Enter second temperature (°C): "))
            pressure2 = float(input("Enter pressure at second temperature (atm): "))
            
            result = calculate_heat_of_vaporization(
                temp1, pressure1, temp2, pressure2
            )
            
            display_results_header()
            for step in result["steps"]:
                print(step)
            
            print(f"\nThe heat of vaporization is {result['heat_of_vaporization']:.2f} kJ/mol")
            
        else:
            print("Invalid choice.")
            
    def _handle_hess_law_calculations(self):
        """
        Handler function for Hess's Law calculations.
        """
        print("\n=== Hess's Law Calculations ===")
        print("This will calculate enthalpy changes using Hess's Law.")
        print("\nYou have two options:")
        print("[1] Enter individual reactions")
        print("[2] Enter a complete problem")
        
        choice = input("\nEnter choice (1-2): ").strip()
        
        if choice == "1":
            # Manual entry of individual reactions
            print("\n--- Enter Known Reactions ---")
            print("For each reaction, use the format: 'reactants → products ΔH° = value'")
            print("Example: 'C (graphite) + O2 (g) → CO2 (g) ΔH° = -393.5 kJ/mol'")
            
            known_reactions = []
            num_known = int(input("\nHow many known reactions do you want to enter? "))
            
            for i in range(num_known):
                reaction_str = input(f"Reaction {i+1}: ")
                try:
                    reaction = parse_reaction_string(reaction_str)
                    known_reactions.append(reaction)
                    print(f"Parsed: {reaction}")
                except ValueError as e:
                    print(f"Error parsing reaction: {e}")
                    print("Please try again.")
                    reaction_str = input(f"Reaction {i+1}: ")
                    reaction = parse_reaction_string(reaction_str)
                    known_reactions.append(reaction)
            
            print("\n--- Enter Target Reaction ---")
            print("Enter the reaction for which you want to calculate the enthalpy change:")
            target_str = input("Target reaction: ")
            target_reaction = parse_reaction_string(target_str)
            
            # Calculate the result
            result = solve_combustion_enthalpy(known_reactions, target_reaction)
            
            display_results_header()
            for step in result["steps"]:
                print(step)
            
        elif choice == "2":
            # Complete problem entry
            print("\n--- Enter Complete Problem ---")
            print("Enter the complete problem with all reactions (one per line).")
            print("The reaction without an enthalpy value will be considered the target.")
            print("Example:")
            print("C (graphite) + O2 (g) → CO2 (g) ΔH° = -393.5 kJ/mol")
            print("H2 (g) + 1/2 O2 (g) → H2O (l) ΔH° = -285.8 kJ/mol")
            print("CH3OH (l) → C (graphite) + 2 H2 (g) + 1/2 O2 (g) ΔH° = 239.0 kJ/mol")
            print("CH3OH (l) + 3/2 O2 (g) → CO2 (g) + 2 H2O (l)")
            
            print("\nEnter your problem (end with a blank line):")
            problem_lines = []
            while True:
                line = input()
                if not line:
                    break
                problem_lines.append(line)
            
            problem_text = "\n".join(problem_lines)
            
            # Solve the problem
            result = solve_enthalpy_problem(problem_text)
            
            display_results_header()
            for step in result["steps"]:
                print(step)
            
        else:
            print("Invalid choice.")
    
    def _handle_equilibrium_constant_calculations(self):
        """
        Handler function for equilibrium constant calculations.
        """
        print("\n=== Equilibrium Constant Calculations ===")
        print("Select the type of calculation:")
        print("[1] Calculate equilibrium constant from ΔH° and ΔS° (or ΔG°)")
        print("[2] Calculate equilibrium constant from formation data")
        
        choice = input("\nEnter choice (1-2): ").strip()
        
        if choice == "1":
            print("\n--- Calculate Equilibrium Constant ---")
            print("You can provide either:")
            print("- ΔG° directly, OR")
            print("- Both ΔH° and ΔS° to calculate ΔG°")
            
            # Get reaction string (optional)
            reaction = input("\nEnter balanced chemical equation (optional): ").strip() or None
            
            # Get temperature
            temp_choice = input("Enter temperature in (C)elsius or (K)elvin? ").lower()
            if temp_choice.startswith('c'):
                temp_c = float(input("Enter temperature (°C): "))
                temp_k = None
            else:
                temp_k = float(input("Enter temperature (K): "))
                temp_c = None
            
            # Get thermodynamic data
            data_choice = input("\nDo you have ΔG° directly? (y/n): ").lower()
            
            if data_choice == 'y':
                delta_g = float(input("Enter ΔG° (kJ/mol): "))
                delta_h = None
                delta_s = None
            else:
                print("Enter ΔH° and ΔS° to calculate ΔG°:")
                delta_h = float(input("Enter ΔH° (kJ/mol): "))
                delta_s = float(input("Enter ΔS° (J/(mol·K)): "))
                delta_g = None
            
            # Solve the problem
            result = solve_equilibrium_constant_problem(
                delta_h_standard=delta_h,
                delta_s_standard=delta_s,
                delta_g_standard=delta_g,
                temperature_c=temp_c,
                temperature_k=temp_k,
                reaction_string=reaction
            )
            
            display_results_header()
            for step in result["steps"]:
                print(step)
            
            print(f"\n=== FINAL RESULTS ===")
            print(f"Temperature: {result['temperature_c']:.2f} °C ({result['temperature_k']:.2f} K)")
            if result['delta_g_standard'] is not None:
                print(f"ΔG°: {result['delta_g_standard']:.3f} kJ/mol")
            print(f"Equilibrium constant (K): {result['equilibrium_constant']:.2e}")
            
            # Interpret the result
            if result['equilibrium_constant'] > 1:
                print("Since K > 1, the reaction favors products at equilibrium.")
            elif result['equilibrium_constant'] < 1:
                print("Since K < 1, the reaction favors reactants at equilibrium.")
            else:
                print("Since K ≈ 1, reactants and products are present in similar amounts at equilibrium.")
                
        elif choice == "2":
            print("\n--- Calculate Equilibrium Constant from Formation Data ---")
            print("This will calculate the equilibrium constant using standard formation enthalpies and entropies.")
            
            # Get reaction string
            reaction = input("\nEnter balanced chemical equation: ").strip()
            
            # Get temperature
            temp_c = float(input("Enter temperature (°C, default=25): ") or "25")
            
            # Parse reactants and products from user input
            print("\n--- Enter Reactants ---")
            reactants = {}
            num_reactants = int(input("How many different reactants? "))
            
            for i in range(num_reactants):
                substance = input(f"Reactant {i+1} name/formula: ").strip()
                coeff = float(input(f"Coefficient for {substance}: "))
                reactants[substance] = coeff
            
            print("\n--- Enter Products ---")
            products = {}
            num_products = int(input("How many different products? "))
            
            for i in range(num_products):
                substance = input(f"Product {i+1} name/formula: ").strip()
                coeff = float(input(f"Coefficient for {substance}: "))
                products[substance] = coeff
            
            reaction_coefficients = {
                'reactants': reactants,
                'products': products
            }
            
            # Get all unique substances
            all_substances = set(reactants.keys()) | set(products.keys())
            
            print(f"\n--- Enter Formation Enthalpies (ΔHf°) ---")
            print("Note: Elements in their standard state have ΔHf° = 0")
            formation_enthalpies = {}
            
            for substance in all_substances:
                default_hint = " (0 for elements)" if substance in ['H2', 'N2', 'O2', 'F2', 'Cl2', 'Br2', 'I2'] else ""
                value = input(f"ΔHf° for {substance} (kJ/mol){default_hint}: ")
                formation_enthalpies[substance] = float(value) if value else 0.0
            
            print(f"\n--- Enter Standard Entropies (S°) ---")
            standard_entropies = {}
            
            for substance in all_substances:
                value = float(input(f"S° for {substance} (J/(mol·K)): "))
                standard_entropies[substance] = value
            
            # Solve the problem
            result = solve_equilibrium_from_formation_data(
                reaction_coefficients=reaction_coefficients,
                formation_enthalpies=formation_enthalpies,
                standard_entropies=standard_entropies,
                temperature_c=temp_c,
                reaction_string=reaction
            )
            
            display_results_header()
            for step in result["steps"]:
                print(step)
            
            print(f"\n=== FINAL RESULTS ===")
            print(f"Reaction: {result['reaction']}")
            print(f"Temperature: {result['temperature_c']:.2f} °C ({result['temperature_k']:.2f} K)")
            print(f"ΔH°rxn: {result['delta_h_standard']:.1f} kJ/mol")
            print(f"ΔS°rxn: {result['delta_s_standard']:.1f} J/(mol·K)")
            print(f"ΔG°rxn: {result['delta_g_standard']:.3f} kJ/mol")
            print(f"Equilibrium constant (K): {result['equilibrium_constant']:.2e}")
            
            # Interpret the result
            if result['equilibrium_constant'] > 1:
                print("Since K > 1, the reaction favors products at equilibrium.")
            elif result['equilibrium_constant'] < 1:
                print("Since K < 1, the reaction favors reactants at equilibrium.")
            else:
                print("Since K ≈ 1, reactants and products are present in similar amounts at equilibrium.")