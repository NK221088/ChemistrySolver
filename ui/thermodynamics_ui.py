"""
Terminal User Interface for Thermodynamics Calculations
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user
from chemistry_solver.thermodynamics import (calculate_heat, calculate_temperature_change, calculate_final_temperature,
                      calculate_molar_heat, solve_thermal_equilibrium, solve_thermal_equilibrium_with_molar_heat,
                      handle_heat_transfer_problem, handle_heat_transfer_with_molar_heat, solve_mixture_problem,
                      calculate_boiling_point_with_pressure, calculate_pressure_with_temperature, 
                      calculate_heat_of_vaporization, parse_reaction_string, solve_enthalpy_problem, 
                      solve_equilibrium_constant_problem, solve_equilibrium_from_formation_data, solve_combustion_enthalpy)
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
            print("Note: You can use fractions like '1/2 O2 (g)' for coefficients")
            
            known_reactions = []
            
            while True:
                try:
                    num_known = int(input("\nHow many known reactions do you want to enter? "))
                    break
                except ValueError:
                    print("Please enter a valid number.")
            
            for i in range(num_known):
                while True:
                    try:
                        reaction_str = input(f"Reaction {i+1}: ")
                        reaction = parse_reaction_string(reaction_str)
                        if reaction.enthalpy is None:
                            print("Warning: No enthalpy value found. Please include ΔH° = value")
                            continue
                        known_reactions.append(reaction)
                        print(f"Parsed: {reaction}")
                        break
                    except ValueError as e:
                        print(f"Error parsing reaction: {e}")
                        print("Please try again with the correct format.")
            
            print("\n--- Enter Target Reaction ---")
            print("Enter the reaction for which you want to calculate the enthalpy change:")
            print("(Don't include the ΔH° value - that's what we're solving for)")
            
            while True:
                try:
                    target_str = input("Target reaction: ")
                    target_reaction = parse_reaction_string(target_str)
                    target_reaction.enthalpy = None  # Make sure it's None
                    break
                except ValueError as e:
                    print(f"Error parsing reaction: {e}")
                    print("Please try again.")
            
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
        Simplified handler function for equilibrium constant calculations.
        """
        print("\n=== Equilibrium Constant Calculations ===")
        print("Select the type of calculation:")
        print("[1] Calculate equilibrium constant from ΔG° directly")
        print("[2] Calculate equilibrium constant from Gibbs free energies of formation")
        
        choice = input("\nEnter choice (1-2): ").strip()
        
        if choice == "1":
            print("\n--- Calculate Equilibrium Constant from ΔG° ---")
            
            # Get reaction string (optional)
            reaction = input("\nEnter balanced chemical equation (optional): ").strip() or None
            
            # Get temperature (optional, defaults to 25°C)
            temp_input = input("Enter temperature in °C (default=25): ").strip()
            temp_c = float(temp_input) if temp_input else 25
            
            # Get ΔG°
            delta_g = float(input("Enter ΔG° (kJ/mol): "))
            
            # Solve the problem
            result = solve_equilibrium_constant_problem(
                delta_g_standard=delta_g,
                temperature_c=temp_c,
                reaction_string=reaction
            )
            
            print("\n" + "="*50)
            for step in result["steps"]:
                print(step)
            
            print(f"\n=== FINAL RESULTS ===")
            print(f"Temperature: {result['temperature_c']:.1f} °C ({result['temperature_k']:.1f} K)")
            print(f"ΔG°: {result['delta_g_standard']:.1f} kJ/mol")
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
            print("This uses standard Gibbs free energies of formation (ΔGf°).")
            
            # Get reaction string and balance it using the existing balancer
            reaction = input("\nEnter chemical equation (unbalanced is fine): ").strip()
            
            # Try to balance the equation (keep your existing balancing logic)
            try:
                from chemistry_solver.balancer import parse_chemical_equation, balance_equation, format_balanced_equation
                
                print("\n--- Balancing Equation ---")
                reactants_parsed, products_parsed = parse_chemical_equation(reaction)
                balanced_reactants, balanced_products = balance_equation(reactants_parsed, products_parsed)
                balanced_equation = format_balanced_equation(balanced_reactants, balanced_products)
                
                print(f"Balanced equation: {balanced_equation}")
                
                # Convert to the format needed for calculations
                reactants = {formula: coeff for coeff, formula in balanced_reactants}
                products = {formula: coeff for coeff, formula in balanced_products}
                
                print(f"Reactants: {reactants}")
                print(f"Products: {products}")
                
                # Ask user to confirm
                confirm = input("\nIs this balanced equation correct? (y/n): ").lower()
                if confirm != 'y':
                    print("Please re-enter the equation or balance it manually.")
                    return
                    
            except ImportError as e:
                print(f"Balancer module not found: {e}")
                print("Please enter reactants and products manually:")
                reactants, products = self._get_manual_coefficients()
                balanced_equation = reaction
            except Exception as e:
                print(f"Error balancing equation: {e}")
                print("Please enter reactants and products manually:")
                reactants, products = self._get_manual_coefficients()
                balanced_equation = reaction
            
            # Get temperature
            temp_input = input("\nEnter temperature in °C (default=25): ").strip()
            temp_c = float(temp_input) if temp_input else 25
            
            reaction_coefficients = {
                'reactants': reactants,
                'products': products
            }
            
            # Get all unique substances
            all_substances = set(reactants.keys()) | set(products.keys())
            
            # Try to extract ΔGf° data from CSV file first
            formation_free_energies = {}
            
            try:
                thermo_data = self._extract_thermodynamic_data_from_csv(all_substances)
                
                if thermo_data:
                    print(f"\n--- Thermodynamic Data Found in Database ---")
                    print("The following ΔGf° data was found in NBS_Tables-Library.csv:")
                    
                    found_substances = []
                    missing_substances = []
                    
                    for substance in sorted(all_substances):
                        if substance in thermo_data and 'delta_f_g' in thermo_data[substance]:
                            formation_free_energies[substance] = thermo_data[substance]['delta_f_g']
                            found_substances.append(substance)
                            print(f"  {substance}: ΔGf° = {thermo_data[substance]['delta_f_g']:.1f} kJ/mol")
                        else:
                            missing_substances.append(substance)
                    
                    # Handle missing substances manually
                    if missing_substances:
                        print(f"\n--- Manual Entry Required ---")
                        print("The following substances were not found in the database:")
                        
                        for substance in missing_substances:
                            print(f"\nEntering ΔGf° for {substance}:")
                            hint = ""
                            if substance in ['H2', 'N2', 'O2', 'F2', 'Cl2', 'Br2', 'I2', 'C(graphite)', 'S(rhombic)']:
                                hint = " (hint: 0 for elements in standard state)"
                            
                            while True:
                                try:
                                    value_input = input(f"ΔGf° for {substance} (kJ/mol){hint}: ").strip()
                                    formation_free_energies[substance] = float(value_input)
                                    break
                                except ValueError:
                                    print("Please enter a valid number.")
                    
                    # Ask user to confirm the data
                    print(f"\n--- Confirm Gibbs Free Energy Data ---")
                    for substance in sorted(all_substances):
                        print(f"{substance}: ΔGf° = {formation_free_energies[substance]:.1f} kJ/mol")
                    
                    confirm = input("\nIs this thermodynamic data correct? (y/n): ").lower()
                    if confirm != 'y':
                        print("Please enter the data manually:")
                        formation_free_energies = self._get_manual_gibbs_data(all_substances)
                else:
                    raise ValueError("No data found in CSV")
                    
            except Exception as e:
                print(f"Could not extract ΔGf° data from CSV file: {e}")
                print("Please enter Gibbs free energy data manually:")
                formation_free_energies = self._get_manual_gibbs_data(all_substances)
            
            # Solve the problem
            result = solve_equilibrium_from_formation_data(
                reaction_coefficients=reaction_coefficients,
                formation_free_energies=formation_free_energies,
                temperature_c=temp_c,
                reaction_string=balanced_equation if 'balanced_equation' in locals() else reaction
            )
            
            print("\n" + "="*50)
            for step in result["steps"]:
                print(step)
            
            print(f"\n=== FINAL RESULTS ===")
            print(f"Reaction: {result['reaction']}")
            print(f"Temperature: {result['temperature_c']:.1f} °C ({result['temperature_k']:.1f} K)")
            print(f"ΔG°rxn: {result['delta_g_standard']:.1f} kJ/mol")
            print(f"Equilibrium constant (K): {result['equilibrium_constant']:.2e}")
            
            # Interpret the result
            if result['equilibrium_constant'] > 1:
                print("Since K > 1, the reaction favors products at equilibrium.")
            elif result['equilibrium_constant'] < 1:
                print("Since K < 1, the reaction favors reactants at equilibrium.")
            else:
                print("Since K ≈ 1, reactants and products are present in similar amounts at equilibrium.")
        
        else:
            print("Invalid choice. Please select 1 or 2.")
    def _extract_thermodynamic_data_from_csv(self, substances):
        """
        Extract thermodynamic data from the CSV format for given substances.
        Returns dict with structure: {formula: {'delta_f_h': value, 's_standard': value, 'delta_f_g': value}}
        
        CSV format expected:
        Substance,ΔfH°_kJ_mol-1,ΔfG°_kJ_mol-1,S°_J_mol-1_K-1
        """
        import pandas as pd
        import os
        import re
        
        csv_file = rf"ui\NBS_Tables-Library.csv"
        
        if not os.path.exists(csv_file):
            raise FileNotFoundError(f"CSV file {csv_file} not found")
        
        try:
            # Read the CSV file with header in first row (default behavior)
            df = pd.read_csv(csv_file)
            
            # Clean up column names by stripping whitespace
            df.columns = df.columns.str.strip()
            
            # Print available columns for debugging
            print(f"Available columns: {df.columns.tolist()}")
            
            # Define column mappings for the new format
            substance_col = 'Substance'
            enthalpy_col = '?fH�_kJ_mol-1'  # Formation enthalpy
            gibbs_col = '?fG�_kJ_mol-1'     # Formation Gibbs free energy
            entropy_col = 'S�_J_mol-1_K-1'  # Standard entropy
            
            # Alternative column names in case of encoding issues
            alt_enthalpy_cols = ['ΔfH°_kJ_mol-1', 'delta_fH_kJ_mol-1', 'DfH_kJ_mol-1']
            alt_gibbs_cols = ['ΔfG°_kJ_mol-1', 'delta_fG_kJ_mol-1', 'DfG_kJ_mol-1']
            alt_entropy_cols = ['S°_J_mol-1_K-1', 'S_J_mol-1_K-1', 'entropy_J_mol-1_K-1']
            
            # Find the correct column names if the default ones don't exist
            if enthalpy_col not in df.columns:
                for alt_col in alt_enthalpy_cols:
                    if alt_col in df.columns:
                        enthalpy_col = alt_col
                        break
                else:
                    # Try to find by partial match
                    for col in df.columns:
                        if 'fH' in col or 'enthalpy' in col.lower():
                            enthalpy_col = col
                            break
            
            if gibbs_col not in df.columns:
                for alt_col in alt_gibbs_cols:
                    if alt_col in df.columns:
                        gibbs_col = alt_col
                        break
                else:
                    # Try to find by partial match
                    for col in df.columns:
                        if 'fG' in col or ('gibbs' in col.lower() and 'formation' in col.lower()):
                            gibbs_col = col
                            break
            
            if entropy_col not in df.columns:
                for alt_col in alt_entropy_cols:
                    if alt_col in df.columns:
                        entropy_col = alt_col
                        break
                else:
                    # Try to find by partial match  
                    for col in df.columns:
                        if 'S°' in col or 'entropy' in col.lower() or ('J_mol' in col and 'K' in col):
                            entropy_col = col
                            break
            
            # Verify we found required columns
            if substance_col not in df.columns:
                raise ValueError(f"Substance column '{substance_col}' not found in CSV")
            
            # For flexibility, we'll work with whatever columns are available
            available_data_types = []
            if enthalpy_col in df.columns:
                available_data_types.append(f"Enthalpy: '{enthalpy_col}'")
            if gibbs_col in df.columns:
                available_data_types.append(f"Gibbs: '{gibbs_col}'")
            if entropy_col in df.columns:
                available_data_types.append(f"Entropy: '{entropy_col}'")
            
            print(f"Using columns - Substance: '{substance_col}', " + ", ".join(available_data_types))
            
            # Clean the dataframe - remove rows where Substance is NaN
            df_clean = df.dropna(subset=[substance_col])
            
            # Extract data for requested substances
            thermo_data = {}
            
            for substance in substances:
                # Look for exact matches first (case-sensitive)
                matches = df_clean[df_clean[substance_col].str.strip() == substance.strip()]
                
                if matches.empty:
                    # Try case-insensitive match
                    matches = df_clean[df_clean[substance_col].str.strip().str.upper() == substance.strip().upper()]
                
                if matches.empty:
                    # Try partial match, but be more specific to avoid false matches
                    # Only do partial matching for substances with parentheses (like states)
                    if '(' in substance:
                        # For substances with states, try matching without being too broad
                        matches = df_clean[df_clean[substance_col].str.contains(re.escape(substance.strip()), case=False, na=False, regex=True)]
                    else:
                        # For simple substances, try exact word boundary matching to avoid H2 matching H2CO3
                        pattern = r'\b' + re.escape(substance.strip()) + r'\b'
                        matches = df_clean[df_clean[substance_col].str.contains(pattern, case=False, na=False, regex=True)]
                
                if not matches.empty:
                    # If multiple matches, take the first one
                    row = matches.iloc[0]
                    
                    try:
                        substance_data = {}
                        
                        # Extract enthalpy (ΔfH°) if available
                        if enthalpy_col in df.columns:
                            delta_f_h_raw = row[enthalpy_col]
                            if pd.isna(delta_f_h_raw) or str(delta_f_h_raw).strip() in ['-', '', 'NaN']:
                                delta_f_h = 0.0
                            else:
                                delta_f_h = float(str(delta_f_h_raw).replace(',', '').strip())
                            substance_data['delta_f_h'] = delta_f_h
                        
                        # Extract Gibbs free energy (ΔfG°) if available
                        if gibbs_col in df.columns:
                            delta_f_g_raw = row[gibbs_col]
                            if pd.isna(delta_f_g_raw) or str(delta_f_g_raw).strip() in ['-', '', 'NaN']:
                                delta_f_g = 0.0
                            else:
                                delta_f_g = float(str(delta_f_g_raw).replace(',', '').strip())
                            substance_data['delta_f_g'] = delta_f_g
                        
                        # Extract entropy (S°) if available
                        if entropy_col in df.columns:
                            s_standard_raw = row[entropy_col]
                            if pd.isna(s_standard_raw) or str(s_standard_raw).strip() in ['-', '', 'NaN']:
                                s_standard = 0.0
                            else:
                                s_standard = float(str(s_standard_raw).replace(',', '').strip())
                            substance_data['s_standard'] = s_standard
                        
                        if substance_data:  # Only add if we found some data
                            thermo_data[substance] = substance_data
                            
                            data_str = []
                            if 'delta_f_h' in substance_data:
                                data_str.append(f"ΔfH° = {substance_data['delta_f_h']} kJ/mol")
                            if 'delta_f_g' in substance_data:
                                data_str.append(f"ΔfG° = {substance_data['delta_f_g']} kJ/mol")
                            if 's_standard' in substance_data:
                                data_str.append(f"S° = {substance_data['s_standard']} J/(mol·K)")
                            
                            print(f"Found data for {substance}: {', '.join(data_str)}")
                            
                    except (ValueError, TypeError) as e:
                        print(f"Could not parse data for {substance}: {e}")
                        if enthalpy_col in df.columns:
                            print(f"Raw enthalpy value: {row[enthalpy_col]}")
                        if gibbs_col in df.columns:
                            print(f"Raw Gibbs value: {row[gibbs_col]}")
                        if entropy_col in df.columns:
                            print(f"Raw entropy value: {row[entropy_col]}")
                        continue
                else:
                    print(f"No data found for substance: {substance}")
                    print(f"Available substances in CSV: {df_clean[substance_col].head(10).tolist()}")
            
            return thermo_data
            
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            print(f"Available columns: {df.columns.tolist() if 'df' in locals() else 'Could not read CSV'}")
            raise Exception(f"Error reading CSV file: {e}")
    def _get_manual_coefficients(self):
        """
        Helper method to manually collect reactants and products when balancer fails.
        """
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
        
        return reactants, products

    def _get_manual_thermodynamic_data(self, substances):
        """
        Helper method to manually collect thermodynamic data when CSV extraction fails.
        """
        print(f"\n--- Enter Formation Enthalpies (ΔHf°) ---")
        print("Note: Elements in their standard state have ΔHf° = 0")
        formation_enthalpies = {}
        
        for substance in sorted(substances):
            default_hint = " (0 for elements)" if substance in ['H2', 'N2', 'O2', 'F2', 'Cl2', 'Br2', 'I2'] else ""
            value = input(f"ΔHf° for {substance} (kJ/mol){default_hint}: ")
            formation_enthalpies[substance] = float(value) if value else 0.0
        
        print(f"\n--- Enter Standard Entropies (S°) ---")
        standard_entropies = {}
        
        for substance in sorted(substances):
            value = float(input(f"S° for {substance} (J/(mol·K)): "))
            standard_entropies[substance] = value
        
        return formation_enthalpies, standard_entropies