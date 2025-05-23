"""
Terminal User Interface for Colligative Properties Calculations
Updated to work with the refactored ColligativePropertyCalculator
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user, display_steps
from chemistry_solver.colligative_properties import (
    ColligativePropertyCalculator,
    calculate_molecular_weight,
    compare_colligative_properties,
    solve_multiple_choice_problem
)
from chemistry_solver.colligative_constants import SOLVENT_CONSTANTS


class ColligativePropertiesUI:
    """UI class for colligative properties calculations."""
    
    def __init__(self):
        self.title = "COLLIGATIVE PROPERTIES CALCULATOR"
    
    def run(self):
        """Run the colligative properties UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-9): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_molecular_weight_calculation()
            elif choice == "2":
                self._handle_solution_comparison()
            elif choice == "3":
                self._handle_multiple_choice_solver()
            elif choice == "4":
                self._handle_freezing_point_depression()
            elif choice == "5":
                self._handle_boiling_point_elevation()
            elif choice == "6":
                self._handle_osmotic_pressure()
            elif choice == "7":
                self._handle_vapor_pressure_lowering()
            elif choice == "8":
                self._handle_property_from_molality()
            elif choice == "9":
                self._display_constants()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the colligative properties module menu."""
        menu = """
        🧪 COLLIGATIVE PROPERTIES CALCULATOR
        
        [1] Calculate Molecular Weight (Unified Calculator)
        [2] Compare Solutions (Colligative Properties)
        [3] Multiple Choice Problem Solver
        [4] Freezing Point Depression (Legacy)
        [5] Boiling Point Elevation (Legacy)
        [6] Osmotic Pressure (Legacy)
        [7] Vapor Pressure Lowering (Legacy)
        [8] Calculate Property from Molality
        [9] Display Solvent Constants
        [0] Return to main menu
        """
        print(menu)
    
    def _get_solvent_selection(self):
        """Helper to select solvent from available constants."""
        print("\nAvailable solvents:")
        solvents = list(SOLVENT_CONSTANTS.keys())
        for i, solvent in enumerate(solvents, 1):
            constants = SOLVENT_CONSTANTS[solvent]
            print(f"  [{i}] {solvent.title()} (Kf={constants['Kf']}, Kb={constants['Kb']})")
        
        while True:
            try:
                choice = input(f"Select solvent (1-{len(solvents)}) or enter custom name: ").strip()
                if choice.isdigit():
                    idx = int(choice) - 1
                    if 0 <= idx < len(solvents):
                        return solvents[idx]
                else:
                    return choice.lower()
            except (ValueError, IndexError):
                print("Invalid selection. Please try again.")
    
    def _get_ionization_factor(self):
        """Helper method to get van't Hoff factor."""
        print("\nIonization Factor (van't Hoff factor):")
        print("  1 = Non-electrolyte (e.g., glucose, urea)")
        print("  2 = Strong 1:1 electrolyte (e.g., NaCl, KBr)")
        print("  3 = Strong 1:2 electrolyte (e.g., CaCl₂, Na₂SO₄)")
        print("  4 = Strong 1:3 electrolyte (e.g., AlCl₃)")
        
        while True:
            try:
                factor = float(input("Enter van't Hoff factor (i): "))
                if factor > 0:
                    return factor
                else:
                    print("van't Hoff factor must be positive.")
            except ValueError:
                print("Please enter a valid number.")
    
    def _handle_molecular_weight_calculation(self):
        """Handle unified molecular weight calculation."""
        print("\n" + "="*70)
        print("MOLECULAR WEIGHT CALCULATOR")
        print("="*70)
        print("Calculate molecular weight from colligative property data")
        print()
        
        # Select method
        print("Select calculation method:")
        print("1. Freezing Point Depression")
        print("2. Boiling Point Elevation") 
        print("3. Osmotic Pressure")
        print("4. Vapor Pressure Lowering")
        
        method_choice = input("Choice (1-4): ").strip()
        method_map = {
            "1": "freezing_point",
            "2": "boiling_point", 
            "3": "osmotic_pressure",
            "4": "vapor_pressure"
        }
        
        if method_choice not in method_map:
            print("Invalid choice.")
            return
        
        method = method_map[method_choice]
        
        try:
            # Get common parameters
            solvent = self._get_solvent_selection()
            
            if method in ["freezing_point", "boiling_point"]:
                self._handle_temperature_method(method, solvent)
            elif method == "osmotic_pressure":
                self._handle_osmotic_method(solvent)
            elif method == "vapor_pressure":
                self._handle_vapor_pressure_method(solvent)
                
        except ValueError as e:
            print(f"❌ Input Error: {str(e)}")
        except Exception as e:
            print(f"❌ Error: {str(e)}")
    
    def _handle_temperature_method(self, method, solvent):
        """Handle freezing point or boiling point methods."""
        property_name = "freezing point depression" if method == "freezing_point" else "boiling point elevation"
        symbol = "ΔTf" if method == "freezing_point" else "ΔTb"
        
        print(f"\n{property_name.upper()} METHOD")
        print(f"Formula: {symbol} = K{'f' if method == 'freezing_point' else 'b'} × m × i")
        
        # Get temperature change
        delta_T = float(input(f"Enter {property_name} ({symbol} in °C): "))
        
        # Get mass data
        solute_mass = float(input("Enter mass of solute (g): "))
        solvent_mass = float(input("Enter mass of solvent (g): "))
        
        # Get ionization factor
        ionization_factor = self._get_ionization_factor()
        
        # Get answer choices if available
        has_choices = input("\nDo you have multiple choice answers? (y/n): ").lower().startswith('y')
        answer_choices = None
        if has_choices:
            choices_input = input("Enter choices separated by commas: ")
            answer_choices = [float(x.strip()) for x in choices_input.split(',')]
        
        # Calculate
        result = calculate_molecular_weight(
            method=method,
            solvent=solvent,
            delta_T=delta_T,
            solute_mass=solute_mass,
            solvent_mass=solvent_mass,
            ionization_factor=ionization_factor,
            answer_choices=answer_choices
        )
        
        self._display_calculation_result(result, f"Molecular Weight from {property_name.title()}")
    
    def _handle_osmotic_method(self, solvent):
        """Handle osmotic pressure method."""
        print("\nOSMOTIC PRESSURE METHOD")
        print("Formula: π = iMRT (M = molarity)")
        
        # Get osmotic pressure
        pressure_input = input("Enter osmotic pressure (include units - atm or torr): ").strip().lower()
        if 'torr' in pressure_input:
            pressure_torr = float(pressure_input.replace('torr', '').strip())
            osmotic_pressure_atm = pressure_torr / 760
            print(f"Converted: {pressure_torr} torr = {osmotic_pressure_atm:.6f} atm")
        else:
            osmotic_pressure_atm = float(pressure_input.replace('atm', '').strip())
        
        # Get other parameters
        temperature_c = float(input("Enter temperature (°C): "))
        solution_volume_L = float(input("Enter solution volume (L): "))
        solute_mass = float(input("Enter mass of solute (g): "))
        ionization_factor = self._get_ionization_factor()
        
        # Get answer choices
        has_choices = input("\nDo you have multiple choice answers? (y/n): ").lower().startswith('y')
        answer_choices = None
        if has_choices:
            choices_input = input("Enter choices separated by commas: ")
            answer_choices = [float(x.strip()) for x in choices_input.split(',')]
        
        # Calculate
        result = calculate_molecular_weight(
            method="osmotic_pressure",
            solvent=solvent,
            delta_T=osmotic_pressure_atm,  # Using delta_T parameter for pressure
            solute_mass=solute_mass,
            solvent_mass=0,  # Not used for osmotic pressure
            ionization_factor=ionization_factor,
            answer_choices=answer_choices,
            temperature_c=temperature_c,
            solution_volume_L=solution_volume_L
        )
        
        self._display_calculation_result(result, "Molecular Weight from Osmotic Pressure")
    
    def _handle_vapor_pressure_method(self, solvent):
        """Handle vapor pressure lowering method."""
        print("\nVAPOR PRESSURE LOWERING METHOD")
        print("Formula: P_solution = P_pure × (1 - X_solute)")
        
        # Get pressure data
        P_pure = float(input("Enter vapor pressure of pure solvent: "))
        P_solution = float(input("Enter vapor pressure of solution: "))
        
        # Get mass data
        solute_mass = float(input("Enter mass of solute (g): "))
        solvent_formula = input("Enter solvent formula (e.g., H2O, C6H6): ").strip()
        solvent_mass = float(input("Enter mass of solvent (g): "))
        
        # Get answer choices
        has_choices = input("\nDo you have multiple choice answers? (y/n): ").lower().startswith('y')
        answer_choices = None
        if has_choices:
            choices_input = input("Enter choices separated by commas: ")
            answer_choices = [float(x.strip()) for x in choices_input.split(',')]
        
        # Calculate
        result = calculate_molecular_weight(
            method="vapor_pressure",
            solvent=solvent,
            delta_T=0,  # Not used
            solute_mass=solute_mass,
            solvent_mass=solvent_mass,
            answer_choices=answer_choices,
            P_pure=P_pure,
            P_solution=P_solution,
            solvent_formula=solvent_formula
        )
        
        self._display_calculation_result(result, "Molecular Weight from Vapor Pressure Lowering")
    
    def _handle_solution_comparison(self):
        """Handle solution comparison for colligative properties."""
        print("\n" + "="*70)
        print("SOLUTION COMPARISON")
        print("="*70)
        print("Compare colligative properties of multiple solutions")
        print()
        
        try:
            # Select solvent
            solvent = self._get_solvent_selection()
            
            # Select property type
            print("\nSelect property to compare:")
            print("1. Freezing Point Depression")
            print("2. Boiling Point Elevation")
            
            prop_choice = input("Choice (1-2): ").strip()
            if prop_choice == "1":
                property_type = "freezing_point_depression"
            elif prop_choice == "2":
                property_type = "boiling_point_elevation"
            else:
                print("Invalid choice.")
                return
            
            # Get number of solutions
            num_solutions = int(input("How many solutions to compare? "))
            
            solutions = []
            for i in range(num_solutions):
                print(f"\nSolution {i+1}:")
                name = input(f"  Name (optional, press Enter for default): ").strip()
                if not name:
                    name = f"Solution {i+1}"
                
                molality = float(input(f"  Molality (mol/kg): "))
                
                # Ask about ionization
                is_electrolyte = input(f"  Is this an electrolyte? (y/n): ").lower().startswith('y')
                ionization_factor = 1
                if is_electrolyte:
                    ionization_factor = self._get_ionization_factor()
                
                solutions.append({
                    'name': f"{name} ({molality} m)",
                    'molality': molality,
                    'ionization_factor': ionization_factor
                })
            
            # Compare solutions
            result = compare_colligative_properties(
                solutions=solutions,
                property_type=property_type,
                solvent=solvent
            )
            
            if result['success']:
                print(f"\n{'='*70}")
                print("SOLUTION COMPARISON RESULTS")
                print(f"{'='*70}")
                
                display_steps(result['steps'])
                
                print(f"\n🏆 HIGHEST EFFECT: {result['highest']['name']}")
                print(f"   Property Value: {result['highest']['property_value']:.3f}°C")
                
            else:
                print(f"❌ Error: {result['error']}")
                
        except ValueError as e:
            print(f"❌ Input Error: {str(e)}")
        except Exception as e:
            print(f"❌ Error: {str(e)}")
    
    def _handle_multiple_choice_solver(self):
        """Handle multiple choice problem solver."""
        print("\n" + "="*70)
        print("MULTIPLE CHOICE PROBLEM SOLVER")
        print("="*70)
        print("Structured solver for multiple choice problems")
        print()
        
        try:
            print("Select problem type:")
            print("1. Molecular Weight Calculation")
            print("2. Solution Comparison")
            
            problem_choice = input("Choice (1-2): ").strip()
            
            if problem_choice == "1":
                self._handle_mc_molecular_weight()
            elif problem_choice == "2":
                self._handle_mc_solution_comparison()
            else:
                print("Invalid choice.")
                
        except Exception as e:
            print(f"❌ Error: {str(e)}")
    
    def _handle_mc_molecular_weight(self):
        """Handle multiple choice molecular weight problems."""
        print("\nMOLECULAR WEIGHT MULTIPLE CHOICE")
        
        # Get method
        print("Select method:")
        print("1. Freezing Point Depression")
        print("2. Boiling Point Elevation")
        print("3. Osmotic Pressure")
        print("4. Vapor Pressure Lowering")
        
        method_choice = input("Choice (1-4): ").strip()
        method_map = {
            "1": "freezing_point",
            "2": "boiling_point",
            "3": "osmotic_pressure", 
            "4": "vapor_pressure"
        }
        
        if method_choice not in method_map:
            print("Invalid choice.")
            return
        
        method = method_map[method_choice]
        
        # Build problem data based on method
        problem_data = {
            'type': 'molecular_weight',
            'method': method,
            'solvent': self._get_solvent_selection()
        }
        
        # Get method-specific parameters
        if method in ["freezing_point", "boiling_point"]:
            property_name = "freezing point depression" if method == "freezing_point" else "boiling point elevation"
            problem_data['delta_T'] = float(input(f"Enter {property_name} (°C): "))
            problem_data['solute_mass'] = float(input("Enter solute mass (g): "))
            problem_data['solvent_mass'] = float(input("Enter solvent mass (g): "))
            problem_data['ionization_factor'] = self._get_ionization_factor()
            
        elif method == "osmotic_pressure":
            pressure_atm = float(input("Enter osmotic pressure (atm): "))
            problem_data['delta_T'] = pressure_atm  # Reusing parameter
            problem_data['solute_mass'] = float(input("Enter solute mass (g): "))
            problem_data['temperature_c'] = float(input("Enter temperature (°C): "))
            problem_data['solution_volume_L'] = float(input("Enter solution volume (L): "))
            problem_data['ionization_factor'] = self._get_ionization_factor()
            
        elif method == "vapor_pressure":
            problem_data['P_pure'] = float(input("Enter pure solvent vapor pressure: "))
            problem_data['P_solution'] = float(input("Enter solution vapor pressure: "))
            problem_data['solute_mass'] = float(input("Enter solute mass (g): "))
            problem_data['solvent_formula'] = input("Enter solvent formula: ").strip()
            problem_data['solvent_mass'] = float(input("Enter solvent mass (g): "))
        
        # Get answer choices
        choices_input = input("Enter answer choices separated by commas: ")
        problem_data['answer_choices'] = [float(x.strip()) for x in choices_input.split(',')]
        
        # Solve problem
        result = solve_multiple_choice_problem(problem_data)
        
        self._display_calculation_result(result, "Multiple Choice Problem Solution")
    
    def _handle_mc_solution_comparison(self):
        """Handle multiple choice solution comparison problems."""
        print("\nSOLUTION COMPARISON MULTIPLE CHOICE")
        
        # Get basic setup
        solvent = self._get_solvent_selection()
        
        print("Select property:")
        print("1. Freezing Point Depression")
        print("2. Boiling Point Elevation")
        
        prop_choice = input("Choice (1-2): ").strip()
        property_type = "freezing_point_depression" if prop_choice == "1" else "boiling_point_elevation"
        
        # Get solutions
        num_solutions = int(input("Number of solutions: "))
        solutions = []
        
        for i in range(num_solutions):
            print(f"\nSolution {i+1}:")
            molality = float(input("  Molality: "))
            ionization_factor = float(input("  van't Hoff factor (1 for non-electrolyte): "))
            
            solutions.append({
                'molality': molality,
                'ionization_factor': ionization_factor,
                'name': f"Solution {i+1} ({molality}m, i={ionization_factor})"
            })
        
        # Solve
        problem_data = {
            'type': 'compare_solutions',
            'solvent': solvent,
            'property_type': property_type,
            'solutions': solutions
        }
        
        result = solve_multiple_choice_problem(problem_data)
        
        if result['success']:
            print(f"\n{'='*70}")
            print("COMPARISON RESULTS")
            print(f"{'='*70}")
            
            display_steps(result['steps'])
            
            if 'highest' in result:
                print(f"\n🏆 ANSWER: {result['highest']['name']}")
                print(f"   Has the highest {property_type.replace('_', ' ')}")
        else:
            print(f"❌ Error: {result['error']}")
    
    def _handle_property_from_molality(self):
        """Calculate colligative property values from molality."""
        print("\n" + "="*70)
        print("CALCULATE PROPERTY FROM MOLALITY")
        print("="*70)
        print("Calculate colligative property values when you know the molality")
        print()
        
        try:
            # Get parameters
            solvent = self._get_solvent_selection()
            calculator = ColligativePropertyCalculator(solvent)
            
            molality = float(input("Enter molality (mol/kg): "))
            ionization_factor = self._get_ionization_factor()
            
            print("\nSelect property to calculate:")
            print("1. Freezing Point Depression")
            print("2. Boiling Point Elevation")
            
            prop_choice = input("Choice (1-2): ").strip()
            
            if prop_choice == "1":
                property_type = "freezing_point_depression"
                symbol = "ΔTf"
            elif prop_choice == "2":
                property_type = "boiling_point_elevation"
                symbol = "ΔTb"
            else:
                print("Invalid choice.")
                return
            
            # Calculate
            property_value = calculator.calculate_property_from_molality(
                molality, property_type, ionization_factor
            )
            
            # Display result
            print(f"\n{'='*50}")
            print(f"CALCULATION RESULT")
            print(f"{'='*50}")
            print(f"Solvent: {solvent.title()}")
            print(f"Molality: {molality} mol/kg")
            print(f"van't Hoff factor: {ionization_factor}")
            print(f"")
            print(f"{symbol} = K × m × i")
            constant = calculator.constants["Kf"] if property_type == "freezing_point_depression" else calculator.constants["Kb"]
            print(f"{symbol} = {constant} × {molality} × {ionization_factor}")
            print(f"{symbol} = {property_value:.4f} °C")
            
        except ValueError as e:
            print(f"❌ Input Error: {str(e)}")
        except Exception as e:
            print(f"❌ Error: {str(e)}")
    
    # Legacy methods for backward compatibility
    def _handle_freezing_point_depression(self):
        """Legacy freezing point depression interface."""
        print("\n⚠️  Using legacy interface. Consider using option [1] for unified calculator.")
        self._handle_temperature_method("freezing_point", "water")
    
    def _handle_boiling_point_elevation(self):
        """Legacy boiling point elevation interface."""
        print("\n⚠️  Using legacy interface. Consider using option [1] for unified calculator.")
        self._handle_temperature_method("boiling_point", "water")
    
    def _handle_osmotic_pressure(self):
        """Legacy osmotic pressure interface."""
        print("\n⚠️  Using legacy interface. Consider using option [1] for unified calculator.")
        self._handle_osmotic_method("water")
    
    def _handle_vapor_pressure_lowering(self):
        """Legacy vapor pressure lowering interface."""
        print("\n⚠️  Using legacy interface. Consider using option [1] for unified calculator.")
        self._handle_vapor_pressure_method("water")
    
    def _display_calculation_result(self, result, title):
        """Display calculation results in a formatted way."""
        print(f"\n{'='*70}")
        print(f"{title.upper()}")
        print(f"{'='*70}")
        
        if result['success']:
            # Display key results
            if 'molecular_weight' in result:
                print(f"📊 Molecular Weight: {result['molecular_weight']:.2f} g/mol")
            
            if 'closest_answer' in result:
                print(f"🎯 Closest Answer Choice: {result['closest_answer']}")
                if 'answer_difference' in result:
                    print(f"   Difference: ±{result['answer_difference']:.4f}")
            
            # Display detailed steps
            print(f"\n{'DETAILED CALCULATION STEPS'}")
            print(f"{'-'*70}")
            display_steps(result['steps'])
            
        else:
            print(f"❌ ERROR: {result['error']}")
    
    def _display_constants(self):
        """Display available solvent constants."""
        print("\n" + "="*80)
        print("SOLVENT CONSTANTS DATABASE")
        print("="*80)
        
        print(f"\n{'Solvent':<15} {'Kf (°C/m)':<12} {'Kb (°C/m)':<12} {'Notes'}")
        print("-" * 80)
        
        for solvent, constants in SOLVENT_CONSTANTS.items():
            kf = constants.get('Kf', 'N/A')
            kb = constants.get('Kb', 'N/A')
            
            # Format values
            kf_str = f"{kf:.3f}" if isinstance(kf, (int, float)) else str(kf)
            kb_str = f"{kb:.3f}" if isinstance(kb, (int, float)) else str(kb)
            
            print(f"{solvent.title():<15} {kf_str:<12} {kb_str:<12}")
        
        print(f"\n📚 Usage Notes:")
        print(f"• Kf = Freezing point depression constant")
        print(f"• Kb = Boiling point elevation constant") 
        print(f"• Formulas: ΔTf = Kf×m×i and ΔTb = Kb×m×i")
        print(f"• m = molality (mol solute/kg solvent)")
        print(f"• i = van't Hoff factor (ionization factor)")
        print(f"• Always verify constants from reliable sources for precise work")