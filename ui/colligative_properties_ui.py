"""
Terminal User Interface for Colligative Properties Calculations
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user, display_steps
from chemistry_solver.colligative_properties import (
    solve_freezing_point_problem,
    solve_boiling_point_problem,
    solve_osmotic_pressure_problem,
    solve_vapor_pressure_problem,
    solve_molecular_weight_multiple_choice,
    calculate_molecular_weight_from_freezing_point,
    calculate_molecular_weight_from_boiling_point,
    FREEZING_POINT_CONSTANTS,
    BOILING_POINT_CONSTANTS
)

class ColligativePropertiesUI:
    """UI class for colligative properties calculations."""
    
    def __init__(self):
        self.title = "COLLIGATIVE PROPERTIES CALCULATOR"
    
    def run(self):
        """Run the colligative properties UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-8): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_freezing_point_depression()
            elif choice == "2":
                self._handle_boiling_point_elevation()
            elif choice == "3":
                self._handle_osmotic_pressure()
            elif choice == "4":
                self._handle_vapor_pressure_lowering()
            elif choice == "5":
                self._handle_molecular_weight_from_freezing_point()  # NEW
            elif choice == "6":
                self._handle_molecular_weight_from_boiling_point()  # NEW
            elif choice == "7":
                self._handle_multiple_choice_solver()  # NEW
            elif choice == "8":
                self._display_constants()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the colligative properties module menu."""
        menu = """
        [1] Freezing Point Depression (Calculate MW from temperatures)
        [2] Boiling Point Elevation (Calculate MW from temperatures)
        [3] Osmotic Pressure (Calculate MW)
        [4] Vapor Pressure Lowering (Calculate MW)
        [5] Molecular Weight from Freezing Point Data
        [6] Molecular Weight from Boiling Point Data
        [7] Multiple Choice Problem Solver
        [8] Display Constants
        [0] Return to main menu
        """
        print(menu)
    
    def _get_solvent_info(self, constants_dict, constant_name):
        """Helper method to get solvent information and constant."""
        solvent = input("Enter solvent name or formula (e.g., water, H2O): ").strip()
        
        # Ask if user wants to use known constants
        use_lookup = input("Use known constant from database? (y/n): ").lower().startswith('y')
        constant_value = None
        
        if use_lookup:
            print(f"\nAvailable solvents with known {constant_name} constants:")
            solvent_list = list(constants_dict.items())
            for i, (solvent_name, constant) in enumerate(solvent_list, 1):
                print(f"  [{i}] {solvent_name.title()} ({constant_name} = {constant} °C/m)")
            
            choice = input(f"\nSelect a solvent (1-{len(solvent_list)}) or press Enter to continue with '{solvent}': ").strip()
            if choice.isdigit():
                idx = int(choice) - 1
                if 0 <= idx < len(solvent_list):
                    selected_solvent, constant_value = solvent_list[idx]
                    solvent = selected_solvent
                    print(f"Selected: {selected_solvent.title()} with {constant_name} = {constant_value} °C/m")
        
        # If no constant selected, ask user to provide it
        if constant_value is None:
            constant_input = input(f"Enter {constant_name} constant for {solvent} (°C/m) [leave blank for auto-lookup]: ").strip()
            if constant_input:
                constant_value = float(constant_input)
        
        return solvent, constant_value
    
    def _get_ionization_factor(self):
        """Helper method to get van't Hoff factor."""
        ionization = input("Is the solute an electrolyte? (y/n): ").lower().startswith('y')
        ionization_factor = 1
        if ionization:
            while True:
                try:
                    ionization_factor = float(input("Enter van't Hoff factor (i) [typical values: 1=nonelectrolyte, 2=NaCl, 3=CaCl2]: "))
                    if ionization_factor > 0:
                        break
                    else:
                        print("van't Hoff factor must be positive.")
                except ValueError:
                    print("Please enter a valid number.")
        return ionization_factor
    
    def _display_result_summary(self, result, calculation_type):
        """Display a summary of calculation results."""
        print(f"\n{'='*60}")
        print(f"CALCULATION SUMMARY - {calculation_type.upper()}")
        print(f"{'='*60}")
        
        if result['success']:
            # Display key results
            print(f"Molecular Weight: {result['molecular_weight']:.2f} g/mol")
            if 'closest_answer' in result:
                print(f"Closest Common Answer: {result['closest_answer']} g/mol")
            
            # Display specific properties based on calculation type
            if 'delta_T' in result:
                temp_change = "Depression" if calculation_type == "freezing point" else "Elevation"
                print(f"Temperature {temp_change} (ΔT): {result['delta_T']:.4f} °C")
            if 'molality' in result:
                print(f"Molality: {result['molality']:.6f} mol/kg")
            if 'moles_solute' in result:
                print(f"Moles of Solute: {result['moles_solute']:.6f} mol")
            if 'mole_fraction_solute' in result:
                print(f"Mole Fraction of Solute: {result['mole_fraction_solute']:.6f}")
            if 'delta_P' in result:
                print(f"Vapor Pressure Lowering (ΔP): {result['delta_P']:.4f}")
            
            print(f"\n{'DETAILED CALCULATION STEPS'}")
            print(f"{'-'*60}")
            display_steps(result['steps'])
        else:
            print(f"❌ ERROR: {result['error']}")
    
    def _handle_freezing_point_depression(self):
        """Handle freezing point depression calculations."""
        print("\n" + "="*70)
        print("FREEZING POINT DEPRESSION CALCULATOR")
        print("="*70)
        print("This calculator determines molecular weight from freezing point depression.")
        print("Formula: ΔTf = Kf × m × i")
        print()
        
        try:
            # Get temperature data
            T_pure = float(input("Enter freezing point of pure solvent (°C): "))
            T_solution = float(input("Enter freezing point of solution (°C): "))
            
            if T_solution >= T_pure:
                print("⚠️  Warning: Solution freezing point should be lower than pure solvent.")
            
            # Get solvent information
            solvent, K_f = self._get_solvent_info(FREEZING_POINT_CONSTANTS, "Kf")
            
            # Get mass data
            solute_mass = float(input("Enter mass of solute (g): "))
            solvent_mass = float(input("Enter mass of solvent (g): "))
            
            # Get ionization factor
            ionization_factor = self._get_ionization_factor()
            
            # Solve the problem
            result = solve_freezing_point_problem(
                T_pure=T_pure,
                T_solution=T_solution,
                solvent=solvent,
                K_f=K_f,
                solute_mass=solute_mass,
                solvent_mass=solvent_mass,
                ionization_factor=ionization_factor
            )
            
            self._display_result_summary(result, "freezing point depression")
                
        except ValueError as e:
            print(f"❌ Input Error: Please enter valid numerical values. ({str(e)})")
        except Exception as e:
            print(f"❌ Error: {str(e)}")
    
    def _handle_boiling_point_elevation(self):
        """Handle boiling point elevation calculations."""
        print("\n" + "="*70)
        print("BOILING POINT ELEVATION CALCULATOR")
        print("="*70)
        print("This calculator determines molecular weight from boiling point elevation.")
        print("Formula: ΔTb = Kb × m × i")
        print()
        
        try:
            # Get temperature data
            T_pure = float(input("Enter boiling point of pure solvent (°C): "))
            T_solution = float(input("Enter boiling point of solution (°C): "))
            
            if T_solution <= T_pure:
                print("⚠️  Warning: Solution boiling point should be higher than pure solvent.")
            
            # Get solvent information
            solvent, K_b = self._get_solvent_info(BOILING_POINT_CONSTANTS, "Kb")
            
            # Get mass data
            solute_mass = float(input("Enter mass of solute (g): "))
            solvent_mass = float(input("Enter mass of solvent (g): "))
            
            # Get ionization factor
            ionization_factor = self._get_ionization_factor()
            
            # Solve the problem
            result = solve_boiling_point_problem(
                T_pure=T_pure,
                T_solution=T_solution,
                solvent=solvent,
                K_b=K_b,
                solute_mass=solute_mass,
                solvent_mass=solvent_mass,
                ionization_factor=ionization_factor
            )
            
            self._display_result_summary(result, "boiling point elevation")
                
        except ValueError as e:
            print(f"❌ Input Error: Please enter valid numerical values. ({str(e)})")
        except Exception as e:
            print(f"❌ Error: {str(e)}")
    
    def _handle_osmotic_pressure(self):
        """Handle osmotic pressure calculations."""
        print("\n" + "="*70)
        print("OSMOTIC PRESSURE CALCULATOR")
        print("="*70)
        print("This calculator determines molecular weight from osmotic pressure.")
        print("Formula: π = iMRT (where M = molarity)")
        print()
        
        try:
            # Get pressure and temperature data
            pressure_input = input("Enter osmotic pressure (atm) or (torr): ").strip()
            
            # Handle different pressure units
            if 'torr' in pressure_input.lower():
                pressure_torr = float(pressure_input.replace('torr', '').strip())
                osmotic_pressure_atm = pressure_torr / 760  # Convert torr to atm
                print(f"Converted: {pressure_torr} torr = {osmotic_pressure_atm:.6f} atm")
            else:
                osmotic_pressure_atm = float(pressure_input.replace('atm', '').strip())
            
            temperature_c = float(input("Enter temperature (°C): "))
            solution_volume_L = float(input("Enter solution volume (L): "))
            solute_mass = float(input("Enter mass of solute (g): "))
            
            # Get ionization factor
            ionization_factor = self._get_ionization_factor()
            
            # Solve the problem
            result = solve_osmotic_pressure_problem(
                osmotic_pressure_atm=osmotic_pressure_atm,
                temperature_c=temperature_c,
                solution_volume_L=solution_volume_L,
                solute_mass=solute_mass,
                ionization_factor=ionization_factor
            )
            
            self._display_result_summary(result, "osmotic pressure")
                
        except ValueError as e:
            print(f"❌ Input Error: Please enter valid numerical values. ({str(e)})")
        except Exception as e:
            print(f"❌ Error: {str(e)}")
    
    def _handle_vapor_pressure_lowering(self):
        """Handle vapor pressure lowering calculations."""
        print("\n" + "="*70)
        print("VAPOR PRESSURE LOWERING CALCULATOR")
        print("="*70)
        print("This calculator determines molecular weight from vapor pressure lowering.")
        print("Formula: P_solution = P_pure × (1 - X_solute) [Raoult's Law]")
        print("Note: This method assumes ideal solution behavior.")
        print()
        
        try:
            # Get pressure data
            P_pure = float(input("Enter vapor pressure of pure solvent: "))
            P_solution = float(input("Enter vapor pressure of solution: "))
            
            if P_solution >= P_pure:
                print("⚠️  Warning: Solution vapor pressure should be lower than pure solvent.")
            
            # Get mass data
            solute_mass = float(input("Enter mass of solute (g): "))
            solvent = input("Enter solvent name or formula (e.g., water, H2O): ").strip()
            solvent_mass = float(input("Enter mass of solvent (g): "))
            
            # Note about electrolytes
            print("\n📝 Note: Vapor pressure lowering calculations assume non-volatile, non-electrolyte solutes.")
            
            # Solve the problem
            result = solve_vapor_pressure_problem(
                P_pure=P_pure,
                P_solution=P_solution,
                solute_mass=solute_mass,
                solvent=solvent,
                solvent_mass=solvent_mass
            )
            
            self._display_result_summary(result, "vapor pressure lowering")
                
        except ValueError as e:
            print(f"❌ Input Error: Please enter valid numerical values. ({str(e)})")
        except Exception as e:
            print(f"❌ Error: {str(e)}")
    
    def _display_constants(self):
        """Display available constants for different solvents."""
        print("\n" + "="*80)
        print("COLLIGATIVE PROPERTY CONSTANTS DATABASE")
        print("="*80)
        
        # Display Freezing Point Depression Constants
        print("\n🧊 FREEZING POINT DEPRESSION CONSTANTS (Kf)")
        print("-" * 60)
        print(f"{'Solvent':<20} {'Kf (°C/m)':<15} {'Normal F.P. (°C)':<20}")
        print("-" * 60)
        
        # Add normal freezing points for context
        normal_fp = {
            "water": 0.0,
            "benzene": 5.5,
            "cyclohexane": 6.6,
            "camphor": 179.8,
            "acetic_acid": 16.6,
            "naphthalene": 80.2
        }
        
        for solvent, constant in FREEZING_POINT_CONSTANTS.items():
            fp = normal_fp.get(solvent, "N/A")
            fp_str = f"{fp}" if isinstance(fp, str) else f"{fp:.1f}"
            print(f"{solvent.title():<20} {constant:<15.3f} {fp_str:<20}")
        
        # Display Boiling Point Elevation Constants
        print("\n🔥 BOILING POINT ELEVATION CONSTANTS (Kb)")
        print("-" * 60)
        print(f"{'Solvent':<20} {'Kb (°C/m)':<15} {'Normal B.P. (°C)':<20}")
        print("-" * 60)
        
        # Add normal boiling points for context
        normal_bp = {
            "water": 100.0,
            "benzene": 80.1,
            "chloroform": 61.2,
            "ethanol": 78.4,
            "acetic_acid": 118.1
        }
        
        for solvent, constant in BOILING_POINT_CONSTANTS.items():
            bp = normal_bp.get(solvent, "N/A")
            bp_str = f"{bp}" if isinstance(bp, str) else f"{bp:.1f}"
            print(f"{solvent.title():<20} {constant:<15.3f} {bp_str:<20}")
        
        print(f"\n📚 Usage Notes:")
        print(f"• Kf and Kb values are specific to each solvent")
        print(f"• These constants are used in the equations: ΔTf = Kf×m×i and ΔTb = Kb×m×i")
        print(f"• m = molality (mol solute/kg solvent), i = van't Hoff factor")
        print(f"• For precise work, always verify constants from reliable sources")
        
# Add these imports to the top of your ui/colligative_properties_ui.py file:
from chemistry_solver.colligative_properties import (
    solve_freezing_point_problem,
    solve_boiling_point_problem,
    solve_osmotic_pressure_problem,
    solve_vapor_pressure_problem,
    solve_molecular_weight_multiple_choice,  # NEW
    calculate_molecular_weight_from_freezing_point,  # NEW
    calculate_molecular_weight_from_boiling_point,  # NEW
    FREEZING_POINT_CONSTANTS,
    BOILING_POINT_CONSTANTS
)

# MODIFY the _display_menu method to include the new options:
def _display_menu(self):
    """Display the colligative properties module menu."""
    menu = """
    [1] Freezing Point Depression (Calculate MW from temperatures)
    [2] Boiling Point Elevation (Calculate MW from temperatures)
    [3] Osmotic Pressure (Calculate MW)
    [4] Vapor Pressure Lowering (Calculate MW)
    [5] Molecular Weight from Freezing Point Data
    [6] Molecular Weight from Boiling Point Data
    [7] Multiple Choice Problem Solver
    [8] Display Constants
    [0] Return to main menu
    """
    print(menu)

# MODIFY the run method to handle the new choices:
def run(self):
    """Run the colligative properties UI."""
    display_title(self.title)
    
    while True:
        self._display_menu()
        choice = input("\nEnter choice (0-8): ").strip()
        
        if choice == "0":
            # Return to main menu
            return
        elif choice == "1":
            self._handle_freezing_point_depression()
        elif choice == "2":
            self._handle_boiling_point_elevation()
        elif choice == "3":
            self._handle_osmotic_pressure()
        elif choice == "4":
            self._handle_vapor_pressure_lowering()
        elif choice == "5":
            self._handle_molecular_weight_from_freezing_point()  # NEW
        elif choice == "6":
            self._handle_molecular_weight_from_boiling_point()  # NEW
        elif choice == "7":
            self._handle_multiple_choice_solver()  # NEW
        elif choice == "8":
            self._display_constants()
        else:
            print("Invalid choice. Please try again.")
        
        wait_for_user()

# ADD these new methods to the ColligativePropertiesUI class:

def _handle_molecular_weight_from_freezing_point(self):
    """Handle molecular weight calculation from freezing point depression data."""
    print("\n" + "="*70)
    print("MOLECULAR WEIGHT FROM FREEZING POINT DEPRESSION")
    print("="*70)
    print("Calculate molecular weight when you know the freezing point depression directly.")
    print("Formula: ΔTf = Kf × m × i, where m = moles/kg solvent")
    print()
    
    try:
        # Get freezing point depression directly
        delta_T = float(input("Enter freezing point depression (ΔTf in °C): "))
        
        # Get constant
        print("\nSelect freezing point depression constant:")
        print("1. Use known solvent constant")
        print("2. Enter custom Kf value")
        
        kf_choice = input("Choice (1-2): ").strip()
        
        if kf_choice == "1":
            # Show available solvents
            solvents = list(FREEZING_POINT_CONSTANTS.items())
            print("\nAvailable solvents:")
            for i, (solvent, kf) in enumerate(solvents, 1):
                print(f"  [{i}] {solvent.title()}: Kf = {kf} °C/m")
            
            solvent_choice = int(input(f"Select solvent (1-{len(solvents)}): ")) - 1
            if 0 <= solvent_choice < len(solvents):
                solvent_name, K_f = solvents[solvent_choice]
                print(f"Selected: {solvent_name.title()} (Kf = {K_f} °C/m)")
            else:
                raise ValueError("Invalid solvent selection")
        else:
            K_f = float(input("Enter Kf constant (°C/m): "))
            solvent_name = input("Enter solvent name for reference: ").strip()
        
        # Get mass data
        solute_mass = float(input("Enter mass of solute (g): "))
        solvent_mass = float(input("Enter mass of solvent (g): "))
        
        # Get ionization factor
        ionization_factor = self._get_ionization_factor()
        
        # Ask if they want to compare with answer choices
        has_choices = input("Do you have multiple choice answers to compare? (y/n): ").lower().startswith('y')
        answer_choices = None
        
        if has_choices:
            choices_input = input("Enter answer choices separated by commas (e.g., 0.04, 0.06, 0.08): ")
            answer_choices = [float(x.strip()) for x in choices_input.split(',')]
        
        # Calculate result
        result = calculate_molecular_weight_from_freezing_point(
            delta_T=delta_T,
            K_f=K_f,
            solute_mass=solute_mass,
            solvent_mass=solvent_mass,
            ionization_factor=ionization_factor,
            answer_choices=answer_choices
        )
        
        self._display_result_summary(result, "molecular weight from freezing point")
        
    except ValueError as e:
        print(f"❌ Input Error: Please enter valid numerical values. ({str(e)})")
    except Exception as e:
        print(f"❌ Error: {str(e)}")

def _handle_molecular_weight_from_boiling_point(self):
    """Handle molecular weight calculation from boiling point elevation data."""
    print("\n" + "="*70)
    print("MOLECULAR WEIGHT FROM BOILING POINT ELEVATION")
    print("="*70)
    print("Calculate molecular weight when you know the boiling point elevation directly.")
    print("Formula: ΔTb = Kb × m × i, where m = moles/kg solvent")
    print()
    
    try:
        # Get boiling point elevation directly
        delta_T = float(input("Enter boiling point elevation (ΔTb in °C): "))
        
        # Get constant
        print("\nSelect boiling point elevation constant:")
        print("1. Use known solvent constant")
        print("2. Enter custom Kb value")
        
        kb_choice = input("Choice (1-2): ").strip()
        
        if kb_choice == "1":
            # Show available solvents
            solvents = list(BOILING_POINT_CONSTANTS.items())
            print("\nAvailable solvents:")
            for i, (solvent, kb) in enumerate(solvents, 1):
                print(f"  [{i}] {solvent.title()}: Kb = {kb} °C/m")
            
            solvent_choice = int(input(f"Select solvent (1-{len(solvents)}): ")) - 1
            if 0 <= solvent_choice < len(solvents):
                solvent_name, K_b = solvents[solvent_choice]
                print(f"Selected: {solvent_name.title()} (Kb = {K_b} °C/m)")
            else:
                raise ValueError("Invalid solvent selection")
        else:
            K_b = float(input("Enter Kb constant (°C/m): "))
            solvent_name = input("Enter solvent name for reference: ").strip()
        
        # Get mass data
        solute_mass = float(input("Enter mass of solute (g): "))
        solvent_mass = float(input("Enter mass of solvent (g): "))
        
        # Get ionization factor
        ionization_factor = self._get_ionization_factor()
        
        # Ask if they want to compare with answer choices
        has_choices = input("Do you have multiple choice answers to compare? (y/n): ").lower().startswith('y')
        answer_choices = None
        
        if has_choices:
            choices_input = input("Enter answer choices separated by commas (e.g., 0.04, 0.06, 0.08): ")
            answer_choices = [float(x.strip()) for x in choices_input.split(',')]
        
        # Calculate result
        result = calculate_molecular_weight_from_boiling_point(
            delta_T=delta_T,
            K_b=K_b,
            solute_mass=solute_mass,
            solvent_mass=solvent_mass,
            ionization_factor=ionization_factor,
            answer_choices=answer_choices
        )
        
        self._display_result_summary(result, "molecular weight from boiling point")
        
    except ValueError as e:
        print(f"❌ Input Error: Please enter valid numerical values. ({str(e)})")
    except Exception as e:
        print(f"❌ Error: {str(e)}")

    def _handle_multiple_choice_solver(self):
        """Handle structured multiple choice problem solving."""
        print("\n" + "="*70)
        print("MULTIPLE CHOICE PROBLEM SOLVER")
        print("="*70)
        print("Solve colligative properties problems with multiple choice answers.")
        print()
        
        try:
            # Select method
            print("Select calculation method:")
            print("1. Freezing Point Depression")
            print("2. Boiling Point Elevation")
            
            method_choice = input("Choice (1-2): ").strip()
            
            if method_choice == "1":
                method = "freezing_point"
                constant_dict = FREEZING_POINT_CONSTANTS
                constant_name = "Kf"
                temp_change_name = "freezing point depression"
            elif method_choice == "2":
                method = "boiling_point"
                constant_dict = BOILING_POINT_CONSTANTS
                constant_name = "Kb"
                temp_change_name = "boiling point elevation"
            else:
                print("Invalid choice.")
                return
            
            # Get temperature change
            delta_T = float(input(f"Enter {temp_change_name} (°C): "))
            
            # Get constant
            print(f"\nSelect {constant_name} constant:")
            print("1. Use known solvent")
            print("2. Enter custom value")
            
            const_choice = input("Choice (1-2): ").strip()
            
            if const_choice == "1":
                solvents = list(constant_dict.items())
                print(f"\nAvailable solvents:")
                for i, (solvent, const) in enumerate(solvents, 1):
                    print(f"  [{i}] {solvent.title()}: {constant_name} = {const} °C/m")
                
                solvent_choice = int(input(f"Select solvent (1-{len(solvents)}): ")) - 1
                if 0 <= solvent_choice < len(solvents):
                    _, constant = solvents[solvent_choice]
                else:
                    raise ValueError("Invalid solvent selection")
            else:
                constant = float(input(f"Enter {constant_name} constant (°C/m): "))
            
            # Get mass data
            solute_mass = float(input("Enter mass of solute (g): "))
            solvent_mass = float(input("Enter mass of solvent (g): "))
            
            # Get ionization factor
            ionization_factor = self._get_ionization_factor()
            
            # Get answer choices
            choices_input = input("Enter answer choices separated by commas: ")
            answer_choices = [float(x.strip()) for x in choices_input.split(',')]
            
            # Prepare problem data
            problem_data = {
                'method': method,
                'delta_T': delta_T,
                'constant': constant,
                'solute_mass': solute_mass,
                'solvent_mass': solvent_mass,
                'ionization_factor': ionization_factor,
                'answer_choices': answer_choices
            }
            
            # Solve the problem
            result = solve_molecular_weight_multiple_choice(problem_data)
            
            if result['success']:
                print(f"\n{'='*70}")
                print(f"MULTIPLE CHOICE PROBLEM SOLUTION")
                print(f"{'='*70}")
                
                # Display the calculation steps
                display_steps(result['steps'])
                
                # Highlight the answer
                if 'closest_answer' in result:
                    print(f"\n🎯 ANSWER: {result['closest_answer']}")
                    if 'answer_difference' in result:
                        print(f"   Difference from calculated: {result['answer_difference']:.4f}")
                
                # Show all answer choices for reference
                print(f"\nAnswer choices were: {answer_choices}")
            else:
                print(f"❌ Error: {result['error']}")
            
        except ValueError as e:
            print(f"❌ Input Error: Please enter valid numerical values. ({str(e)})")
        except Exception as e:
            print(f"❌ Error: {str(e)}")

    # MODIFY the _display_result_summary method to handle new result types:
    def _display_result_summary(self, result, calculation_type):
        """Display a summary of calculation results."""
        print(f"\n{'='*60}")
        print(f"CALCULATION SUMMARY - {calculation_type.upper()}")
        print(f"{'='*60}")
        
        if result['success']:
            # Display key results
            if 'molecular_weight' in result:
                print(f"Molecular Weight: {result['molecular_weight']:.2f} g/mol")
            
            if 'closest_answer' in result:
                print(f"Closest Answer Choice: {result['closest_answer']}")
                if 'answer_difference' in result:
                    print(f"Difference: {result['answer_difference']:.4f}")
            
            # Display specific properties based on calculation type
            if 'delta_T' in result:
                temp_change = "Depression" if "freezing" in calculation_type else "Elevation"
                print(f"Temperature {temp_change} (ΔT): {result['delta_T']:.4f} °C")
            if 'molality' in result:
                print(f"Molality: {result['molality']:.6f} mol/kg")
            if 'moles_solute' in result:
                print(f"Moles of Solute: {result['moles_solute']:.6f} mol")
            if 'mole_fraction_solute' in result:
                print(f"Mole Fraction of Solute: {result['mole_fraction_solute']:.6f}")
            if 'delta_P' in result:
                print(f"Vapor Pressure Lowering (ΔP): {result['delta_P']:.4f}")
            
            print(f"\n{'DETAILED CALCULATION STEPS'}")
            print(f"{'-'*60}")
            display_steps(result['steps'])
        else:
            print(f"❌ ERROR: {result['error']}")