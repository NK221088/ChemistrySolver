"""
Terminal User Interface for Colligative Properties Calculations
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user, display_steps
from chemistry_solver.colligative_properties import (
    solve_freezing_point_problem,
    solve_boiling_point_problem,
    solve_osmotic_pressure_problem,
    solve_vapor_pressure_problem,
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
            choice = input("\nEnter choice (0-5): ").strip()
            
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
                self._display_constants()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the colligative properties module menu."""
        menu = """
        [1] Freezing Point Depression
        [2] Boiling Point Elevation
        [3] Osmotic Pressure
        [4] Vapor Pressure Lowering
        [5] Display Constants
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