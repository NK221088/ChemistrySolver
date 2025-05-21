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
    
    def _handle_freezing_point_depression(self):
        """Handle freezing point depression calculations."""
        print("\n===== FREEZING POINT DEPRESSION CALCULATOR =====")
        print("\nThis calculator determines molecular weight from freezing point depression.")
        
        try:
            # Get input parameters
            T_pure = float(input("Enter freezing point of pure solvent (°C): "))
            T_solution = float(input("Enter freezing point of solution (°C): "))
            solvent = input("Enter solvent name or formula (e.g., water, H2O): ")
            
            # Ask if user wants to provide K_f or use a lookup
            use_lookup = input("Use known freezing point constant? (y/n): ").lower() == 'y'
            K_f = None
            if use_lookup:
                print("\nAvailable solvents with known constants:")
                for i, (solvent_name, constant) in enumerate(FREEZING_POINT_CONSTANTS.items(), 1):
                    print(f"  [{i}] {solvent_name} (K_f = {constant} °C/m)")
                
                choice = input("\nSelect a solvent (number) or press Enter to use input name: ")
                if choice.strip():
                    idx = int(choice) - 1
                    if 0 <= idx < len(FREEZING_POINT_CONSTANTS):
                        solvent = list(FREEZING_POINT_CONSTANTS.keys())[idx]
                        K_f = list(FREEZING_POINT_CONSTANTS.values())[idx]
            
            if K_f is None:
                K_f_input = input(f"Enter the freezing point constant K_f for {solvent} (°C/m) [leave blank to use lookup]: ")
                if K_f_input.strip():
                    K_f = float(K_f_input)
            
            solute_mass = float(input("Enter mass of solute (g): "))
            solvent_mass = float(input("Enter mass of solvent (g): "))
            
            # Ask for van't Hoff factor if needed
            ionization = input("Is the solute an electrolyte? (y/n): ").lower() == 'y'
            ionization_factor = 1
            if ionization:
                ionization_factor = float(input("Enter van't Hoff factor (i): "))
            
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
            
            display_results_header()
            if result['success']:
                print(f"Freezing Point Depression (ΔT): {result['delta_T']:.4f} °C")
                print(f"Molality: {result['molality']:.6f} mol/kg")
                print(f"Moles of Solute: {result['moles_solute']:.6f} mol")
                print(f"Molecular Weight: {result['molecular_weight']:.2f} g/mol")
                
                # Display calculation steps
                display_steps(result['steps'])
            else:
                print(f"Error: {result['error']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_boiling_point_elevation(self):
        """Handle boiling point elevation calculations."""
        print("\n===== BOILING POINT ELEVATION CALCULATOR =====")
        print("\nThis calculator determines molecular weight from boiling point elevation.")
        
        try:
            # Get input parameters
            T_pure = float(input("Enter boiling point of pure solvent (°C): "))
            T_solution = float(input("Enter boiling point of solution (°C): "))
            solvent = input("Enter solvent name or formula (e.g., water, H2O): ")
            
            # Ask if user wants to provide K_b or use a lookup
            use_lookup = input("Use known boiling point constant? (y/n): ").lower() == 'y'
            K_b = None
            if use_lookup:
                print("\nAvailable solvents with known constants:")
                for i, (solvent_name, constant) in enumerate(BOILING_POINT_CONSTANTS.items(), 1):
                    print(f"  [{i}] {solvent_name} (K_b = {constant} °C/m)")
                
                choice = input("\nSelect a solvent (number) or press Enter to use input name: ")
                if choice.strip():
                    idx = int(choice) - 1
                    if 0 <= idx < len(BOILING_POINT_CONSTANTS):
                        solvent = list(BOILING_POINT_CONSTANTS.keys())[idx]
                        K_b = list(BOILING_POINT_CONSTANTS.values())[idx]
            
            if K_b is None:
                K_b_input = input(f"Enter the boiling point constant K_b for {solvent} (°C/m) [leave blank to use lookup]: ")
                if K_b_input.strip():
                    K_b = float(K_b_input)
            
            solute_mass = float(input("Enter mass of solute (g): "))
            solvent_mass = float(input("Enter mass of solvent (g): "))
            
            # Ask for van't Hoff factor if needed
            ionization = input("Is the solute an electrolyte? (y/n): ").lower() == 'y'
            ionization_factor = 1
            if ionization:
                ionization_factor = float(input("Enter van't Hoff factor (i): "))
            
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
            
            display_results_header()
            if result['success']:
                print(f"Boiling Point Elevation (ΔT): {result['delta_T']:.4f} °C")
                print(f"Molality: {result['molality']:.6f} mol/kg")
                print(f"Moles of Solute: {result['moles_solute']:.6f} mol")
                print(f"Molecular Weight: {result['molecular_weight']:.2f} g/mol")
                
                # Display calculation steps
                display_steps(result['steps'])
            else:
                print(f"Error: {result['error']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_osmotic_pressure(self):
        """Handle osmotic pressure calculations."""
        print("\n===== OSMOTIC PRESSURE CALCULATOR =====")
        print("\nThis calculator determines molecular weight from osmotic pressure.")
        
        try:
            # Get input parameters
            osmotic_pressure_atm = float(input("Enter osmotic pressure (atm): "))
            temperature_c = float(input("Enter temperature (°C): "))
            solution_volume_L = float(input("Enter solution volume (L): "))
            solute_mass = float(input("Enter mass of solute (g): "))
            
            # Ask for van't Hoff factor if needed
            ionization = input("Is the solute an electrolyte? (y/n): ").lower() == 'y'
            ionization_factor = 1
            if ionization:
                ionization_factor = float(input("Enter van't Hoff factor (i): "))
            
            # Solve the problem
            result = solve_osmotic_pressure_problem(
                osmotic_pressure_atm=osmotic_pressure_atm,
                temperature_c=temperature_c,
                solution_volume_L=solution_volume_L,
                solute_mass=solute_mass,
                ionization_factor=ionization_factor
            )
            
            display_results_header()
            if result['success']:
                print(f"Moles of Solute: {result['moles_solute']:.6f} mol")
                print(f"Molecular Weight: {result['molecular_weight']:.2f} g/mol")
                
                # Display calculation steps
                display_steps(result['steps'])
            else:
                print(f"Error: {result['error']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_vapor_pressure_lowering(self):
        """Handle vapor pressure lowering calculations."""
        print("\n===== VAPOR PRESSURE LOWERING CALCULATOR =====")
        print("\nThis calculator determines molecular weight from vapor pressure lowering.")
        
        try:
            # Get input parameters
            P_pure = float(input("Enter vapor pressure of pure solvent: "))
            P_solution = float(input("Enter vapor pressure of solution: "))
            solute_mass = float(input("Enter mass of solute (g): "))
            solvent = input("Enter solvent name or formula (e.g., water, H2O): ")
            solvent_mass = float(input("Enter mass of solvent (g): "))
            
            # Solve the problem
            result = solve_vapor_pressure_problem(
                P_pure=P_pure,
                P_solution=P_solution,
                solute_mass=solute_mass,
                solvent=solvent,
                solvent_mass=solvent_mass
            )
            
            display_results_header()
            if result['success']:
                print(f"Vapor Pressure Lowering (ΔP): {result['delta_P']:.4f}")
                print(f"Mole Fraction of Solute: {result['mole_fraction_solute']:.6f}")
                print(f"Moles of Solvent: {result['moles_solvent']:.6f} mol")
                print(f"Moles of Solute: {result['moles_solute']:.6f} mol")
                print(f"Molecular Weight: {result['molecular_weight']:.2f} g/mol")
                
                # Display calculation steps
                display_steps(result['steps'])
            else:
                print(f"Error: {result['error']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _display_constants(self):
        """Display available constants for different solvents."""
        print("\n===== COLLIGATIVE PROPERTY CONSTANTS =====")
        
        # Display Freezing Point Depression Constants
        print("\nFreezing Point Depression Constants (K_f):")
        print("-" * 50)
        print(f"{'Solvent':<20} {'K_f (°C/m)':<15}")
        print("-" * 50)
        for solvent, constant in FREEZING_POINT_CONSTANTS.items():
            print(f"{solvent:<20} {constant:<15.3f}")
        
        # Display Boiling Point Elevation Constants
        print("\nBoiling Point Elevation Constants (K_b):")
        print("-" * 50)
        print(f"{'Solvent':<20} {'K_b (°C/m)':<15}")
        print("-" * 50)
        for solvent, constant in BOILING_POINT_CONSTANTS.items():
            print(f"{solvent:<20} {constant:<15.3f}")