"""
Terminal User Interface for Stoichiometry Calculations
"""
from ui.terminal_ui import display_title, display_results_header, display_steps, wait_for_user
from chemistry_solver.stoichiometry import (
    solve_stoichiometry_problem,
    solve_multireactant_problem,
    solve_gas_stoichiometry_problem
)


class StoichiometryUI:
    """UI class for stoichiometry calculations."""
    
    def __init__(self):
        self.title = "STOICHIOMETRY CALCULATOR"
    
    def run(self):
        """Run the stoichiometry UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-4): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_basic_stoichiometry()
            elif choice == "2":
                self._handle_limiting_reactant()
            elif choice == "3":
                self._handle_gas_stoichiometry()
            elif choice == "4":
                self._handle_multireactant_fixed_masses()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the stoichiometry module menu."""
        menu = """
        [1] Basic Stoichiometry (Mass-to-Mass)
        [2] Limiting Reactant Analysis
        [3] Gas Stoichiometry
        [4] Fixed Mass Multireactant Problem
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_basic_stoichiometry(self):
        """Handle basic stoichiometry calculations."""
        print("\n===== BASIC STOICHIOMETRY CALCULATOR =====")
        
        try:
            print("\nEnter the chemical equation (use -> or → for the reaction arrow):")
            equation = input("Equation: ")
            
            given_compound = input("\nEnter the formula of the compound you're starting with: ")
            given_mass = float(input("Enter the mass (in grams) of this compound: "))
            
            target_compound = input("\nEnter the formula of the compound you want to find: ")
            
            result = solve_stoichiometry_problem(equation, given_compound, given_mass, target_compound)
            
            display_results_header()
            print(f"Balanced Equation: {result['balanced_equation']}")
            print(f"Target Mass: {result['target_mass']:.4f} g")
            print(f"Target Moles: {result['target_moles']:.6f} mol")
            
            display_steps(result['steps'])
                
        except ValueError as e:
            print(f"Error: {str(e)}")
        except Exception as e:
            print(f"An unexpected error occurred: {str(e)}")
    
    def _handle_limiting_reactant(self):
        """Handle limiting reactant calculations."""
        print("\n===== LIMITING REACTANT ANALYSIS =====")
        
        try:
            print("\nEnter the chemical equation (use -> or → for the reaction arrow):")
            equation = input("Equation: ")
            
            reactant_data = {}
            print("\nEnter reactant information (leave empty to finish):")
            
            while True:
                reactant = input("\nEnter reactant formula (or press Enter to finish): ")
                if not reactant:
                    break
                    
                mass = float(input(f"Enter mass of {reactant} in grams: "))
                reactant_data[reactant] = mass
            
            if not reactant_data:
                print("Error: No reactants provided.")
                return
            
            target_compound = input("\nEnter the formula of the product you want to find: ")
            
            result = solve_multireactant_problem(equation, reactant_data, target_compound)
            
            display_results_header()
            print(f"Balanced Equation: {result['balanced_equation']}")
            print(f"Limiting Reactant: {result['limiting_reactant']}")
            print(f"Maximum Product Mass: {result['target_mass']:.4f} g")
            
            display_steps(result['steps'])
                
        except ValueError as e:
            print(f"Error: {str(e)}")
        except Exception as e:
            print(f"An unexpected error occurred: {str(e)}")
    
    def _handle_multireactant_fixed_masses(self):
        """Handle a problem with specific fixed masses of multiple reactants."""
        print("\n===== FIXED MASS MULTIREACTANT PROBLEM =====")
        
        try:
            print("\nEnter the chemical equation (use -> or → for the reaction arrow):")
            equation = input("Equation: ")
            
            # Get masses for reactants
            reactant_data = {}
            
            # First reactant
            reactant1 = input("\nEnter first reactant formula: ")
            mass1 = float(input(f"Enter mass of {reactant1} in grams: "))
            reactant_data[reactant1] = mass1
            
            # Second reactant
            reactant2 = input("\nEnter second reactant formula: ")
            mass2 = float(input(f"Enter mass of {reactant2} in grams: "))
            reactant_data[reactant2] = mass2
            
            # Ask for additional reactants
            add_more = input("\nDo you want to add more reactants? (y/n): ").lower()
            while add_more == 'y':
                reactant = input("Enter reactant formula: ")
                mass = float(input(f"Enter mass of {reactant} in grams: "))
                reactant_data[reactant] = mass
                add_more = input("Do you want to add more reactants? (y/n): ").lower()
            
            # Get target product
            target_compound = input("\nEnter the formula of the product you want to find: ")
            
            result = solve_multireactant_problem(equation, reactant_data, target_compound)
            
            display_results_header()
            print(f"Balanced Equation: {result['balanced_equation']}")
            print(f"Limiting Reactant: {result['limiting_reactant']}")
            print(f"Maximum Product Mass: {result['target_mass']:.4f} g")
            
            display_steps(result['steps'])
                
        except ValueError as e:
            print(f"Error: {str(e)}")
        except Exception as e:
            print(f"An unexpected error occurred: {str(e)}")
    
    def _handle_gas_stoichiometry(self):
        """Handle gas stoichiometry calculations."""
        print("\n===== GAS STOICHIOMETRY CALCULATOR =====")
        
        try:
            print("\nEnter the chemical equation (use -> or → for the reaction arrow):")
            equation = input("Equation: ")
            
            given_compound = input("\nEnter the formula of the compound you're starting with: ")
            given_mass = float(input("Enter the mass (in grams) of this compound: "))
            
            target_gas = input("\nEnter the formula of the gas product: ")
            
            temperature = float(input("\nEnter the temperature (in °C, default 0): ") or "0")
            pressure = float(input("Enter the pressure (in atm, default 1.0): ") or "1.0")
            
            result = solve_gas_stoichiometry_problem(
                equation, given_compound, given_mass, target_gas, temperature, pressure
            )
            
            display_results_header()
            print(f"Balanced Equation: {result['balanced_equation']}")
            print(f"Target Gas: {target_gas}")
            print(f"Gas Mass: {result['target_mass']:.4f} g")
            print(f"Gas Moles: {result['target_moles']:.6f} mol")
            print(f"Gas Volume: {result['gas_volume']:.4f} L (at {temperature}°C and {pressure} atm)")
            
            display_steps(result['steps'])
                
        except ValueError as e:
            print(f"Error: {str(e)}")
        except Exception as e:
            print(f"An unexpected error occurred: {str(e)}")