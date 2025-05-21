"""
Terminal User Interface for Molar Mass Calculations
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user
from chemistry_solver.molar_mass import calculate_molar_mass
from chemistry_solver.name_to_formula import get_formula_from_name

class MolarMassUI:
    """UI class for molar mass calculations."""
    
    def __init__(self):
        self.title = "MOLAR MASS CALCULATOR"
    
    def run(self):
        """Run the molar mass UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-3): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_formula_calculation()
            elif choice == "2":
                self._handle_name_to_formula()
            elif choice == "3":
                self._handle_full_calculation()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the molar mass module menu."""
        menu = """
        [1] Calculate molar mass from formula
        [2] Convert chemical name to formula
        [3] Calculate molar mass from name or formula
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_formula_calculation(self):
        """Handle calculation of molar mass from formula."""
        print("\n===== MOLAR MASS CALCULATOR =====")
        
        try:
            formula = input("Enter chemical formula: ")
            
            result = calculate_molar_mass(formula)
            
            display_results_header()
            if result['success']:
                print(f"Formula: {result['formula']}")
                print(f"Hill Notation: {result['hill_notation']}")
                print(f"Molar Mass: {result['molar_mass']:.4f} g/mol")
                print(f"Monoisotopic Mass: {result['monoisotopic_mass']:.4f} u")
                print("\nElement Composition:")
                print("-" * 40)
                print(f"{'Element':<10} {'Count':<8} {'Atomic Mass':<15} {'Contribution':<15}")
                print("-" * 40)
                for e in result['composition']:
                    print(f"{e['element']:<10} {e['count']:<8} {e['atomic_mass']:.4f} u      {e['contribution']:.4f} g/mol")
            else:
                print(f"Error: {result['error']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_name_to_formula(self):
        """Handle conversion of chemical name to formula."""
        print("\n===== CHEMICAL NAME TO FORMULA CONVERTER =====")
        
        try:
            name = input("Enter chemical name: ")
            
            result = get_formula_from_name(name)
            
            display_results_header()
            if result['success']:
                print(f"Name: {result['name']}")
                print(f"Formula: {result['formula']}")
                print(f"IUPAC Name: {result['iupac_name']}")
                print(f"Molecular Weight: {result['weight']:.4f} g/mol")
            else:
                print(f"Error: {result['error']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_full_calculation(self):
        """Handle full workflow - name to formula to molar mass."""
        print("\n===== MOLAR MASS FROM NAME OR FORMULA =====")
        
        try:
            input_type = input("Do you want to enter a chemical name or formula? (name/formula): ").lower()
            
            formula = None
            if input_type == 'name':
                print("\nConverting chemical name to formula...")
                name = input("Enter chemical name: ")
                name_result = get_formula_from_name(name)
                
                if name_result['success']:
                    print(f"\nName: {name_result['name']}")
                    print(f"Formula: {name_result['formula']}")
                    print(f"IUPAC Name: {name_result['iupac_name']}")
                    print(f"Molecular Weight: {name_result['weight']:.4f} g/mol")
                    
                    use_formula = input("\nUse this formula for molar mass calculation? (y/n): ").lower()
                    if use_formula == 'y':
                        formula = name_result['formula']
                    else:
                        formula = input("\nFallback: Please enter chemical formula directly: ")
                else:
                    print(f"\nError: {name_result['error']}")
                    formula = input("\nFallback: Please enter chemical formula directly: ")
            else:
                formula = input("Enter chemical formula: ")
            
            if formula:
                result = calculate_molar_mass(formula)
                
                display_results_header()
                if result['success']:
                    print(f"Formula: {result['formula']}")
                    print(f"Hill Notation: {result['hill_notation']}")
                    print(f"Molar Mass: {result['molar_mass']:.4f} g/mol")
                    print(f"Monoisotopic Mass: {result['monoisotopic_mass']:.4f} u")
                    print("\nElement Composition:")
                    print("-" * 40)
                    print(f"{'Element':<10} {'Count':<8} {'Atomic Mass':<15} {'Contribution':<15}")
                    print("-" * 40)
                    for e in result['composition']:
                        print(f"{e['element']:<10} {e['count']:<8} {e['atomic_mass']:.4f} u      {e['contribution']:.4f} g/mol")
                else:
                    print(f"Error: {result['error']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")