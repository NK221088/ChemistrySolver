"""
Enhanced Terminal User Interface for Molar Mass and Chemistry Calculations
Properly integrated with the enhanced molar mass calculator module
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user
from chemistry_solver.molar_mass import (
    calculate_molar_mass, 
    count_atoms_in_formula,
    calculate_atoms_from_mass,
    calculate_mass_from_atoms,
    solve_isopropanol_problem,
    AVOGADROS_NUMBER
)

# Try to import name_to_formula, but don't fail if it doesn't exist
try:
    from chemistry_solver.name_to_formula import get_formula_from_name
    NAME_TO_FORMULA_AVAILABLE = True
except ImportError:
    NAME_TO_FORMULA_AVAILABLE = False

class MolarMassUI:
    """Enhanced UI class for comprehensive molar mass and chemistry calculations."""
    
    def __init__(self):
        self.title = "MOLAR MASS & CHEMISTRY CALCULATOR"
    
    def run(self):
        """Run the molar mass and chemistry UI."""
        display_title(self.title)
        
        while True:
            self._display_main_menu()
            choice = input("\nEnter choice (0-7): ").strip()
            
            if choice == "0":
                return
            elif choice == "1":
                self._handle_molar_mass()
            elif choice == "2":
                self._handle_atom_counting()
            elif choice == "3":
                self._handle_mass_to_atoms()
            elif choice == "4":
                self._handle_atoms_to_mass()
            elif choice == "5":
                self._handle_problem_solver()
            elif choice == "6":
                self._handle_isopropanol_example()
            elif choice == "7":
                if NAME_TO_FORMULA_AVAILABLE:
                    self._handle_name_to_formula()
                else:
                    print("Name to formula converter not available.")
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_main_menu(self):
        """Display the main menu."""
        menu = """
Available Functions:
  [1] Calculate molar mass from formula
  [2] Count atoms in formula
  [3] Calculate number of atoms from mass
  [4] Calculate mass from number of atoms
  [5] Chemistry problem solver (custom problems)
  [6] Isopropanol example problem
  [7] Convert chemical name to formula"""
        
        if not NAME_TO_FORMULA_AVAILABLE:
            menu += " (unavailable)"
        
        menu += "\n  [0] Return to main menu"
        
        print(menu)
    
    def _handle_molar_mass(self):
        """Handle molar mass calculations."""
        print("\n===== MOLAR MASS CALCULATOR =====")
        
        try:
            formula = input("Enter chemical formula (e.g., H2O, C6H12O6, CaCl2): ").strip()
            if not formula:
                print("No formula entered.")
                return
                
            result = calculate_molar_mass(formula)
            
            display_results_header()
            if result['success']:
                print(f"Formula: {result['formula']}")
                print(f"Hill Notation: {result['hill_notation']}")
                print(f"Molar Mass: {result['molar_mass']:.4f} g/mol")
                print(f"Monoisotopic Mass: {result['monoisotopic_mass']:.4f} u")
                print("\nElement Composition:")
                print("-" * 65)
                print(f"{'Element':<10} {'Count':<8} {'Atomic Mass':<15} {'Contribution':<15}")
                print("-" * 65)
                for e in result['composition']:
                    print(f"{e['element']:<10} {e['count']:<8} {e['atomic_mass']:.4f} u      {e['contribution']:.4f} g/mol")
            else:
                print(f"Error: {result['error']}")
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_atom_counting(self):
        """Handle atom counting in formulas."""
        print("\n===== ATOM COUNTER =====")
        
        try:
            formula = input("Enter chemical formula (e.g., C3H8O): ").strip()
            if not formula:
                print("No formula entered.")
                return
                
            element = input("Enter specific element (or press Enter for all atoms): ").strip()
            
            if not element:
                element = None
            
            result = count_atoms_in_formula(formula, element)
            
            display_results_header()
            if result['success']:
                if element:
                    print(f"Formula: {formula}")
                    print(f"Element: {result['element']}")
                    print(f"Number of {result['element']} atoms per molecule: {result['count']}")
                else:
                    print(f"Formula: {formula}")
                    print(f"Total atoms per molecule: {result['total_atoms']}")
                    print("\nAtom breakdown:")
                    for elem, count in result['atom_breakdown'].items():
                        print(f"  {elem}: {count}")
            else:
                print(f"Error: {result['error']}")
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_mass_to_atoms(self):
        """Handle conversion from mass to number of atoms."""
        print("\n===== MASS TO ATOMS CALCULATOR =====")
        
        try:
            formula = input("Enter chemical formula: ").strip()
            if not formula:
                print("No formula entered.")
                return
                
            mass_input = input("Enter mass in grams: ").strip()
            if not mass_input:
                print("No mass entered.")
                return
                
            mass = float(mass_input)
            element = input("Enter specific element (or press Enter for total atoms): ").strip()
            
            if not element:
                element = None
            
            result = calculate_atoms_from_mass(formula, mass, element)
            
            display_results_header()
            if result['success']:
                print(f"Formula: {result['formula']}")
                print(f"Mass: {result['mass_grams']} g")
                print(f"Molar Mass: {result['molar_mass']:.4f} g/mol")
                print(f"Moles: {result['moles']:.6f} mol")
                print(f"{result['element'].capitalize()} atoms per molecule: {result['atoms_per_molecule']}")
                print(f"Total {result['element']} atoms: {result['total_atoms']:.3e}")
                print(f"\nAnswer: {result['total_atoms']:.2e} atoms")
            else:
                print(f"Error: {result['error']}")
        except ValueError:
            print("Error: Invalid mass value. Please enter a number.")
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_atoms_to_mass(self):
        """Handle conversion from number of atoms to mass."""
        print("\n===== ATOMS TO MASS CALCULATOR =====")
        
        try:
            formula = input("Enter chemical formula: ").strip()
            if not formula:
                print("No formula entered.")
                return
                
            atoms_input = input("Enter number of atoms (can use scientific notation like 1.5e23): ").strip()
            if not atoms_input:
                print("No atom count entered.")
                return
                
            num_atoms = float(atoms_input)
            element = input("Enter specific element (or press Enter for total atoms): ").strip()
            
            if not element:
                element = None
            
            result = calculate_mass_from_atoms(formula, num_atoms, element)
            
            display_results_header()
            if result['success']:
                print(f"Formula: {result['formula']}")
                print(f"Number of {result['element']} atoms: {result['num_atoms']:.3e}")
                print(f"{result['element'].capitalize()} atoms per molecule: {result['atoms_per_molecule']}")
                print(f"Number of molecules: {result['molecules']:.3e}")
                print(f"Moles: {result['moles']:.6f} mol")
                print(f"Molar Mass: {result['molar_mass']:.4f} g/mol")
                print(f"\nMass: {result['mass_grams']:.4f} g")
            else:
                print(f"Error: {result['error']}")
        except ValueError:
            print("Error: Invalid number format. Please enter a valid number.")
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_problem_solver(self):
        """Handle custom chemistry problems."""
        print("\n===== CHEMISTRY PROBLEM SOLVER =====")
        print("This section helps solve problems involving Avogadro's number")
        print("Example: 'How many hydrogen atoms are in 36.25g of isopropanol?'")
        
        try:
            problem_type = input("\nProblem type:\n[1] Find atoms from mass\n[2] Find mass from atoms\nChoice: ").strip()
            
            if problem_type == "1":
                formula = input("Enter chemical formula: ").strip()
                if not formula:
                    print("No formula entered.")
                    return
                    
                mass_input = input("Enter mass in grams: ").strip()
                if not mass_input:
                    print("No mass entered.")
                    return
                    
                mass = float(mass_input)
                element = input("Enter element to count: ").strip()
                if not element:
                    print("No element entered.")
                    return
                
                result = calculate_atoms_from_mass(formula, mass, element)
                
                display_results_header()
                if result['success']:
                    print(f"PROBLEM: How many {element} atoms are in {mass}g of {formula}?")
                    print(f"\nSOLUTION STEPS:")
                    print(f"1. Formula: {formula}")
                    print(f"2. Molar mass: {result['molar_mass']:.4f} g/mol")
                    print(f"3. Moles = {mass}g ÷ {result['molar_mass']:.4f} g/mol = {result['moles']:.6f} mol")
                    print(f"4. {element} atoms per molecule: {result['atoms_per_molecule']}")
                    print(f"5. Total atoms = {result['moles']:.6f} mol × {AVOGADROS_NUMBER:.3e} × {result['atoms_per_molecule']}")
                    print(f"\nFINAL ANSWER: {result['total_atoms']:.2e} {element} atoms")
                else:
                    print(f"Error: {result['error']}")
            
            elif problem_type == "2":
                formula = input("Enter chemical formula: ").strip()
                if not formula:
                    print("No formula entered.")
                    return
                    
                atoms_input = input("Enter number of atoms: ").strip()
                if not atoms_input:
                    print("No atom count entered.")
                    return
                    
                num_atoms = float(atoms_input)
                element = input("Enter element: ").strip()
                if not element:
                    element = None
                
                result = calculate_mass_from_atoms(formula, num_atoms, element)
                
                display_results_header()
                if result['success']:
                    print(f"PROBLEM: What mass contains {num_atoms:.2e} {element or 'total'} atoms of {formula}?")
                    print(f"\nSOLUTION STEPS:")
                    print(f"1. Formula: {formula}")
                    print(f"2. {result['element'].capitalize()} atoms per molecule: {result['atoms_per_molecule']}")
                    print(f"3. Molecules = {num_atoms:.2e} ÷ {result['atoms_per_molecule']} = {result['molecules']:.3e}")
                    print(f"4. Moles = {result['molecules']:.3e} ÷ {AVOGADROS_NUMBER:.3e} = {result['moles']:.6f} mol")
                    print(f"5. Mass = {result['moles']:.6f} mol × {result['molar_mass']:.4f} g/mol")
                    print(f"\nFINAL ANSWER: {result['mass_grams']:.4f} g")
                else:
                    print(f"Error: {result['error']}")
            else:
                print("Invalid choice.")
                
        except ValueError:
            print("Error: Invalid number format.")
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_isopropanol_example(self):
        """Handle the specific isopropanol example problem."""
        print("\n===== ISOPROPANOL EXAMPLE PROBLEM =====")
        print("Problem: How many hydrogen atoms are in 36.25g of isopropanol?")
        
        try:
            result = solve_isopropanol_problem()
            
            display_results_header()
            if result['success']:
                print(f"Problem: {result['problem']}")
                print(f"\nSOLUTION STEPS:")
                print(f"1. Isopropanol formula: {result['formula']} (C₃H₈O)")
                print(f"2. Molar mass: {result['molar_mass']:.4f} g/mol")
                print(f"3. Moles = {result['mass_grams']}g ÷ {result['molar_mass']:.4f} g/mol = {result['moles']:.6f} mol")
                print(f"4. Hydrogen atoms per molecule: {result['hydrogen_atoms_per_molecule']}")
                print(f"5. Total H atoms = {result['moles']:.6f} mol × {AVOGADROS_NUMBER:.3e} × {result['hydrogen_atoms_per_molecule']}")
                print(f"\nFINAL ANSWER: {result['scientific_notation']} H atoms")
                
                # Show answer choices comparison
                print(f"\nAnswer choices comparison:")
                choices = [
                    ("A", 2.90e12),
                    ("B", 2.40e12), 
                    ("C", 2.90e24),
                    ("D", 5.80e24),
                    ("E", 4.80e12)
                ]
                
                for letter, value in choices:
                    diff = abs(result['total_hydrogen_atoms'] - value)
                    match = " ✓ CORRECT" if diff < value * 0.01 else ""  # Within 1%
                    print(f"  {letter}: {value:.2e}{match}")
            else:
                print(f"Error: {result['error']}")
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_name_to_formula(self):
        """Handle chemical name to formula conversion."""
        print("\n===== CHEMICAL NAME TO FORMULA CONVERTER =====")
        
        try:
            name = input("Enter chemical name: ").strip()
            if not name:
                print("No name entered.")
                return
                
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

# For backwards compatibility with the old class name
EnhancedChemistryUI = MolarMassUI