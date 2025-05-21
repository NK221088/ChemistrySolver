"""
Terminal User Interface for Electron Configuration Problems
"""
import time
from ui.terminal_ui import display_title, display_results_header, display_steps, wait_for_user
from chemistry_solver.electron_config import (
    solve_electron_problem, 
    solve_molecular_orbital_problem,
    get_element_by_symbol,
    get_element_by_number
)

class ElectronConfigUI:
    """UI class for electron configuration problems."""
    
    def __init__(self):
        self.title = "ELECTRON CONFIGURATION PROBLEM SOLVER"
    
    def run(self):
        """Run the electron configuration UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-5): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._solve_unpaired_electrons()
            elif choice == "2":
                self._solve_orbital_count()
            elif choice == "3":
                self._solve_electron_config()
            elif choice == "4":
                self._solve_quantum_numbers()
            elif choice == "5":
                self._solve_molecular_orbitals()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the electron configuration module menu."""
        menu = """
        [1] Find number of unpaired electrons in an atom
        [2] Count electrons in specific orbital types (s, p, d, f)
        [3] Generate electron configuration
        [4] Analyze quantum numbers for specific electron
        [5] Solve molecular orbital problems
        [0] Return to main menu
        """
        print(menu)
    
    def _get_element_input(self):
        """Get element input from user and validate."""
        while True:
            element_input = input("\nEnter element symbol or atomic number: ").strip()
            
            # Check if input is atomic number
            if element_input.isdigit():
                atomic_num = int(element_input)
                result = get_element_by_number(atomic_num)
                if result["success"]:
                    return result["symbol"], atomic_num
                print(f"Error: {result['error']}")
            
            # Check if input is element symbol
            else:
                result = get_element_by_symbol(element_input)
                if result["success"]:
                    return result["symbol"], result["atomic_number"]
                print(f"Error: {result['error']}")
    
    def _solve_unpaired_electrons(self):
        """Solve problem for unpaired electrons."""
        print("\n===== UNPAIRED ELECTRONS ANALYZER =====")
        
        element_symbol, atomic_number = self._get_element_input()
        
        print(f"\nAnalyzing unpaired electrons for {element_symbol} (atomic number {atomic_number})...")
        time.sleep(0.5)  # Brief pause for effect
        
        result = solve_electron_problem("unpaired", element_symbol)
        
        if result["success"]:
            display_results_header()
            print(f"Element: {result['element']} (atomic number {atomic_number})")
            print(f"Electron configuration: {result['configuration']}")
            print(f"Total unpaired electrons: {result['unpaired_electrons']}")
            
            if "diagrams" in result:
                print("\nOrbital diagrams:")
                for diagram in result["diagrams"]:
                    print(diagram)
            
            # Filtering steps - show from step 3 onwards (skip first two)
            steps_to_display = result["steps"][2:] if len(result["steps"]) > 2 else []
            display_steps(steps_to_display)
        else:
            print(f"Error: {result['error']}")
    
    def _solve_orbital_count(self):
        """Solve problem for counting electrons in specific orbital types."""
        print("\n===== ORBITAL ELECTRON COUNTER =====")
        
        element_symbol, atomic_number = self._get_element_input()
        
        orbital_types = ["s", "p", "d", "f"]
        print("\nChoose orbital type:")
        for i, orbital in enumerate(orbital_types, 1):
            print(f"[{i}] {orbital}-orbitals")
        
        while True:
            choice = input("\nEnter choice (1-4): ").strip()
            if choice.isdigit() and 1 <= int(choice) <= 4:
                orbital_type = orbital_types[int(choice) - 1]
                break
            print("Invalid choice. Please enter a number between 1 and 4.")
        
        print(f"\nCounting electrons in {orbital_type}-orbitals for {element_symbol}...")
        time.sleep(0.5)  # Brief pause for effect
        
        result = solve_electron_problem("orbital_count", element_symbol, orbital_type=orbital_type)
        
        if result["success"]:
            display_results_header()
            print(f"Element: {result['element']} (atomic number {atomic_number})")
            print(f"Electron configuration: {result['configuration']}")
            print(f"Total electrons in {orbital_type}-orbitals: {result[f'{orbital_type}_electrons']}")
            
            # Filtering steps - show from step 3 onwards (skip first two)
            steps_to_display = result["steps"][2:] if len(result["steps"]) > 2 else []
            display_steps(steps_to_display)
        else:
            print(f"Error: {result['error']}")
    
    def _solve_electron_config(self):
        """Generate electron configuration for an element."""
        print("\n===== ELECTRON CONFIGURATION GENERATOR =====")
        
        element_symbol, atomic_number = self._get_element_input()
        
        print(f"\nGenerating electron configuration for {element_symbol}...")
        time.sleep(0.5)  # Brief pause for effect
        
        result = solve_electron_problem("electron_config", element_symbol)
        
        if result["success"]:
            display_results_header()
            print(f"Element: {result['element']} (atomic number {atomic_number})")
            print(f"Full electron configuration: {result['full_configuration']}")
            print(f"Shorthand electron configuration: {result['shorthand_configuration']}")
            
            if result["is_anomalous"]:
                print("\nNote: This element has an anomalous configuration due to stability.")
            
            # Filtering steps - show from step 2 onwards (skip first step)
            steps_to_display = result["steps"][1:] if len(result["steps"]) > 1 else []
            display_steps(steps_to_display)
        else:
            print(f"Error: {result['error']}")
    
    def _solve_quantum_numbers(self):
        """Analyze quantum numbers for a specific electron."""
        print("\n===== QUANTUM NUMBER ANALYZER =====")
        
        element_symbol, atomic_number = self._get_element_input()
        
        while True:
            electron_input = input(f"\nWhich electron to analyze (1-{atomic_number}, or press Enter for last electron): ").strip()
            if not electron_input:
                electron_number = atomic_number
                break
            elif electron_input.isdigit() and 1 <= int(electron_input) <= atomic_number:
                electron_number = int(electron_input)
                break
            print(f"Invalid input. Please enter a number between 1 and {atomic_number}.")
        
        print(f"\nAnalyzing quantum numbers for electron #{electron_number} in {element_symbol}...")
        time.sleep(0.5)  # Brief pause for effect
        
        result = solve_electron_problem("quantum_numbers", element_symbol, electron_number=electron_number)
        
        if result["success"]:
            display_results_header()
            print(f"Element: {result['element']} (atomic number {atomic_number})")
            print(f"Analyzing electron #{electron_number}")
            print(f"This electron is in the {result['orbital']} orbital")
            
            print("\nQuantum numbers:")
            print(f"  Principal quantum number (n): {result['n']}")
            print(f"  Azimuthal quantum number (l): {result['l']}")
            print(f"  Magnetic quantum number (ml): possible values {result['ml_possibilities']}")
            print(f"  Spin quantum number (ms): possible values {result['ms_possibilities']}")
            
            # Filtering steps - show from step 3 onwards (skip first two)
            steps_to_display = result["steps"][2:] if len(result["steps"]) > 2 else []
            display_steps(steps_to_display)
        else:
            print(f"Error: {result['error']}")
    
    def _solve_molecular_orbitals(self):
        """Solve molecular orbital problems."""
        print("\n===== MOLECULAR ORBITAL ANALYZER =====")
        
        # List of common diatomic molecules
        common_molecules = ["H2", "He2", "Li2", "B2", "C2", "N2", "O2", "F2", "Ne2", "NO", "CO"]
        
        print("\nCommon diatomic molecules:")
        for i, molecule in enumerate(common_molecules, 1):
            print(f"[{i}] {molecule}")
        
        while True:
            choice = input("\nEnter molecule number or custom formula: ").strip()
            if choice.isdigit() and 1 <= int(choice) <= len(common_molecules):
                formula = common_molecules[int(choice) - 1]
                break
            elif choice:  # Custom input
                formula = choice
                break
            print("Invalid input. Please try again.")
        
        print("\nChoose problem type:")
        print("[1] Calculate bond order")
        print("[2] Determine paramagnetic/diamagnetic nature")
        
        while True:
            problem_choice = input("\nEnter choice (1-2): ").strip()
            if problem_choice == "1":
                problem_type = "bond_order"
                break
            elif problem_choice == "2":
                problem_type = "paramagnetic"
                break
            print("Invalid choice. Please enter 1 or 2.")
        
        print(f"\nAnalyzing {formula} for {problem_type.replace('_', ' ')}...")
        time.sleep(0.5)  # Brief pause for effect
        
        result = solve_molecular_orbital_problem(formula, problem_type)
        
        if result["success"]:
            display_results_header()
            print(f"Molecule: {result['molecule']}")
            print(f"Molecular orbital configuration: {result['mo_configuration']}")
            
            if problem_type == "bond_order":
                print(f"Bond order: {result['bond_order']}")
            else:
                paramag_result = "Paramagnetic" if result['paramagnetic'] else "Diamagnetic"
                print(f"Magnetic property: {paramag_result}")
            
            # Display all steps since there's no need to skip
            display_steps(result["steps"])
        else:
            print(f"Error: {result['error']}")