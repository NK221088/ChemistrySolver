"""
Terminal User Interface for Thermodynamics Calculations
"""
from ui.terminal_ui import display_title, display_results_header, display_steps, wait_for_user
from chemistry_solver.enthalpy import solve_enthalpy_problem

class ThermodynamicsUI:
    """UI class for thermodynamics calculations."""
    
    def __init__(self):
        self.title = "THERMODYNAMICS CALCULATOR"
    
    def run(self):
        """Run the thermodynamics UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-3): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_hess_law_calculation()
            elif choice == "2":
                self._handle_enthalpy_calculation()
            elif choice == "3":
                self._handle_gibbs_energy_calculation()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the thermodynamics module menu."""
        menu = """
        [1] Hess's Law Calculations
        [2] Enthalpy of Reaction
        [3] Gibbs Free Energy
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_hess_law_calculation(self):
        """Handle Hess's Law calculation problems."""
        print("\n===== HESS'S LAW CALCULATOR =====")
        print("\nEnter a Hess's Law problem. Format as follows:")
        print("1. List known reactions with enthalpies")
        print("2. Include target reaction without enthalpy")
        print("Example:")
        print("""
C (graphite) + O2 (g) → CO2 (g) ΔH° = -393.5 kJ/mol
H2 (g) + 0.5 O2 (g) → H2O (l) ΔH° = -285.8 kJ/mol
C (graphite) + 2 H2 (g) + 0.5 O2 (g) → CH3OH (l) ΔH° = -238.7 kJ/mol
CH3OH (l) + 1.5 O2 (g) → CO2 (g) + 2 H2O (l)
        """)
        
        print("\nEnter problem (press Enter twice when done):")
        problem_lines = []
        while True:
            line = input()
            if not line and problem_lines:
                break
            problem_lines.append(line)
        
        problem_text = "\n".join(problem_lines)
        
        try:
            result = solve_enthalpy_problem(problem_text)
            
            display_results_header()
            
            if result['enthalpy'] is not None:
                display_steps(result['steps'])
                print(f"\nFinal Answer: {result['enthalpy']:.1f} kJ/mol")
            else:
                print("\nCould not solve the problem. Please check your input.")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_enthalpy_calculation(self):
        """Handle enthalpy calculation from formation enthalpies."""
        print("\n===== ENTHALPY OF REACTION CALCULATOR =====")
        print("\nThis feature will be implemented in a future update.")
        # Future implementation
    
    def _handle_gibbs_energy_calculation(self):
        """Handle Gibbs free energy calculations."""
        print("\n===== GIBBS FREE ENERGY CALCULATOR =====")
        print("\nThis feature will be implemented in a future update.")
        # Future implementation