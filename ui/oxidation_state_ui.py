"""
Terminal User Interface for Oxidation State Calculator
"""
from ui.terminal_ui import display_title, display_results_header, display_steps, wait_for_user
from chemistry_solver.oxidation_state import calculate_oxidation_number, parse_formula

class OxidationStateUI:
    """UI class for oxidation state calculations."""
    
    def __init__(self):
        self.title = "OXIDATION STATE CALCULATOR"
    
    def run(self):
        """Run the oxidation state calculator UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-3): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_oxidation_calculation()
            elif choice == "2":
                self._handle_formula_analysis()
            elif choice == "3":
                self._handle_common_oxidation_states()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the oxidation state module menu."""
        menu = """
        [1] Calculate oxidation number
        [2] Analyze chemical formula
        [3] View common oxidation states
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_oxidation_calculation(self):
        """Handle oxidation number calculation."""
        print("\n===== OXIDATION NUMBER CALCULATOR =====")
        
        try:
            compound = input("Enter chemical compound (e.g., CrO2Cl2): ").strip()
            element = input("Enter element to find oxidation state for: ").strip()
            
            result = calculate_oxidation_number(compound, element)
            
            display_results_header()
            print(f"Compound: {result['compound']}")
            print(f"Element: {result['element']}")
            print(f"Oxidation Number: {result['oxidation_number']}")
            
            display_steps(result['steps'])
            
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_formula_analysis(self):
        """Handle chemical formula analysis."""
        print("\n===== CHEMICAL FORMULA ANALYSIS =====")
        
        try:
            formula = input("Enter chemical formula to analyze (e.g., Fe2(SO4)3): ").strip()
            
            # Currently only handling simple formulas without parentheses
            elements = parse_formula(formula)
            
            display_results_header()
            print(f"Formula: {formula}")
            print("\nElement Composition:")
            for element, count in elements.items():
                print(f"  {element}: {count}")
            
            print("\nNote: This analysis does not fully support complex formulas with parentheses.")
            
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_common_oxidation_states(self):
        """Display information about common oxidation states."""
        print("\n===== COMMON OXIDATION STATES =====")
        
        common_states = """
Common Element Oxidation States:

Alkali Metals (Group 1):
- Li, Na, K, Rb, Cs: +1

Alkaline Earth Metals (Group 2):
- Be, Mg, Ca, Sr, Ba: +2

Transition Metals (common states):
- Fe: +2, +3
- Cu: +1, +2
- Cr: +2, +3, +6
- Mn: +2, +3, +4, +7
- Co: +2, +3
- Ni: +2
- Zn: +2

Other Main Group Elements:
- B: +3
- Al: +3
- C: -4, +2, +4
- Si: +4
- N: -3, +3, +5
- P: -3, +3, +5
- O: -2 (except in peroxides: -1, superoxides: -1/2)
- S: -2, +4, +6
- F: -1
- Cl, Br, I: -1, +1, +3, +5, +7

Special Cases:
- H: +1 (in most compounds), -1 (in metal hydrides)
- O: -2 (typical), -1 (in peroxides), -1/2 (in superoxides)

Note: This is not an exhaustive list. Elements can have other oxidation states in specific compounds.
        """
        
        print(common_states)