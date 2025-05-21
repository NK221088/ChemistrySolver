"""
Terminal User Interface for Acid-Base Analysis
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user
from chemistry_solver.acid_base import identify_acid_base, analyze_compound_list

class AcidBaseUI:
    """UI class for acid-base chemistry analysis."""
    
    def __init__(self):
        self.title = "ACID-BASE ANALYZER"
    
    def run(self):
        """Run the acid-base UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-3): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_single_compound()
            elif choice == "2":
                self._handle_multiple_compounds()
            elif choice == "3":
                self._handle_acid_base_info()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the acid-base module menu."""
        menu = """
        [1] Identify acid/base for single compound
        [2] Analyze multiple compounds
        [3] Learn about acid-base theory
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_single_compound(self):
        """Handle acid/base identification for a single compound."""
        print("\n===== ACID/BASE IDENTIFIER =====")
        
        try:
            compound = input("Enter chemical formula: ")
            
            result = identify_acid_base(compound)
            
            display_results_header()
            print(f"Compound: {compound}")
            print(f"Classification: {result['classification']}")
            print(f"Explanation: {result['explanation']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_multiple_compounds(self):
        """Handle analysis of multiple compounds."""
        print("\n===== COMPOUND LIST ANALYZER =====")
        
        try:
            compounds_input = input("Enter chemical formulas (comma-separated): ")
            compounds = [comp.strip() for comp in compounds_input.split(',')]
            
            result = analyze_compound_list(compounds)
            
            display_results_header()
            print(f"Analyzing compounds: {', '.join(compounds)}")
            print("-" * 50)
            
            for compound_result in result["compounds"]:
                print(f"Compound: {compound_result['compound']}")
                print(f"Classification: {compound_result['classification']}")
                print(f"Explanation: {compound_result['explanation']}")
                print("-" * 50)
            
            if result["all_acids"]:
                print("RESULT: All compounds in this list are acids.")
            else:
                print("RESULT: Not all compounds in this list are acids.")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_acid_base_info(self):
        """Display educational information about acid-base chemistry."""
        print("\n===== ACID-BASE THEORY =====")
        
        info = """
Acid-Base Definitions:

1. Brønsted-Lowry Definition:
   - Acid: A substance that donates a proton (H+)
   - Base: A substance that accepts a proton (H+)

2. Lewis Definition:
   - Acid: An electron pair acceptor
   - Base: An electron pair donor

3. Arrhenius Definition:
   - Acid: A substance that increases H+ concentration in water
   - Base: A substance that increases OH- concentration in water

Common Acid-Base Properties:

- Acids typically:
  • Taste sour (e.g., citrus fruits)
  • Turn blue litmus paper red
  • Have pH values less than 7
  • React with metals to produce hydrogen gas
  • React with carbonates to produce carbon dioxide

- Bases typically:
  • Taste bitter
  • Feel slippery or soapy
  • Turn red litmus paper blue
  • Have pH values greater than 7
  • React with acids to form salts and water (neutralization)

Strong vs. Weak Acids/Bases:
- Strong acids/bases completely dissociate in solution
- Weak acids/bases partially dissociate in solution
        """
        
        print(info)