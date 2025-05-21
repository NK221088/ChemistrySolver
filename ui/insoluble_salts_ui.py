"""
Terminal User Interface for Insoluble Salts & Qualitative Analysis
"""
from ui.terminal_ui import display_title, display_results_header, display_steps, wait_for_user
from chemistry_solver.insoluble_salts import (
    solve_qualitative_analysis_problem,
    analyze_specific_scenario,
    get_available_cations,
    get_available_reagents
)

class InsolubleSaltsUI:
    """UI class for insoluble salts and qualitative analysis."""
    
    def __init__(self):
        self.title = "INSOLUBLE SALTS & QUALITATIVE ANALYSIS"
    
    def run(self):
        """Run the insoluble salts and qualitative analysis UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-3): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_qualitative_analysis()
            elif choice == "2":
                self._handle_predefined_scenarios()
            elif choice == "3":
                self._handle_solubility_rules()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the insoluble salts module menu."""
        menu = """
        [1] Qualitative Analysis Problem Solver
        [2] Analyze Predefined Scenarios
        [3] Solubility Rules Reference
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_qualitative_analysis(self):
        """Handle qualitative analysis problems."""
        print("\n===== QUALITATIVE ANALYSIS PROBLEM SOLVER =====")
        
        try:
            # Print available cations and reagents
            print("\nAvailable cations:")
            available_cations = get_available_cations()
            print(", ".join(available_cations))
            
            print("\nAvailable reagents:")
            available_reagents = get_available_reagents()
            print(", ".join(available_reagents))
            
            # Get candidate cations
            cation_input = input("\nEnter possible cations (comma-separated): ")
            cation_candidates = [c.strip() for c in cation_input.split(",")]
            
            # Get reagents that cause precipitation
            precip_input = input("Enter reagents that cause precipitation (comma-separated): ")
            if precip_input.strip():
                precipitates_with = [r.strip() for r in precip_input.split(",")]
            else:
                precipitates_with = []
            
            # Get reagents that DON'T cause precipitation
            no_precip_input = input("Enter reagents that DON'T cause precipitation (comma-separated): ")
            if no_precip_input.strip():
                no_precipitate_with = [r.strip() for r in no_precip_input.split(",")]
            else:
                no_precipitate_with = []
            
            # Solve the problem
            result = solve_qualitative_analysis_problem(
                cation_candidates,
                precipitates_with,
                no_precipitate_with
            )
            
            # Display results
            display_results_header()
            
            if "error" in result:
                print(f"Error: {result['error']}")
                if "available_cations" in result:
                    print("Available cations:", ", ".join(result["available_cations"]))
                if "available_reagents" in result:
                    print("Available reagents:", ", ".join(result["available_reagents"]))
                return
            
            # Display analysis steps
            display_steps(result["steps"])
            
            print("\nConclusion:")
            print(result["conclusion"])
            
            if result["identified_cations"]:
                print("\nIdentified cation(s):", ", ".join(result["identified_cations"]))
            else:
                print("\nNo cation could be identified with the given constraints.")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_predefined_scenarios(self):
        """Handle predefined analysis scenarios."""
        print("\n===== PREDEFINED ANALYSIS SCENARIOS =====")
        
        scenarios = {
            "1": {
                "id": "W20_8",
                "description": "Determine if unknown solution contains Ag+, Ba2+, or Pb2+"
            }
            # Additional scenarios can be added here
        }
        
        print("\nAvailable scenarios:")
        for key, scenario in scenarios.items():
            print(f"  [{key}] {scenario['description']}")
        
        choice = input("\nSelect scenario (or 0 to go back): ")
        
        if choice == "0":
            return
        elif choice in scenarios:
            scenario_id = scenarios[choice]["id"]
            result = analyze_specific_scenario(scenario_id)
            
            display_results_header()
            
            if "error" in result:
                print(f"Error: {result['error']}")
                return
                
            # Display analysis steps
            display_steps(result["steps"])
            
            print("\nConclusion:")
            print(result["conclusion"])
            
            if result["identified_cations"]:
                print("\nIdentified cation(s):", ", ".join(result["identified_cations"]))
            else:
                print("\nNo cation could be identified with the given constraints.")
        else:
            print("Invalid choice.")
    
    def _handle_solubility_rules(self):
        """Display information about solubility rules."""
        print("\n===== SOLUBILITY RULES REFERENCE =====")
        
        solubility_rules = """
General Solubility Rules:
------------------------
1. Most salts containing Group 1 elements (Li+, Na+, K+, etc.) are soluble.
2. Most salts containing ammonium ion (NH4+) are soluble.
3. Most nitrates (NO3-), acetates (CH3COO-), and perchlorates (ClO4-) are soluble.
4. Most chlorides (Cl-), bromides (Br-), and iodides (I-) are soluble.
   Exceptions: AgCl, PbCl2, Hg2Cl2, HgCl2 are insoluble or slightly soluble.
5. Most sulfates (SO4^2-) are soluble.
   Exceptions: CaSO4, BaSO4, PbSO4, SrSO4, Ag2SO4 are insoluble or slightly soluble.
6. Most hydroxides (OH-) are insoluble.
   Exceptions: Group 1 hydroxides and Ba(OH)2 are soluble.
7. Most sulfides (S^2-), carbonates (CO3^2-), phosphates (PO4^3-),
   and chromates (CrO4^2-) are insoluble.
   Exceptions: Group 1 and NH4+ salts of these anions are soluble.

Common Precipitation Reagents:
----------------------------
- AgNO3: Test for halides (Cl-, Br-, I-)
- BaCl2: Test for sulfates (SO4^2-)
- NaOH: Test for metal hydroxides
- H2SO4: Test for Ba^2+, Pb^2+
- HCl: Test for Ag+, Pb^2+
- (NH4)2CO3: Test for carbonates
        """
        
        print(solubility_rules)