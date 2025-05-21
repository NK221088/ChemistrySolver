"""
Terminal User Interface for Redox Reactions Analysis
"""
from ui.terminal_ui import display_title, display_results_header, display_steps, wait_for_user
from chemistry_solver.redox_reactions import (
    parse_redox_reaction, 
    determine_redox_favorability, 
    determine_favorability_from_potentials,
    get_half_reaction_by_element
)

class RedoxUI:
    """UI class for redox reactions analysis."""
    
    def __init__(self):
        self.title = "REDOX REACTIONS ANALYZER"
    
    def run(self):
        """Run the redox reactions UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-4): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_redox_analysis()
            elif choice == "2":
                self._handle_favorability_calculator()
            elif choice == "3":
                self._handle_half_reaction_lookup()
            elif choice == "4":
                self._handle_about_redox()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the redox reactions module menu."""
        menu = """
        [1] Analyze redox reaction
        [2] Calculate reaction favorability from potentials
        [3] Look up half-reactions by element
        [4] About redox reactions
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_redox_analysis(self):
        """Handle redox reaction analysis."""
        print("\n===== REDOX REACTION ANALYZER =====")
        
        try:
            print("Enter a redox reaction (use -> for reaction arrow)")
            print("Example: Zn(s) + Cu2+(aq) -> Zn2+(aq) + Cu(s)")
            reaction = input("Reaction: ")
            
            result = determine_redox_favorability(reaction)
            
            display_results_header()
            if result['favorable'] is not None:
                print(f"Reaction is {result['message']}")
                print(f"Oxidation half-reaction: {result['oxidation_half']}")
                print(f"Reduction half-reaction: {result['reduction_half']}")
                
                if 'e_cell' in result:
                    print(f"Standard cell potential: {result['e_cell']:.3f} V")
                
                display_steps(result['steps'])
            else:
                print(f"Unable to automatically analyze: {result.get('message', 'Insufficient information')}")
                print("\nWould you like to enter the half-reactions manually? (y/n)")
                choice = input("> ").strip().lower()
                
                if choice == 'y':
                    print("\nEnter oxidation half-reaction:")
                    oxidation_half = input("> ")
                    print("Enter reduction half-reaction:")
                    reduction_half = input("> ")
                    
                    print("\nDo you know the standard reduction potentials for these half-reactions? (y/n)")
                    pot_choice = input("> ").strip().lower()
                    
                    if pot_choice == 'y':
                        try:
                            print("\nEnter standard reduction potential for oxidation half-reaction (V):")
                            oxidation_potential = float(input("> "))
                            print("Enter standard reduction potential for reduction half-reaction (V):")
                            reduction_potential = float(input("> "))
                            
                            # Calculate favorability
                            favor_result = determine_favorability_from_potentials(oxidation_potential, reduction_potential)
                            
                            print("\n" + "-"*50)
                            print("MANUAL ANALYSIS RESULTS".center(50))
                            print("-"*50)
                            
                            print(f"Oxidation half-reaction: {oxidation_half}")
                            print(f"Reduction half-reaction: {reduction_half}")
                            print(f"Reaction is {favor_result['message']}")
                            print(f"Standard cell potential: {favor_result['e_cell']:.3f} V")
                            
                            display_steps(favor_result['steps'])
                        except ValueError:
                            print("Error: Please enter valid numerical values for potentials.")
                    else:
                        print("\nWithout reduction potentials, favorability cannot be determined.")
                        print("You can use option [2] from the menu to calculate favorability if you find the potentials.")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_favorability_calculator(self):
        """Handle redox favorability calculator."""
        print("\n===== REDOX FAVORABILITY CALCULATOR =====")
        
        try:
            print("Enter standard reduction potentials for half-reactions")
            oxidation_potential = float(input("Oxidation half-reaction potential (V): "))
            reduction_potential = float(input("Reduction half-reaction potential (V): "))
            
            result = determine_favorability_from_potentials(oxidation_potential, reduction_potential)
            
            display_results_header()
            print(f"Reaction is {result['message']}")
            print(f"Standard cell potential: {result['e_cell']:.3f} V")
            
            display_steps(result['steps'])
                
        except ValueError:
            print("Error: Please enter valid numerical values for potentials.")
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_half_reaction_lookup(self):
        """Handle half-reaction lookup by element."""
        print("\n===== HALF-REACTION LOOKUP =====")
        print("This tool helps find standard reduction potentials for half-reactions")
        print("For complex redox reactions like SO2 + MnO4- + H2O -> SO42- + Mn2+, look up")
        print("the individual elements (e.g., 'Mn' for MnO4- and 'S' for SO2/SO42-)")
        
        element = input("\nEnter element symbol (e.g., Cu, Fe, Zn, Mn, S): ")
        state = input("Enter state (optional - s, aq, g): ")
        
        half_reactions = get_half_reaction_by_element(element, state)
        
        display_results_header()
        if half_reactions:
            print(f"Found {len(half_reactions)} half-reactions for {element}:")
            for i, reaction in enumerate(half_reactions, 1):
                print(f"\n{i}. {reaction['half_reaction']}")
                print(f"   Standard Reduction Potential: {reaction['potential']} V")
                
            print("\nTo use these in analysis:")
            print("1. Identify the relevant half-reactions for your redox process")
            print("2. Use option [2] in the menu to manually calculate favorability")
            print("   using the potentials shown above")
        else:
            print(f"No half-reactions found for {element} in state {state if state else 'any'}.")
            print("Try searching for the base element (e.g., 'Mn' instead of 'MnO4-')")
            print("or check if the element symbol is correct.")
    
    def _handle_about_redox(self):
        """Display information about redox reactions."""
        print("\n===== ABOUT REDOX REACTIONS =====")
        
        info = """
Redox (reduction-oxidation) reactions are chemical reactions where electrons 
are transferred between species, resulting in changes to their oxidation states.

Key Concepts:
- Oxidation: Loss of electrons (increase in oxidation state)
- Reduction: Gain of electrons (decrease in oxidation state)
- Oxidizing agent: The species that causes oxidation (and is itself reduced)
- Reducing agent: The species that causes reduction (and is itself oxidized)

Standard Reduction Potentials:
- Measured in volts (V) relative to the standard hydrogen electrode (SHE)
- Higher (more positive) potential indicates greater tendency to be reduced
- Cell potential (E°cell) = E°reduction - E°oxidation
- A positive E°cell indicates a spontaneous reaction

Activity Series:
- Metals higher in the series are more reactive (better reducing agents)
- A metal can displace ions of metals lower in the series from solution

Balancing Complex Redox Reactions:
1. Split the reaction into half-reactions (oxidation and reduction)
2. Balance atoms other than O and H
3. Balance O atoms by adding H2O
4. Balance H atoms by adding H+
5. Balance charge by adding electrons
6. Multiply half-reactions to equalize electrons transferred
7. Add half-reactions and cancel out common terms

Example: SO2 + MnO4- + H2O -> SO42- + Mn2+ + H+
Oxidation: SO2 + 2H2O -> SO42- + 4H+ + 2e-
Reduction: MnO4- + 8H+ + 5e- -> Mn2+ + 4H2O

To balance: Multiply oxidation by 5, reduction by 2, then add.

For analyzing such reactions:
- Look up standard reduction potentials for each half-reaction
- Calculate E°cell = E°reduction - E°oxidation
- Positive E°cell indicates spontaneous reaction
        """
        
        print(info)