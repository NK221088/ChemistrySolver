"""
Terminal User Interface for Redox Reactions Analysis
"""
from ui.terminal_ui import display_title, display_results_header, display_steps, wait_for_user
from chemistry_solver.redox_reactions import (
    parse_redox_reaction, 
    determine_redox_favorability, 
    determine_favorability_from_potentials,
    get_half_reaction_by_element,
    find_standard_reduction_potential,
    calculate_nernst_equation,
    get_common_redox_reaction,
    get_redox_pair,
    enhanced_determine_redox_favorability,
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
            choice = input("\nEnter choice (0-6): ").strip()
            
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
                self._handle_common_redox_lookup()
            elif choice == "5":
                self._handle_nernst_equation()
            elif choice == "6":
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
        [4] Look up common redox reactions
        [5] Calculate potential using Nernst equation
        [6] About redox reactions
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_redox_analysis(self):
        """Handle redox reaction analysis."""
        print("\n===== REDOX REACTION ANALYZER =====")
        
        try:
            print("Enter a redox reaction (use -> for reaction arrow)")
            print("Example 1: Zn(s) + Cu2+(aq) -> Zn2+(aq) + Cu(s)")
            print("Example 2: SO2 + MnO4- + H2O -> SO4^2- + Mn2+ + H+")
            reaction = input("Reaction: ")
            
            # Use the enhanced parser for complex reactions
            result = enhanced_determine_redox_favorability(reaction)
            
            display_results_header()
            if result['favorable'] is not None:
                print(f"Reaction is {result['message']}")
                
                # Check if oxidation_half is available in the result
                if 'oxidation_half' in result and result['oxidation_half']:
                    print(f"Oxidation half-reaction: {result['oxidation_half']}")
                
                # Check if reduction_half is available in the result
                if 'reduction_half' in result and result['reduction_half']:
                    print(f"Reduction half-reaction: {result['reduction_half']}")
                
                # Check if balanced equation is available
                if 'balanced_equation' in result and result['balanced_equation']:
                    print(f"Balanced equation: {result['balanced_equation']}")
                
                # Display oxidizing and reducing agents if available
                if 'oxidizing_agent' in result and result['oxidizing_agent']:
                    print(f"Oxidizing agent: {result['oxidizing_agent']}")
                
                if 'reducing_agent' in result and result['reducing_agent']:
                    print(f"Reducing agent: {result['reducing_agent']}")
                
                # Check if e_cell is available in the result
                if 'e_cell' in result:
                    print(f"Standard cell potential: {result['e_cell']:.3f} V")
                
                # Display notes if available
                if 'notes' in result and result['notes']:
                    print(f"\nNotes: {result['notes']}")
                
                display_steps(result['steps'])
            else:
                print(f"Unable to automatically analyze: {result.get('message', 'Insufficient information')}")
                
                # For complex reactions that need manual completion
                if 'oxidation_half' in result and result['oxidation_half'].startswith('[Please'):
                    print("\nThis appears to be a complex redox reaction that requires manual balancing.")
                    print("Would you like help with balancing this reaction manually? (y/n)")
                    help_choice = input("> ").strip().lower()
                    
                    if help_choice == 'y':
                        print("\n===== MANUAL REDOX REACTION BALANCING =====")
                        print("Follow these steps to balance your complex redox reaction:")
                        print("1. Identify what's being oxidized and reduced (look at oxidation states)")
                        print("2. Write separate half-reactions for oxidation and reduction")
                        print("3. Balance atoms other than O and H in each half-reaction")
                        print("4. Balance O by adding H2O")
                        print("5. Balance H by adding H+")
                        print("6. Balance charge by adding electrons")
                        print("7. Multiply half-reactions to equalize electrons transferred")
                        print("8. Add the half-reactions together and cancel out common terms")
                        
                        print("\nFor your reaction:", reaction)
                        print("Let's identify the species involved:")
                        
                        # Try to identify reactants and products
                        parts = reaction.split("->")
                        if len(parts) == 2:
                            reactants = parts[0].strip().split("+")
                            products = parts[1].strip().split("+")
                            
                            print("\nReactants:", ", ".join(r.strip() for r in reactants))
                            print("Products:", ", ".join(p.strip() for p in products))
                            
                            # For permanganate reactions with SO2 (like your example)
                            if any("MnO4" in r for r in reactants) and any("SO2" in r for r in reactants):
                                print("\nThis appears to be the oxidation of sulfur dioxide by permanganate.")
                                print("Oxidation half-reaction: SO2 + 2H2O → SO4²⁻ + 4H⁺ + 2e⁻")
                                print("Reduction half-reaction: MnO4⁻ + 8H⁺ + 5e⁻ → Mn²⁺ + 4H2O")
                                print("\nTo balance electrons: multiply oxidation by 5 and reduction by 2")
                                print("5(SO2 + 2H2O → SO4²⁻ + 4H⁺ + 2e⁻)")
                                print("2(MnO4⁻ + 8H⁺ + 5e⁻ → Mn²⁺ + 4H2O)")
                                print("\nAdd these equations:")
                                print("5SO2 + 10H2O + 2MnO4⁻ + 16H⁺ → 5SO4²⁻ + 20H⁺ + 2Mn²⁺ + 8H2O")
                                print("\nSimplify (cancel 16H⁺ and 8H2O):")
                                print("5SO2 + 2H2O + 2MnO4⁻ → 5SO4²⁻ + 4H⁺ + 2Mn²⁺")
                                print("\nFinal balanced equation: 5SO2 + 2MnO4⁻ + 2H2O → 5SO4²⁻ + 4H⁺ + 2Mn²⁺")
                
                print("\nWould you like to enter the half-reactions manually? (y/n)")
                choice = input("> ").strip().lower()
                
                if choice == 'y':
                    print("\nEnter oxidation half-reaction:")
                    oxidation_half = input("> ")
                    print("Enter reduction half-reaction:")
                    reduction_half = input("> ")
                    
                    # Try to find standard reduction potentials
                    oxidation_potential = find_standard_reduction_potential(
                        oxidation_half.replace("->", "←").replace("+", "-").replace("2e⁻", "")
                    )
                    reduction_potential = find_standard_reduction_potential(reduction_half)
                    
                    if oxidation_potential is not None and reduction_potential is not None:
                        print(f"\nFound standard reduction potentials:")
                        print(f"Oxidation half-reaction potential: {oxidation_potential} V")
                        print(f"Reduction half-reaction potential: {reduction_potential} V")
                        
                        # Calculate favorability automatically
                        favor_result = determine_favorability_from_potentials(oxidation_potential, reduction_potential)
                        
                        print("\n" + "-"*50)
                        print("AUTOMATIC ANALYSIS RESULTS".center(50))
                        print("-"*50)
                        
                        print(f"Oxidation half-reaction: {oxidation_half}")
                        print(f"Reduction half-reaction: {reduction_half}")
                        print(f"Reaction is {favor_result['message']}")
                        print(f"Standard cell potential: {favor_result['e_cell']:.3f} V")
                        
                        display_steps(favor_result['steps'])
                    else:
                        print("\nCouldn't find standard reduction potentials automatically.")
                        print("Do you know the standard reduction potentials for these half-reactions? (y/n)")
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
            
            # Add option to check a specific half-reaction
            print("\nWould you like to look up a specific half-reaction? (y/n)")
            choice = input("> ").strip().lower()
            
            if choice == 'y':
                print("\nEnter the half-reaction:")
                half_reaction = input("> ")
                potential = find_standard_reduction_potential(half_reaction)
                
                if potential is not None:
                    print(f"\nStandard reduction potential for '{half_reaction}': {potential} V")
                else:
                    print(f"\nCould not find standard reduction potential for '{half_reaction}'")
        else:
            print(f"No half-reactions found for {element} in state {state if state else 'any'}.")
            print("Try searching for the base element (e.g., 'Mn' instead of 'MnO4-')")
            print("or check if the element symbol is correct.")
    
    def _handle_common_redox_lookup(self):
        """Handle lookup of common redox reactions and pairs."""
        print("\n===== COMMON REDOX REACTIONS LOOKUP =====")
        print("This tool helps find information about common redox reactions and pairs")
        
        print("\nLookup options:")
        print("[1] Common redox reaction by name")
        print("[2] Common redox pair by name")
        print("[0] Back to main redox menu")
        
        choice = input("\nEnter choice (0-2): ").strip()
        
        if choice == "0":
            return
        elif choice == "1":
            print("\nEnter the name of a common redox reaction:")
            print("Examples: Metal displacement, Combustion, Silver mirror test, etc.")
            name = input("> ")
            
            reaction = get_common_redox_reaction(name)
            if reaction:
                print("\n" + "-"*50)
                print(f"COMMON REDOX REACTION: {name}".center(50))
                print("-"*50)
                
                print(f"Reaction: {reaction.get('reaction', 'N/A')}")
                print(f"Balanced equation: {reaction.get('balanced_equation', 'N/A')}")
                print(f"Oxidation half: {reaction.get('oxidation_half', 'N/A')}")
                print(f"Reduction half: {reaction.get('reduction_half', 'N/A')}")
                
                if 'notes' in reaction:
                    print(f"\nNotes: {reaction['notes']}")
            else:
                print(f"No common redox reaction found with name '{name}'")
                
        elif choice == "2":
            print("\nEnter the name of a common redox pair:")
            print("Examples: Copper-Zinc, Iron-Magnesium, etc.")
            name = input("> ")
            
            pair = get_redox_pair(name)
            if pair:
                print("\n" + "-"*50)
                print(f"COMMON REDOX PAIR: {name}".center(50))
                print("-"*50)
                
                print(f"Net reaction: {pair.get('net_reaction', 'N/A')}")
                print(f"Oxidation: {pair.get('oxidation', 'N/A')}")
                print(f"Reduction: {pair.get('reduction', 'N/A')}")
                print(f"Standard cell potential: {pair.get('e_cell', 'N/A')} V")
                print(f"Favorable: {'Yes' if pair.get('favorable', False) else 'No'}")
                
                if 'notes' in pair:
                    print(f"\nNotes: {pair['notes']}")
            else:
                print(f"No common redox pair found with name '{name}'")
        else:
            print("Invalid choice.")
    
    def _handle_nernst_equation(self):
        """Handle Nernst equation calculations."""
        print("\n===== NERNST EQUATION CALCULATOR =====")
        print("This tool calculates cell potential under non-standard conditions")
        
        try:
            e_cell_standard = float(input("Enter standard cell potential E°cell (V): "))
            n = int(input("Enter number of electrons transferred: "))
            concentration_ratio = float(input("Enter concentration ratio [products]/[reactants]: "))
            temp_choice = input("Use standard temperature (25°C)? (y/n): ").strip().lower()
            
            if temp_choice == 'n':
                temp_celsius = float(input("Enter temperature (°C): "))
                temperature = temp_celsius + 273.15
            else:
                temperature = 298.15  # 25°C in Kelvin
            
            e_cell = calculate_nernst_equation(e_cell_standard, n, concentration_ratio, temperature)
            
            display_results_header()
            print(f"Standard cell potential (E°cell): {e_cell_standard:.3f} V")
            print(f"Cell potential under given conditions (Ecell): {e_cell:.3f} V")
            
            # Show calculation steps
            steps = [
                f"1. Standard cell potential (E°cell): {e_cell_standard:.3f} V",
                f"2. Number of electrons transferred (n): {n}",
                f"3. Concentration ratio [products]/[reactants]: {concentration_ratio}",
                f"4. Temperature: {temperature-273.15:.1f}°C ({temperature:.2f} K)",
                f"5. Using the Nernst equation: Ecell = E°cell - ((RT)/(nF))·ln(Q)",
                f"   where R = 8.314 J/(mol·K), F = 96485 C/mol, Q = concentration ratio",
                f"6. Ecell = {e_cell_standard:.3f} - ((8.314 × {temperature:.2f})/(${n} × 96485))·2.303·log10({concentration_ratio})",
                f"7. Ecell = {e_cell:.3f} V"
            ]
            
            display_steps(steps)
            
            # Provide interpretation
            if e_cell > 0:
                print("\nThe reaction is favorable under these conditions (Ecell > 0).")
            else:
                print("\nThe reaction is not favorable under these conditions (Ecell ≤ 0).")
            
        except ValueError:
            print("Error: Please enter valid numerical values.")
        except Exception as e:
            print(f"Error: {str(e)}")
    
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

The Nernst Equation:
- Used to calculate cell potential under non-standard conditions
- Ecell = E°cell - ((RT)/(nF))·ln(Q)
  where:
  - R = gas constant (8.314 J/(mol·K))
  - T = temperature in Kelvin
  - n = number of electrons transferred
  - F = Faraday constant (96485 C/mol)
  - Q = reaction quotient ([products]/[reactants])

For analyzing reactions:
- Look up standard reduction potentials for each half-reaction
- Calculate E°cell = E°reduction - E°oxidation
- Positive E°cell indicates spontaneous reaction
        """
        
        print(info)