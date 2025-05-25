"""
Terminal User Interface for Redox Reactions Analysis
Updated to work with the refactored redox_reactions module
"""
from ui.terminal_ui import display_title, display_results_header, display_steps, wait_for_user
from chemistry_solver.redox_reactions import (
    FavorabilityAnalyzer,
    RedoxParser,
    OxidationStateAnalyzer,
    PotentialCalculator,
    ComplexReactionHandler,
    ActivitySeriesAnalyzer,
    # Backward compatibility functions
    parse_redox_reaction,
    determine_redox_favorability,
    identify_redox_from_oxidation_states,
    get_half_reaction_by_element,
    find_standard_reduction_potential,
    calculate_nernst_equation
)
from chemistry_solver.redox_data import (
    get_common_redox_reaction,
    get_redox_pair,
    COMMON_REDOX_REACTIONS,
    COMMON_REDOX_PAIRS
)


class RedoxUI:
    """UI class for redox reactions analysis."""
    
    def __init__(self):
        self.title = "REDOX REACTIONS ANALYZER"
        self.favorability_analyzer = FavorabilityAnalyzer()
        self.redox_parser = RedoxParser()
        self.oxidation_analyzer = OxidationStateAnalyzer()
        self.potential_calculator = PotentialCalculator()
        self.activity_analyzer = ActivitySeriesAnalyzer()
    
    def run(self):
        """Run the redox reactions UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-8): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_redox_analysis()
            elif choice == "2":
                self._handle_favorability_calculator()
            elif choice == "3":
                self._handle_oxidation_state_analysis()
            elif choice == "4":
                self._handle_half_reaction_lookup()
            elif choice == "5":
                self._handle_common_redox_lookup()
            elif choice == "6":
                self._handle_nernst_equation()
            elif choice == "7":
                self._handle_complex_reaction_patterns()
            elif choice == "8":
                self._handle_about_redox()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the redox reactions module menu."""
        menu = """
        [1] Analyze redox reaction (comprehensive)
        [2] Calculate reaction favorability from potentials
        [3] Analyze oxidation state changes
        [4] Look up half-reactions by element
        [5] Look up common redox reactions
        [6] Calculate potential using Nernst equation
        [7] Complex reaction pattern recognition
        [8] About redox reactions
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_redox_analysis(self):
        """Handle comprehensive redox reaction analysis."""
        print("\n===== COMPREHENSIVE REDOX REACTION ANALYZER =====")
        
        try:
            print("Enter a redox reaction (use -> for reaction arrow)")
            print("Example 1: Zn(s) + Cu2+(aq) -> Zn2+(aq) + Cu(s)")
            print("Example 2: SO2 + MnO4- + H2O -> SO42- + Mn2+ + H+")
            print("Example 3: 2H2O2 -> 2H2O + O2")
            reaction = input("Reaction: ")
            
            # Use the new FavorabilityAnalyzer
            result = self.favorability_analyzer.determine_favorability(reaction)
            
            display_results_header()
            
            # Display main result
            if result.favorable is not None:
                print(f"Reaction is {result.message}")
                
                # Display half-reactions if available
                if result.oxidation_half and not result.oxidation_half.startswith("[Please"):
                    print(f"Oxidation half-reaction: {result.oxidation_half}")
                
                if result.reduction_half and not result.reduction_half.startswith("[Please"):
                    print(f"Reduction half-reaction: {result.reduction_half}")
                
                # Display balanced equation if available
                if result.balanced_equation:
                    print(f"Balanced equation: {result.balanced_equation}")
                
                # Display oxidizing and reducing agents
                if result.oxidizing_agent:
                    print(f"Oxidizing agent: {result.oxidizing_agent}")
                
                if result.reducing_agent:
                    print(f"Reducing agent: {result.reducing_agent}")
                
                # Display cell potential if available
                if result.e_cell is not None:
                    print(f"Standard cell potential: {result.e_cell:.3f} V")
                
                # Display notes if available
                if result.notes:
                    print(f"\nNotes: {result.notes}")
                
            else:
                print(f"Analysis result: {result.message}")
                
                # Check if this needs manual completion
                if (result.oxidation_half and result.oxidation_half.startswith("[Please")) or \
                   (result.reduction_half and result.reduction_half.startswith("[Please")):
                    self._handle_manual_balancing_help(reaction)
            
            # Display analysis steps
            if result.steps:
                display_steps(result.steps)
            
            # Offer additional analysis
            print("\nWould you like additional oxidation state analysis? (y/n)")
            if input("> ").strip().lower() == 'y':
                self._show_oxidation_state_details(reaction)
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_manual_balancing_help(self, reaction: str):
        """Provide help for manually balancing complex reactions."""
        print("\nThis appears to be a complex redox reaction that requires manual balancing.")
        print("Would you like step-by-step balancing guidance? (y/n)")
        
        if input("> ").strip().lower() == 'y':
            print("\n===== MANUAL REDOX REACTION BALANCING GUIDE =====")
            
            # Check for known complex patterns
            complex_handler = ComplexReactionHandler()
            complex_data = complex_handler.parse_complex_reaction(reaction)
            
            if not complex_data["oxidation_half"].startswith("[Please"):
                print(f"\nRecognized pattern! Here's the solution for: {reaction}")
                print(f"Oxidation half-reaction: {complex_data['oxidation_half']}")
                print(f"Reduction half-reaction: {complex_data['reduction_half']}")
                print(f"Balanced equation: {complex_data['balanced_equation']}")
                print(f"Oxidizing agent: {complex_data['oxidizing_agent']}")
                print(f"Reducing agent: {complex_data['reducing_agent']}")
                print(f"Electron transfer: {complex_data['electron_transfer']}")
                if complex_data.get('notes'):
                    print(f"Notes: {complex_data['notes']}")
            else:
                # Provide general balancing steps
                self._show_general_balancing_steps(reaction)
    
    def _show_general_balancing_steps(self, reaction: str):
        """Show general steps for balancing redox reactions."""
        print("\nGeneral steps for balancing redox reactions:")
        print("1. Identify what's being oxidized and reduced (check oxidation states)")
        print("2. Write separate half-reactions for oxidation and reduction")
        print("3. Balance atoms other than O and H in each half-reaction")
        print("4. Balance O by adding H2O")
        print("5. Balance H by adding H+")
        print("6. Balance charge by adding electrons")
        print("7. Multiply half-reactions to equalize electrons transferred")
        print("8. Add the half-reactions together and cancel out common terms")
        
        # Try to provide specific guidance based on reaction components
        parts = reaction.split("->")
        if len(parts) == 2:
            reactants = [r.strip() for r in parts[0].split("+")]
            products = [p.strip() for p in parts[1].split("+")]
            
            print(f"\nFor your reaction: {reaction}")
            print("Reactants:", ", ".join(reactants))
            print("Products:", ", ".join(products))
            
            # Provide specific hints for common reaction types
            self._provide_specific_hints(reactants, products)
    
    def _provide_specific_hints(self, reactants, products):
        """Provide specific hints based on reaction components."""
        reactant_str = " ".join(reactants).lower()
        product_str = " ".join(products).lower()
        
        if "mno4" in reactant_str and "mn2" in product_str:
            print("\nHint: This involves permanganate (MnO4-) being reduced to Mn2+")
            print("Reduction half: MnO4- + 8H+ + 5e- -> Mn2+ + 4H2O")
            
        if "so2" in reactant_str and "so4" in product_str:
            print("\nHint: This involves SO2 being oxidized to SO42-")
            print("Oxidation half: SO2 + 2H2O -> SO42- + 4H+ + 2e-")
            
        if "cr2o7" in reactant_str and "cr3" in product_str:
            print("\nHint: This involves dichromate (Cr2O7 2-) being reduced to Cr3+")
            print("Reduction half: Cr2O7 2- + 14H+ + 6e- -> 2Cr3+ + 7H2O")
            
        if "h2o2" in reactant_str and ("h2o" in product_str and "o2" in product_str):
            print("\nHint: This is hydrogen peroxide disproportionation")
            print("H2O2 acts as both oxidizing and reducing agent")
    
    def _show_oxidation_state_details(self, reaction: str):
        """Show detailed oxidation state analysis."""
        print("\n===== DETAILED OXIDATION STATE ANALYSIS =====")
        
        try:
            analysis = self.oxidation_analyzer.analyze_reaction(reaction)
            
            if analysis.is_redox:
                print("✓ This IS a redox reaction")
                print("\nOxidation state changes:")
                for item in analysis.oxidation_analysis:
                    print(f"  {item}")
                
                if analysis.oxidized_elements:
                    print(f"\nElements oxidized (lose electrons):")
                    for elem in analysis.oxidized_elements:
                        print(f"  {elem.element}: {elem.from_state} → {elem.to_state} "
                              f"(loses {elem.change} electrons) in {elem.reactant_compound} → {elem.product_compound}")
                
                if analysis.reduced_elements:
                    print(f"\nElements reduced (gain electrons):")
                    for elem in analysis.reduced_elements:
                        print(f"  {elem.element}: {elem.from_state} → {elem.to_state} "
                              f"(gains {elem.change} electrons) in {elem.reactant_compound} → {elem.product_compound}")
                
                # Check electron balance
                if analysis.is_balanced:
                    print(f"\n✓ Electron transfer is balanced:")
                    print(f"  Total electrons lost: {analysis.total_electrons_lost}")
                    print(f"  Total electrons gained: {analysis.total_electrons_gained}")
                else:
                    print(f"\n⚠ Electron transfer may be unbalanced:")
                    print(f"  Total electrons lost: {analysis.total_electrons_lost}")
                    print(f"  Total electrons gained: {analysis.total_electrons_gained}")
            else:
                print("✗ This is NOT a redox reaction (no electron transfer detected)")
                if analysis.oxidation_analysis:
                    print("\nOxidation state analysis:")
                    for item in analysis.oxidation_analysis:
                        print(f"  {item}")
                        
        except Exception as e:
            print(f"Error in oxidation state analysis: {str(e)}")
    
    def _handle_favorability_calculator(self):
        """Handle manual favorability calculation from potentials."""
        print("\n===== REDOX FAVORABILITY CALCULATOR =====")
        
        try:
            print("Enter standard reduction potentials for half-reactions")
            print("Note: Enter the potential for the oxidation half-reaction as written")
            print("      (the calculator will handle the sign conversion)")
            
            oxidation_potential = float(input("Oxidation half-reaction potential (V): "))
            reduction_potential = float(input("Reduction half-reaction potential (V): "))
            
            e_cell = reduction_potential - oxidation_potential
            is_favorable = e_cell > 0
            
            display_results_header()
            print(f"Reaction is {'Favorable' if is_favorable else 'Not favorable'}")
            print(f"Standard cell potential: {e_cell:.3f} V")
            
            steps = [
                f"1. Standard reduction potential for oxidation half-reaction: {oxidation_potential:.3f} V",
                f"2. Standard reduction potential for reduction half-reaction: {reduction_potential:.3f} V",
                f"3. Standard cell potential: E°cell = E°reduction - E°oxidation",
                f"4. E°cell = {reduction_potential:.3f} V - ({oxidation_potential:.3f} V) = {e_cell:.3f} V",
                f"5. Since E°cell {'>' if is_favorable else '≤'} 0, the reaction is {'favorable' if is_favorable else 'not favorable'}"
            ]
            
            display_steps(steps)
            
            # Additional information
            if is_favorable:
                print("\nThe positive cell potential indicates this reaction will occur spontaneously.")
            else:
                print("\nThe negative cell potential indicates this reaction will not occur spontaneously.")
                print("Energy input would be required to drive this reaction forward.")
                
        except ValueError:
            print("Error: Please enter valid numerical values for potentials.")
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_oxidation_state_analysis(self):
        """Handle standalone oxidation state analysis."""
        print("\n===== OXIDATION STATE ANALYSIS =====")
        
        try:
            print("Enter a chemical reaction to analyze oxidation state changes")
            print("Example: Fe + CuSO4 -> FeSO4 + Cu")
            reaction = input("Reaction: ")
            
            analysis = self.oxidation_analyzer.analyze_reaction(reaction)
            
            display_results_header()
            
            if analysis.is_redox:
                print("✓ This IS a redox reaction")
                
                print(f"\nOxidation state changes:")
                for item in analysis.oxidation_analysis:
                    print(f"  {item}")
                
                if analysis.oxidized_elements:
                    print(f"\n🔺 Elements oxidized (lose electrons):")
                    for elem in analysis.oxidized_elements:
                        print(f"  • {elem.element} in {elem.reactant_compound}: "
                              f"{elem.from_state:+d} → {elem.to_state:+d} (loses {elem.change} e⁻)")
                
                if analysis.reduced_elements:
                    print(f"\n🔻 Elements reduced (gain electrons):")
                    for elem in analysis.reduced_elements:
                        print(f"  • {elem.element} in {elem.reactant_compound}: "
                              f"{elem.from_state:+d} → {elem.to_state:+d} (gains {elem.change} e⁻)")
                
                # Electron balance check
                print(f"\nElectron balance check:")
                print(f"  Total electrons lost: {analysis.total_electrons_lost}")
                print(f"  Total electrons gained: {analysis.total_electrons_gained}")
                
                if analysis.is_balanced:
                    print("  ✓ Electrons are balanced")
                else:
                    print("  ⚠ Electrons are not balanced - check reaction coefficients")
                    
            else:
                print("✗ This is NOT a redox reaction")
                print("No changes in oxidation states were detected.")
                
                if analysis.oxidation_analysis:
                    print(f"\nOxidation state analysis:")
                    for item in analysis.oxidation_analysis:
                        print(f"  {item}")
                        
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_half_reaction_lookup(self):
        """Handle half-reaction lookup by element."""
        print("\n===== HALF-REACTION LOOKUP =====")
        print("This tool helps find standard reduction potentials for half-reactions")
        
        element = input("\nEnter element symbol (e.g., Cu, Fe, Zn, Mn, S): ")
        state = input("Enter state (optional - s, aq, g): ")
        
        half_reactions = get_half_reaction_by_element(element, state)
        
        display_results_header()
        if half_reactions:
            print(f"Found {len(half_reactions)} half-reactions for {element}:")
            for i, reaction in enumerate(half_reactions, 1):
                print(f"\n{i}. {reaction['half_reaction']}")
                print(f"   Standard Reduction Potential: {reaction['potential']:+.3f} V")
                
            print("\nTo use these in analysis:")
            print("1. Identify the relevant half-reactions for your redox process")
            print("2. Use option [2] in the menu to calculate favorability")
            print("   using the potentials shown above")
            
            # Add option to check a specific half-reaction
            print("\nWould you like to look up a specific half-reaction? (y/n)")
            choice = input("> ").strip().lower()
            
            if choice == 'y':
                print("\nEnter the half-reaction (use -> for arrow):")
                print("Example: Cu2+ + 2e- -> Cu")
                half_reaction = input("> ")
                potential = self.potential_calculator.find_standard_reduction_potential(half_reaction)
                
                if potential is not None:
                    print(f"\nStandard reduction potential for '{half_reaction}': {potential:+.3f} V")
                else:
                    print(f"\nCould not find standard reduction potential for '{half_reaction}'")
                    print("Try checking the format or spelling of the half-reaction.")
        else:
            print(f"No half-reactions found for {element}" + 
                  (f" in state {state}" if state else "") + ".")
            print("Try searching for the base element or check the spelling.")
    
    def _handle_common_redox_lookup(self):
        """Handle lookup of common redox reactions and pairs."""
        print("\n===== COMMON REDOX REACTIONS LOOKUP =====")
        
        print("\nLookup options:")
        print("[1] Browse all common redox reactions")
        print("[2] Browse all common redox pairs")
        print("[3] Search by specific name")
        print("[0] Back to main redox menu")
        
        choice = input("\nEnter choice (0-3): ").strip()
        
        if choice == "0":
            return
        elif choice == "1":
            self._browse_common_reactions()
        elif choice == "2":
            self._browse_common_pairs()
        elif choice == "3":
            self._search_by_name()
        else:
            print("Invalid choice.")
    
    def _browse_common_reactions(self):
        """Browse all common redox reactions."""
        print("\n===== COMMON REDOX REACTIONS =====")
        
        for i, (name, data) in enumerate(COMMON_REDOX_REACTIONS.items(), 1):
            print(f"\n{i}. {name}")
            print(f"   Reaction: {data.get('reaction', 'N/A')}")
            print(f"   Balanced: {data.get('balanced_equation', 'N/A')}")
            
        print(f"\nFound {len(COMMON_REDOX_REACTIONS)} common redox reactions.")
        print("Enter the number to see details, or 0 to go back:")
        
        try:
            choice = int(input("> "))
            if 1 <= choice <= len(COMMON_REDOX_REACTIONS):
                name = list(COMMON_REDOX_REACTIONS.keys())[choice - 1]
                self._show_reaction_details(name, COMMON_REDOX_REACTIONS[name])
        except (ValueError, IndexError):
            print("Invalid selection.")
    
    def _browse_common_pairs(self):
        """Browse all common redox pairs."""
        print("\n===== COMMON REDOX PAIRS =====")
        
        for i, (name, data) in enumerate(COMMON_REDOX_PAIRS.items(), 1):
            favorable_str = "✓ Favorable" if data.get('favorable', False) else "✗ Not favorable"
            print(f"\n{i}. {name} - {favorable_str}")
            print(f"   E°cell: {data.get('e_cell', 'N/A')} V")
            print(f"   Net reaction: {data.get('net_reaction', 'N/A')}")
            
        print(f"\nFound {len(COMMON_REDOX_PAIRS)} common redox pairs.")
        print("Enter the number to see details, or 0 to go back:")
        
        try:
            choice = int(input("> "))
            if 1 <= choice <= len(COMMON_REDOX_PAIRS):
                name = list(COMMON_REDOX_PAIRS.keys())[choice - 1]
                self._show_pair_details(name, COMMON_REDOX_PAIRS[name])
        except (ValueError, IndexError):
            print("Invalid selection.")
    
    def _search_by_name(self):
        """Search for specific reactions or pairs by name."""
        print("\nEnter search term (partial matches allowed):")
        search_term = input("> ").lower()
        
        # Search reactions
        matching_reactions = {name: data for name, data in COMMON_REDOX_REACTIONS.items() 
                            if search_term in name.lower()}
        
        # Search pairs
        matching_pairs = {name: data for name, data in COMMON_REDOX_PAIRS.items() 
                         if search_term in name.lower()}
        
        if matching_reactions:
            print(f"\nFound {len(matching_reactions)} matching reactions:")
            for name, data in matching_reactions.items():
                print(f"• {name}: {data.get('reaction', 'N/A')}")
        
        if matching_pairs:
            print(f"\nFound {len(matching_pairs)} matching pairs:")
            for name, data in matching_pairs.items():
                favorable_str = "✓" if data.get('favorable', False) else "✗"
                print(f"• {name} {favorable_str}: {data.get('net_reaction', 'N/A')}")
        
        if not matching_reactions and not matching_pairs:
            print(f"No reactions or pairs found matching '{search_term}'")
    
    def _show_reaction_details(self, name: str, data: dict):
        """Show detailed information about a reaction."""
        print(f"\n===== {name.upper()} =====")
        print(f"Reaction: {data.get('reaction', 'N/A')}")
        print(f"Balanced equation: {data.get('balanced_equation', 'N/A')}")
        print(f"Oxidation half: {data.get('oxidation_half', 'N/A')}")
        print(f"Reduction half: {data.get('reduction_half', 'N/A')}")
        
        if 'notes' in data:
            print(f"\nNotes: {data['notes']}")
    
    def _show_pair_details(self, name: str, data: dict):
        """Show detailed information about a redox pair."""
        print(f"\n===== {name.upper()} =====")
        print(f"Net reaction: {data.get('net_reaction', 'N/A')}")
        print(f"Oxidation: {data.get('oxidation', 'N/A')}")
        print(f"Reduction: {data.get('reduction', 'N/A')}")
        print(f"Standard cell potential: {data.get('e_cell', 'N/A')} V")
        print(f"Favorable: {'Yes' if data.get('favorable', False) else 'No'}")
        
        if 'notes' in data:
            print(f"\nNotes: {data['notes']}")
    
    def _handle_nernst_equation(self):
        """Handle Nernst equation calculations."""
        print("\n===== NERNST EQUATION CALCULATOR =====")
        print("Calculate cell potential under non-standard conditions")
        
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
            
            e_cell = self.potential_calculator.calculate_nernst_equation(
                e_cell_standard, n, concentration_ratio, temperature)
            
            display_results_header()
            print(f"Standard cell potential (E°cell): {e_cell_standard:+.3f} V")
            print(f"Cell potential under given conditions (Ecell): {e_cell:+.3f} V")
            
            # Show calculation steps
            rt_nf = (8.314 * temperature) / (n * 96485)
            log_term = 2.303 * __import__('math').log10(concentration_ratio)
            
            steps = [
                f"1. Standard cell potential (E°cell): {e_cell_standard:+.3f} V",
                f"2. Number of electrons transferred (n): {n}",
                f"3. Concentration ratio [products]/[reactants]: {concentration_ratio}",
                f"4. Temperature: {temperature-273.15:.1f}°C ({temperature:.2f} K)",
                f"5. Using the Nernst equation: Ecell = E°cell - ((RT)/(nF))·2.303·log10(Q)",
                f"   where R = 8.314 J/(mol·K), F = 96485 C/mol",
                f"6. (RT)/(nF) = (8.314 × {temperature:.2f})/({n} × 96485) = {rt_nf:.6f}",
                f"7. 2.303·log10({concentration_ratio}) = {log_term:.6f}",
                f"8. Ecell = {e_cell_standard:.3f} - ({rt_nf:.6f} × {log_term:.6f}) = {e_cell:.3f} V"
            ]
            
            display_steps(steps)
            
            # Provide interpretation
            if e_cell > 0:
                print("\n✓ The reaction is favorable under these conditions (Ecell > 0).")
            else:
                print("\n✗ The reaction is not favorable under these conditions (Ecell ≤ 0).")
            
            # Show effect of concentration
            if concentration_ratio > 1:
                print("Higher product concentrations decrease the driving force.")
            elif concentration_ratio < 1:
                print("Higher reactant concentrations increase the driving force.")
            else:
                print("Equal concentrations of products and reactants.")
            
        except ValueError:
            print("Error: Please enter valid numerical values.")
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_complex_reaction_patterns(self):
        """Handle complex reaction pattern recognition."""
        print("\n===== COMPLEX REACTION PATTERN RECOGNITION =====")
        print("This tool recognizes patterns in complex redox reactions")
        
        print("\nSupported patterns:")
        patterns = ComplexReactionHandler.COMPLEX_PATTERNS
        for i, pattern_name in enumerate(patterns.keys(), 1):
            description = pattern_name.replace("_", " ").title()
            print(f"{i}. {description}")
        
        print(f"\nEnter a complex redox reaction to analyze:")
        print("Examples:")
        print("• SO2 + MnO4- + H2O -> SO42- + Mn2+ + H+")
        print("• Cr2O7 2- + Fe2+ + H+ -> Cr3+ + Fe3+ + H2O")
        print("• 2H2O2 -> 2H2O + O2")
        
        reaction = input("\nReaction: ")
        
        complex_handler = ComplexReactionHandler()
        result = complex_handler.parse_complex_reaction(reaction)
        
        display_results_header()
        
        if not result["oxidation_half"].startswith("[Please"):
            print("✓ Pattern recognized!")
            print(f"Oxidation half-reaction: {result['oxidation_half']}")
            print(f"Reduction half-reaction: {result['reduction_half']}")
            print(f"Balanced equation: {result['balanced_equation']}")
            print(f"Oxidizing agent: {result['oxidizing_agent']}")
            print(f"Reducing agent: {result['reducing_agent']}")
            print(f"Electron transfer: {result['electron_transfer']}")
            
            if result.get('notes'):
                print(f"\nNotes: {result['notes']}")
        else:
            print("✗ Pattern not recognized")
            print("This complex reaction needs manual balancing using the half-reaction method.")
            print("Use option [1] for comprehensive analysis with balancing guidance.")
    
    def _handle_about_redox(self):
        """Display comprehensive information about redox reactions."""
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