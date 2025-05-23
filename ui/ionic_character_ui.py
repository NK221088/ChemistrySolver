"""
Terminal User Interface for Ionic Character Calculator
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user
from chemistry_solver.ionic_character import IonicCharacterCalculator

class IonicCharacterUI:
    """UI class for ionic character calculations."""
    
    def __init__(self):
        self.title = "IONIC CHARACTER CALCULATOR"
        self.calculator = IonicCharacterCalculator()
    
    def run(self):
        """Run the ionic character UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-7): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_single_compound()
            elif choice == "2":
                self._handle_compound_comparison()
            elif choice == "3":
                self._handle_order_compounds()
            elif choice == "4":
                self._handle_multiple_choice()
            elif choice == "5":
                self._handle_polyatomic_analysis()
            elif choice == "6":
                self._handle_polyatomic_multiple_choice()
            elif choice == "7":
                self._handle_theory_explanation()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the ionic character module menu."""
        menu = """
        [1] Analyze single compound
        [2] Compare two compounds
        [3] Order compounds by ionic character
        [4] Solve ionic character multiple choice
        [5] Analyze polyatomic ions in compound
        [6] Solve polyatomic ions multiple choice
        [7] Theory and explanation
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_single_compound(self):
        """Handle analysis of a single compound."""
        print("\n===== SINGLE COMPOUND ANALYSIS =====")
        
        try:
            compound = input("Enter compound formula (e.g., HCl, NaCl, CaF2): ").strip()
            
            result = self.calculator.get_electronegativity_difference(compound)
            
            display_results_header()
            if result['success']:
                print(f"Compound: {result['compound']}")
                print(f"Elements: {result['element1']} and {result['element2']}")
                print(f"Electronegativity of {result['element1']}: {result['electronegativity1']}")
                print(f"Electronegativity of {result['element2']}: {result['electronegativity2']}")
                print(f"Electronegativity difference: {result['electronegativity_difference']:.3f}")
                print(f"Bond type: {result['bond_type']}")
                print(f"Percent ionic character: {result['percent_ionic_character']:.1f}%")
            else:
                print(f"Error: {result['error']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_compound_comparison(self):
        """Handle comparison of two compounds."""
        print("\n===== COMPOUND COMPARISON =====")
        
        try:
            compound1 = input("Enter first compound formula: ").strip()
            compound2 = input("Enter second compound formula: ").strip()
            
            result1 = self.calculator.get_electronegativity_difference(compound1)
            result2 = self.calculator.get_electronegativity_difference(compound2)
            
            display_results_header()
            
            if result1['success'] and result2['success']:
                print(f"Compound 1: {result1['compound']}")
                print(f"  Electronegativity difference: {result1['electronegativity_difference']:.3f}")
                print(f"  Bond type: {result1['bond_type']}")
                print(f"  Percent ionic character: {result1['percent_ionic_character']:.1f}%")
                
                print(f"\nCompound 2: {result2['compound']}")
                print(f"  Electronegativity difference: {result2['electronegativity_difference']:.3f}")
                print(f"  Bond type: {result2['bond_type']}")
                print(f"  Percent ionic character: {result2['percent_ionic_character']:.1f}%")
                
                print("\nComparison:")
                if result1['electronegativity_difference'] > result2['electronegativity_difference']:
                    print(f"{result1['compound']} has MORE ionic character than {result2['compound']}")
                elif result1['electronegativity_difference'] < result2['electronegativity_difference']:
                    print(f"{result2['compound']} has MORE ionic character than {result1['compound']}")
                else:
                    print(f"Both compounds have the SAME ionic character")
            else:
                if not result1['success']:
                    print(f"Error with {compound1}: {result1['error']}")
                if not result2['success']:
                    print(f"Error with {compound2}: {result2['error']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_order_compounds(self):
        """Handle ordering compounds by ionic character."""
        print("\n===== ORDER COMPOUNDS BY IONIC CHARACTER =====")
        
        try:
            print("Enter compound formulas separated by commas (e.g., HCl, NaCl, CaF2):")
            compounds_input = input().strip()
            compounds = [c.strip() for c in compounds_input.split(',') if c.strip()]
            
            if len(compounds) < 2:
                print("Please enter at least 2 compounds.")
                return
            
            result = self.calculator.order_by_ionic_character(compounds)
            
            display_results_header()
            if result['success']:
                print("Compounds ordered by INCREASING ionic character:")
                for i, compound in enumerate(result['ordered_compounds'], 1):
                    print(f"{i}. {compound}")
                
                print("\nDetailed analysis:")
                for analysis in result['detailed_analysis']:
                    print(f"{analysis['compound']}: ΔEN = {analysis['electronegativity_difference']:.3f}, "
                          f"{analysis['percent_ionic_character']:.1f}% ionic")
                
                if result['failed_compounds']:
                    print("\nFailed to analyze:")
                    for failed in result['failed_compounds']:
                        print(f"- {failed.get('compound', 'Unknown')}: {failed.get('error', 'Unknown error')}")
            else:
                print(f"Error: {result['error']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_multiple_choice(self):
        """Handle multiple choice problem solving."""
        print("\n===== MULTIPLE CHOICE PROBLEM SOLVER =====")
        print("This will help you solve problems like:")
        print("'Identify the list arranged in increasing order of ionic character'")
        
        try:
            choices = []
            print("\nEnter each choice (compounds separated by commas).")
            print("Type 'none' for 'None of the above' option.")
            print("Type 'done' when finished entering choices.")
            
            choice_num = 1
            while True:
                choice_input = input(f"Choice {choice_num}: ").strip()
                
                if choice_input.lower() == 'done':
                    break
                elif choice_input.lower() in ['none', 'none of the above']:
                    choices.append('none of the above')
                else:
                    compounds = [c.strip() for c in choice_input.split(',') if c.strip()]
                    if compounds:
                        choices.append(compounds)
                
                choice_num += 1
            
            if not choices:
                print("No choices entered.")
                return
            
            result = self.calculator.solve_multiple_choice_problem(choices)
            
            display_results_header()
            print(f"ANSWER: Choice {result['answer']}")
            print(f"\nEXPLANATION:")
            print(result['explanation'])
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_polyatomic_analysis(self):
        """Handle analysis of polyatomic ions in a compound."""
        print("\n===== POLYATOMIC ION ANALYSIS =====")
        
        try:
            compound = input("Enter compound formula (e.g., NH4ClO3, K2Cr2O7): ").strip()
            
            result = self.calculator.parse_ionic_compound(compound)
            
            display_results_header()
            print(f"Compound: {result['compound']}")
            
            print(f"\nCations found:")
            if result['cations']:
                for cation in result['cations']:
                    ion_type = "POLYATOMIC" if cation['type'] == 'polyatomic' else "Monoatomic"
                    print(f"  - {cation['name']} ({cation['formula']}) - {ion_type}")
                    print(f"    Charge: +{cation['charge']}, Count: {cation['count']}")
            else:
                print("  No cations identified")
            
            print(f"\nAnions found:")
            if result['anions']:
                for anion in result['anions']:
                    ion_type = "POLYATOMIC" if anion['type'] == 'polyatomic' else "Monoatomic"
                    print(f"  - {anion['name']} ({anion['formula']}) - {ion_type}")
                    print(f"    Charge: {anion['charge']}, Count: {anion['count']}")
            else:
                print("  No anions identified")
            
            print(f"\nSummary:")
            print(f"  Has polyatomic cation: {'YES' if result['has_polyatomic_cation'] else 'NO'}")
            print(f"  Has polyatomic anion: {'YES' if result['has_polyatomic_anion'] else 'NO'}")
            print(f"  Has BOTH polyatomic cation AND anion: {'YES' if result['has_polyatomic_cation'] and result['has_polyatomic_anion'] else 'NO'}")
            
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_polyatomic_multiple_choice(self):
        """Handle polyatomic ion multiple choice problems."""
        print("\n===== POLYATOMIC ION MULTIPLE CHOICE SOLVER =====")
        print("This solves problems like:")
        print("'Which compound contains a polyatomic cation AND polyatomic anion?'")
        
        try:
            choices = []
            print("\nEnter each choice (one compound per choice).")
            print("Type 'none' for 'None of the above' option.")
            print("Type 'done' when finished entering choices.")
            
            choice_num = 1
            while True:
                choice_input = input(f"Choice {choice_num}: ").strip()
                
                if choice_input.lower() == 'done':
                    break
                elif choice_input.lower() in ['none', 'none of the above']:
                    choices.append('none of the above')
                else:
                    choices.append(choice_input)
                
                choice_num += 1
            
            if not choices:
                print("No choices entered.")
                return
            
            result = self.calculator.solve_polyatomic_multiple_choice(choices)
            
            display_results_header()
            print(f"ANSWER: Choice {result['answer']}")
            print(f"\nEXPLANATION:")
            print(result['explanation'])
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_theory_explanation(self):
        """Display theory and explanation about ionic character."""
        print("\n===== IONIC CHARACTER THEORY =====")
        
        theory = """
IONIC CHARACTER FUNDAMENTALS:

1. ELECTRONEGATIVITY:
   - Measure of an atom's ability to attract electrons in a bond
   - Higher electronegativity = stronger pull on electrons
   - Fluorine has highest electronegativity (3.98)

2. ELECTRONEGATIVITY DIFFERENCE (ΔEN):
   - ΔEN = |EN₁ - EN₂|
   - Determines bond type:
     * ΔEN < 0.5: Nonpolar covalent
     * 0.5 ≤ ΔEN < 1.7: Polar covalent
     * ΔEN ≥ 1.7: Ionic

3. IONIC CHARACTER:
   - Measure of how "ionic" vs "covalent" a bond is
   - Higher ΔEN = More ionic character
   - Can be calculated as percentage using Pauling's equation

4. ORDERING BY IONIC CHARACTER:
   - Calculate ΔEN for each compound
   - Arrange from smallest to largest ΔEN
   - Smallest ΔEN = Most covalent = Least ionic
   - Largest ΔEN = Most ionic = Least covalent

5. COMMON PATTERNS:
   - Metal + Nonmetal = Usually ionic
   - Nonmetal + Nonmetal = Usually covalent
   - Large ΔEN between periods = More ionic
   - Similar elements = More covalent

EXAMPLE:
   HF: ΔEN = |2.20 - 3.98| = 1.78 (Ionic)
   HCl: ΔEN = |2.20 - 3.16| = 0.96 (Polar covalent)
   Order: HCl < HF (increasing ionic character)

6. POLYATOMIC IONS:
   - Polyatomic ions are groups of atoms that carry a charge
   - Common polyatomic cations: NH4+ (ammonium), H3O+ (hydronium)
   - Common polyatomic anions: OH- (hydroxide), NO3- (nitrate), 
     SO4²- (sulfate), CO3²- (carbonate), ClO3- (chlorate)
   - Compounds can have polyatomic cations, polyatomic anions, or both
   
   Examples:
   - NH4Cl: polyatomic cation (NH4+) + monoatomic anion (Cl-)
   - NaNO3: monoatomic cation (Na+) + polyatomic anion (NO3-)
   - NH4ClO3: polyatomic cation (NH4+) + polyatomic anion (ClO3-)
        """
        
        print(theory)