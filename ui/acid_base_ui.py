"""
Terminal User Interface for Acid-Base Analysis
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user
from chemistry_solver.acid_base import (
    identify_acid_base, 
    analyze_compound_list, 
    AcidBaseEquilibrium,
    analyze_weak_acid_question
)
# Import the new solubility equilibrium functions
from chemistry_solver.solubility_pH_extension import (
    solve_hydroxide_salt_pH,
    solve_ksp_from_pH,
    analyze_salt_solubility_questions
)

class AcidBaseUI:
    """UI class for acid-base chemistry analysis."""
    
    def __init__(self):
        self.title = "ACID-BASE ANALYZER"
        self.equilibrium_solver = AcidBaseEquilibrium()
    
    def run(self):
        """Run the acid-base UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-8): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_single_compound()
            elif choice == "2":
                self._handle_multiple_compounds()
            elif choice == "3":
                self._handle_weak_acid_equilibrium()
            elif choice == "4":
                self._handle_weak_base_equilibrium()
            elif choice == "5":
                self._handle_acid_statement_verification()
            elif choice == "6":
                self._handle_acid_base_info()
            elif choice == "7":
                self._handle_salt_pH_calculation()
            elif choice == "8":
                self._handle_ksp_from_pH()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the acid-base module menu."""
        menu = """
        [1] Identify acid/base for single compound
        [2] Analyze multiple compounds
        [3] Calculate weak acid equilibrium
        [4] Calculate weak base equilibrium
        [5] Verify statements about weak acid solutions
        [6] Learn about acid-base theory
        [7] Calculate pH from solubility product (Ksp)
        [8] Calculate solubility product (Ksp) from pH
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
    
    def _handle_weak_acid_equilibrium(self):
        """Handle weak acid equilibrium calculations."""
        print("\n===== WEAK ACID EQUILIBRIUM CALCULATOR =====")
        
        try:
            formula = input("Enter acid formula (e.g., CH3COOH): ")
            
            # Display Ka value if available
            ka_value = self.equilibrium_solver.get_acid_ka(formula)
            if ka_value is not None:
                print(f"Ka value found in database: {ka_value:.2e}")
                use_default = input("Use this Ka value? (y/n): ").strip().lower() == 'y'
            else:
                print(f"No Ka value found for {formula} in database.")
                use_default = False
            
            if not use_default:
                try:
                    ka_input = input("Enter Ka value (scientific notation OK, e.g., 1.8e-5): ")
                    ka_value = float(ka_input)
                except ValueError:
                    print("Invalid Ka value. Using 1.0e-5 as default.")
                    ka_value = 1.0e-5
            
            concentration_input = input("Enter initial acid concentration (mol/L): ")
            try:
                concentration = float(concentration_input)
            except ValueError:
                print("Invalid concentration. Using 0.1 mol/L as default.")
                concentration = 0.1
            
            result = self.equilibrium_solver.solve_weak_acid_equilibrium(
                formula, concentration, ka_value
            )
            
            if "error" in result:
                print(f"Error: {result['error']}")
                return
            
            display_results_header()
            print(f"Equilibrium Analysis for {formula} at {concentration} mol/L")
            print("-" * 50)
            print(f"Ka value: {result['ka']:.2e}")
            print(f"pH: {result['ph']:.4f}")
            print(f"[H+]: {result['h_concentration']:.2e} mol/L")
            print(f"[{formula[:-1]}⁻]: {result['conjugate_base_concentration']:.2e} mol/L")
            print(f"[{formula}]: {result['undissociated_concentration']:.2e} mol/L")
            print(f"Percent dissociation: {result['percent_dissociation']:.2f}%")
            
            if result['is_approximation_valid']:
                print("\nNote: The approximation x << C₀ is valid for this calculation.")
            else:
                print("\nNote: The approximation x << C₀ is NOT valid. Quadratic formula was used.")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_weak_base_equilibrium(self):
        """Handle weak base equilibrium calculations."""
        print("\n===== WEAK BASE EQUILIBRIUM CALCULATOR =====")
        
        try:
            formula = input("Enter base formula (e.g., NH3): ")
            
            # Display Kb value if available
            kb_value = self.equilibrium_solver.get_base_kb(formula)
            if kb_value is not None:
                print(f"Kb value found in database: {kb_value:.2e}")
                use_default = input("Use this Kb value? (y/n): ").strip().lower() == 'y'
            else:
                print(f"No Kb value found for {formula} in database.")
                use_default = False
            
            if not use_default:
                try:
                    kb_input = input("Enter Kb value (scientific notation OK, e.g., 1.8e-5): ")
                    kb_value = float(kb_input)
                except ValueError:
                    print("Invalid Kb value. Using 1.0e-5 as default.")
                    kb_value = 1.0e-5
            
            concentration_input = input("Enter initial base concentration (mol/L): ")
            try:
                concentration = float(concentration_input)
            except ValueError:
                print("Invalid concentration. Using 0.1 mol/L as default.")
                concentration = 0.1
            
            result = self.equilibrium_solver.solve_weak_base_equilibrium(
                formula, concentration, kb_value
            )
            
            if "error" in result:
                print(f"Error: {result['error']}")
                return
            
            display_results_header()
            print(f"Equilibrium Analysis for {formula} at {concentration} mol/L")
            print("-" * 50)
            print(f"Kb value: {result['kb']:.2e}")
            print(f"pH: {result['ph']:.4f}")
            print(f"pOH: {result['poh']:.4f}")
            print(f"[OH⁻]: {result['oh_concentration']:.2e} mol/L")
            print(f"[{formula}H⁺]: {result['conjugate_acid_concentration']:.2e} mol/L")
            print(f"[{formula}]: {result['undissociated_concentration']:.2e} mol/L")
            print(f"Percent ionization: {result['percent_ionization']:.2f}%")
            
            if result['is_approximation_valid']:
                print("\nNote: The approximation x << C₀ is valid for this calculation.")
            else:
                print("\nNote: The approximation x << C₀ is NOT valid. Quadratic formula was used.")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_acid_statement_verification(self):
        """Handle verification of statements about weak acid solutions."""
        print("\n===== WEAK ACID STATEMENT VERIFIER =====")
        print("This tool verifies common statements about weak acid solutions")
        
        try:
            formula = input("Enter acid formula (default: HX): ").strip()
            if not formula:
                formula = "HX"
            
            # Get Ka value
            try:
                ka_input = input("Enter Ka value (default: 1.0e-5): ").strip()
                ka_value = float(ka_input) if ka_input else 1.0e-5
            except ValueError:
                print("Invalid Ka value. Using 1.0e-5 as default.")
                ka_value = 1.0e-5
            
            # Get concentration
            try:
                conc_input = input("Enter concentration in mol/L (default: 1.0): ").strip()
                concentration = float(conc_input) if conc_input else 1.0
            except ValueError:
                print("Invalid concentration. Using 1.0 mol/L as default.")
                concentration = 1.0
            
            result = self.equilibrium_solver.verify_weak_acid_statements(
                formula, concentration, ka_value
            )
            
            if "error" in result:
                print(f"Error: {result['error']}")
                return
            
            display_results_header()
            print(f"Statement Analysis for {formula} (Ka = {ka_value:.2e}, [{formula}] = {concentration} M)")
            print("-" * 50)
            
            # Display equilibrium data
            eq_data = result["equilibrium_data"]
            print(f"pH: {eq_data['ph']:.4f}")
            print(f"[H⁺]: {eq_data['h_concentration']:.2e} mol/L")
            print(f"Percent dissociation: {eq_data['percent_dissociation']:.2f}%")
            print("-" * 50)
            
            # Display statement verification
            print("Statement Analysis:")
            statements = result["statements"]
            
            print(f"pH = 0: {statements['ph_equals_0']}")
            print(f"[X⁻] = initial concentration: {statements['conjugate_base_equals_initial']}")
            print(f"[HX] > [H⁺]: {statements['undissociated_greater_than_h']}")
            print(f"[H⁺] = initial concentration: {statements['h_equals_initial']}")
            print(f"pH = -log(initial concentration): {statements['ph_equals_negative_log_initial']}")
            print(f"100% dissociation: {statements['percent_dissociation_100']}")
            print(f"Less than 5% dissociation: {statements['percent_dissociation_less_than_5']}")
            
            print("-" * 50)
            print("Explanation:")
            print(result["explanation"])
                
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

Equilibrium Concepts:
- Ka: Acid dissociation constant, higher Ka = stronger acid
- Kb: Base dissociation constant, higher Kb = stronger base
- Ka × Kb = Kw = 1.0 × 10^-14 at 25°C (for conjugate acid-base pairs)
- pH = -log[H+], pOH = -log[OH-], pH + pOH = 14 at 25°C

Solubility Product (Ksp):
- Describes the equilibrium between a solid salt and its constituent ions in solution
- For a salt MₓXᵧ: MₓXᵧ(s) ⇌ xMⁿ⁺(aq) + yXᵐ⁻(aq)
- Ksp = [Mⁿ⁺]ˣ[Xᵐ⁻]ʸ
- Used to calculate solubility of sparingly soluble salts
- Smaller Ksp values indicate lower solubility
        """
        
        print(info)

    def _handle_salt_pH_calculation(self):
        """Handle pH calculation from solubility product (Ksp)."""
        print("\n===== SOLUBILITY PRODUCT pH CALCULATOR =====")
        print("This tool calculates pH of a solution saturated with a sparingly soluble salt")
        
        try:
            formula = input("Enter salt formula (e.g., Mg(OH)2, Ca(OH)2): ").strip()
            if not formula:
                print("Error: Formula cannot be empty.")
                return
            
            # Check if formula contains hydroxide
            if "OH" not in formula:
                print("Note: This calculation currently works best for metal hydroxides.")
                print("Continuing with calculation, but results may not be accurate for non-hydroxide salts.")
            
            # Get Ksp value
            try:
                ksp_input = input("Enter Ksp value (scientific notation OK, e.g., 5.61e-12): ").strip()
                ksp_value = float(ksp_input)
            except ValueError:
                print("Invalid Ksp value. Using 1.0e-10 as default.")
                ksp_value = 1.0e-10
            
            # Calculate pH using the imported function
            result = solve_hydroxide_salt_pH(formula, ksp_value)
            
            if "error" in result:
                print(f"Error: {result['error']}")
                return
            
            display_results_header()
            print(f"pH Calculation for Saturated {formula} Solution")
            print("-" * 50)
            print(f"Salt formula: {formula}")
            print(f"Ksp value: {ksp_value:.2e}")
            print(f"Metal: {result['metal']}")
            print(f"Hydroxide groups: {result['hydroxide_groups']}")
            print("-" * 50)
            print(f"Solubility: {result['solubility']:.2e} mol/L")
            print(f"[{result['metal']}²⁺]: {result['metal_ion_concentration']:.2e} mol/L")
            print(f"[OH⁻]: {result['hydroxide_concentration']:.2e} mol/L")
            print(f"pOH: {result['poh']:.4f}")
            print(f"pH: {result['ph']:.4f}")
            
            # If it's the specific problem from the example, show the answer options
            if formula.lower() == "mg(oh)2" and abs(ksp_value - 5.61e-12) < 1e-13:
                print("\nFor the multiple-choice question:")
                options = ["~8.4", "~9.4", "~10.4", "~11.4", "~12.4"]
                closest_option = None
                smallest_diff = float('inf')
                
                print("Options analysis:")
                for i, option_str in enumerate(options):
                    # Extract numeric value
                    import re
                    ph_match = re.search(r"~?\s*(\d+\.\d+)", option_str)
                    if ph_match:
                        option_ph = float(ph_match.group(1))
                        diff = abs(option_ph - result['ph'])
                        
                        if diff < smallest_diff:
                            smallest_diff = diff
                            closest_option = i
                        
                        print(f"  {option_str}: {'✓' if diff < 0.1 else '✗'} (diff: {diff:.4f})")
                
                print(f"\nThe calculated pH is {result['ph']:.4f}, closest to option: {options[closest_option]}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_ksp_from_pH(self):
        """Handle calculation of solubility product (Ksp) from pH."""
        print("\n===== SOLUBILITY PRODUCT FROM pH =====")
        print("This tool calculates Ksp from the pH of a saturated salt solution")
        
        try:
            formula = input("Enter salt formula (e.g., Mg(OH)2, Ca(OH)2): ").strip()
            if not formula:
                print("Error: Formula cannot be empty.")
                return
            
            # Check if formula contains hydroxide
            if "OH" not in formula:
                print("Note: This calculation currently works best for metal hydroxides.")
                print("Continuing with calculation, but results may not be accurate for non-hydroxide salts.")
            
            # Get pH value
            try:
                ph_input = input("Enter pH value of the saturated solution: ").strip()
                ph_value = float(ph_input)
            except ValueError:
                print("Invalid pH value. Using 10.0 as default.")
                ph_value = 10.0
            
            # Calculate Ksp using the imported function
            result = solve_ksp_from_pH(formula, ph_value)
            
            if "error" in result:
                print(f"Error: {result['error']}")
                return
            
            display_results_header()
            print(f"Ksp Calculation for {formula} from pH")
            print("-" * 50)
            print(f"Salt formula: {formula}")
            print(f"Input pH: {ph_value:.4f}")
            print(f"Metal: {result['metal']}")
            print(f"Hydroxide groups: {result['hydroxide_groups']}")
            print("-" * 50)
            print(f"Calculated pOH: {result['poh']:.4f}")
            print(f"[OH⁻]: {result['hydroxide_concentration']:.2e} mol/L")
            print(f"Solubility: {result['solubility']:.2e} mol/L")
            print(f"[{result['metal']}²⁺]: {result['metal_ion_concentration']:.2e} mol/L")
            print(f"Calculated Ksp: {result['ksp']:.2e}")
                
        except Exception as e:
            print(f"Error: {str(e)}")