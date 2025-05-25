"""
Terminal User Interface for Acid-Base Analysis
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user
from chemistry_solver.acid_base import identify_acid_base, AcidBaseEquilibrium

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
            choice = input("\nEnter choice (0-7): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_single_compound()
            elif choice == "2":
                self._handle_multiple_compounds()
            elif choice == "3":
                self._handle_general_ph_solver()
            elif choice == "4":
                self._handle_weak_acid_equilibrium()
            elif choice == "5":
                self._handle_weak_base_equilibrium()
            elif choice == "6":
                self._handle_weak_acid_analysis()
            elif choice == "7":
                self._handle_acid_base_info()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the acid-base module menu."""
        menu = """
        [1] Identify acid/base for single compound
        [2] Analyze multiple compounds
        [3] General pH solver (auto-detect acid/base type)
        [4] Calculate weak acid equilibrium
        [5] Calculate weak base equilibrium
        [6] Analyze weak acid multiple choice questions
        [7] Learn about acid-base theory
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_single_compound(self):
        """Handle acid/base identification for a single compound."""
        print("\n===== ACID/BASE IDENTIFIER =====")
        
        try:
            compound = input("Enter chemical formula: ").strip()
            if not compound:
                print("Error: Formula cannot be empty.")
                return
            
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
            compounds_input = input("Enter chemical formulas (comma-separated): ").strip()
            if not compounds_input:
                print("Error: Input cannot be empty.")
                return
                
            compounds = [comp.strip() for comp in compounds_input.split(',')]
            
            display_results_header()
            print(f"Analyzing compounds: {', '.join(compounds)}")
            print("-" * 50)
            
            acids = []
            bases = []
            neutral = []
            unknown = []
            
            for compound in compounds:
                result = identify_acid_base(compound)
                print(f"Compound: {compound}")
                print(f"Classification: {result['classification']}")
                print(f"Explanation: {result['explanation']}")
                print("-" * 50)
                
                # Categorize results
                classification = result['classification'].lower()
                if 'acid' in classification:
                    acids.append(compound)
                elif 'base' in classification:
                    bases.append(compound)
                elif 'neutral' in classification:
                    neutral.append(compound)
                else:
                    unknown.append(compound)
            
            # Summary
            print("SUMMARY:")
            if acids:
                print(f"Acids: {', '.join(acids)}")
            if bases:
                print(f"Bases: {', '.join(bases)}")
            if neutral:
                print(f"Neutral: {', '.join(neutral)}")
            if unknown:
                print(f"Unknown: {', '.join(unknown)}")
            
            all_acids = len(acids) == len(compounds) and len(compounds) > 0
            print(f"\nAll compounds are acids: {'Yes' if all_acids else 'No'}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_general_ph_solver(self):
        """Handle general pH solving that auto-detects acid/base type."""
        print("\n===== GENERAL pH SOLVER =====")
        print("This tool automatically detects if your compound is a strong/weak acid/base")
        
        try:
            formula = input("Enter chemical formula: ").strip()
            if not formula:
                print("Error: Formula cannot be empty.")
                return
            
            concentration_input = input("Enter concentration (mol/L): ").strip()
            try:
                concentration = float(concentration_input)
                if concentration <= 0:
                    print("Error: Concentration must be positive.")
                    return
            except ValueError:
                print("Invalid concentration. Using 0.1 mol/L as default.")
                concentration = 0.1
            
            # Check if user wants to provide Ka/Kb values
            compound_info = identify_acid_base(formula)
            print(f"\nCompound identified as: {compound_info['classification']}")
            
            ka_value = None
            kb_value = None
            
            # For weak acids, check if we have Ka or ask for it
            if 'weak' in compound_info['classification'].lower() and 'acid' in compound_info['classification'].lower():
                ka_value = self.equilibrium_solver.get_acid_ka(formula)
                if ka_value is not None:
                    print(f"Found Ka value in database: {ka_value:.2e}")
                    use_default = input("Use this Ka value? (y/n): ").strip().lower() == 'y'
                    if not use_default:
                        ka_value = None
                
                if ka_value is None:
                    ka_input = input("Enter Ka value (or press Enter to skip): ").strip()
                    if ka_input:
                        try:
                            ka_value = float(ka_input)
                        except ValueError:
                            print("Invalid Ka value entered.")
            
            # For weak bases, check if we have Kb or ask for it
            elif 'weak' in compound_info['classification'].lower() and 'base' in compound_info['classification'].lower():
                kb_value = self.equilibrium_solver.get_base_kb(formula)
                if kb_value is not None:
                    print(f"Found Kb value in database: {kb_value:.2e}")
                    use_default = input("Use this Kb value? (y/n): ").strip().lower() == 'y'
                    if not use_default:
                        kb_value = None
                
                if kb_value is None:
                    kb_input = input("Enter Kb value (or press Enter to skip): ").strip()
                    if kb_input:
                        try:
                            kb_value = float(kb_input)
                        except ValueError:
                            print("Invalid Kb value entered.")
            
            # Solve the problem
            result = self.equilibrium_solver.solve_ph_problem(
                formula, concentration, ka=ka_value, kb=kb_value
            )
            
            if "error" in result:
                print(f"Error: {result['error']}")
                return
            
            display_results_header()
            self._display_equilibrium_results(result, formula, concentration)
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_weak_acid_equilibrium(self):
        """Handle weak acid equilibrium calculations."""
        print("\n===== WEAK ACID EQUILIBRIUM CALCULATOR =====")
        
        try:
            formula = input("Enter acid formula (e.g., CH3COOH): ").strip()
            if not formula:
                print("Error: Formula cannot be empty.")
                return
            
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
                    ka_input = input("Enter Ka value (scientific notation OK, e.g., 1.8e-5): ").strip()
                    ka_value = float(ka_input)
                except ValueError:
                    print("Invalid Ka value. Using 1.0e-5 as default.")
                    ka_value = 1.0e-5
            
            concentration_input = input("Enter initial acid concentration (mol/L): ").strip()
            try:
                concentration = float(concentration_input)
                if concentration <= 0:
                    print("Error: Concentration must be positive.")
                    return
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
            print(f"pOH: {result['poh']:.4f}")
            print(f"[H+]: {result['h_concentration']:.2e} mol/L")
            print(f"[OH-]: {result['oh_concentration']:.2e} mol/L")
            print(f"[A⁻]: {result['conjugate_base_concentration']:.2e} mol/L")
            print(f"[HA]: {result['undissociated_concentration']:.2e} mol/L")
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
            formula = input("Enter base formula (e.g., NH3): ").strip()
            if not formula:
                print("Error: Formula cannot be empty.")
                return
            
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
                    kb_input = input("Enter Kb value (scientific notation OK, e.g., 1.8e-5): ").strip()
                    kb_value = float(kb_input)
                except ValueError:
                    print("Invalid Kb value. Using 1.0e-5 as default.")
                    kb_value = 1.0e-5
            
            concentration_input = input("Enter initial base concentration (mol/L): ").strip()
            try:
                concentration = float(concentration_input)
                if concentration <= 0:
                    print("Error: Concentration must be positive.")
                    return
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
            print(f"[H+]: {result['h_concentration']:.2e} mol/L")
            print(f"[BH⁺]: {result['conjugate_acid_concentration']:.2e} mol/L")
            print(f"[B]: {result['undissociated_concentration']:.2e} mol/L")
            print(f"Percent ionization: {result['percent_ionization']:.2f}%")
            
            if result['is_approximation_valid']:
                print("\nNote: The approximation x << C₀ is valid for this calculation.")
            else:
                print("\nNote: The approximation x << C₀ is NOT valid. Quadratic formula was used.")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_weak_acid_analysis(self):
        """Handle analysis of weak acid multiple choice questions."""
        print("\n===== WEAK ACID MULTIPLE CHOICE ANALYZER =====")
        print("This tool analyzes multiple choice questions about weak acid solutions")
        
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
            
            # Get options
            print("\nEnter the multiple choice options (press Enter on empty line to finish):")
            options = []
            option_num = 1
            while True:
                option = input(f"Option {option_num}: ").strip()
                if not option:
                    break
                options.append(option)
                option_num += 1
            
            if not options:
                # Use default options for demonstration
                options = [
                    "pH = 0",
                    "[X-] = 1 M", 
                    "[HX] > [H+]",
                    "[H+] = 1 M"
                ]
                print("Using default options for demonstration.")
            
            result = self.equilibrium_solver.analyze_weak_acid_question(
                options, concentration, formula, ka_value
            )
            
            if "error" in result:
                print(f"Error: {result['error']}")
                return
            
            display_results_header()
            print(f"Multiple Choice Analysis for {formula} (Ka = {ka_value:.2e}, [{formula}] = {concentration} M)")
            print("-" * 50)
            
            # Display equilibrium data
            eq_data = result["equilibrium_data"]
            print("Equilibrium Concentrations:")
            print(f"pH: {eq_data['ph']:.4f}")
            print(f"[H⁺]: {eq_data['h_concentration']:.2e} mol/L")
            print(f"[A⁻]: {eq_data['conjugate_base_concentration']:.2e} mol/L")
            print(f"[HA]: {eq_data['undissociated_concentration']:.2e} mol/L")
            print(f"Percent dissociation: {eq_data['percent_dissociation']:.2f}%")
            print("-" * 50)
            
            # Display options analysis
            print("Options Analysis:")
            for i, (option_key, analysis) in enumerate(result["options_analysis"].items()):
                status = "✓ CORRECT" if analysis["is_correct"] else "✗ INCORRECT"
                print(f"{i+1}. {analysis['statement']} - {status}")
                if analysis["explanation"]:
                    print(f"   Explanation: {analysis['explanation']}")
            
            print("-" * 50)
            if result["correct_options"]:
                print(f"Correct option(s): {', '.join(map(str, result['correct_options']))}")
            else:
                print("No clearly correct options found among the given choices.")
            
            print(f"\nSummary: {result['summary']}")
                
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

Common Strong Acids:
- HCl, HBr, HI, HNO3, H2SO4, HClO4, HClO3

Common Strong Bases:
- Group 1 hydroxides: NaOH, KOH, LiOH, RbOH, CsOH
- Some Group 2 hydroxides: Ca(OH)2, Sr(OH)2, Ba(OH)2

Equilibrium Concepts:
- Ka: Acid dissociation constant, higher Ka = stronger acid
- Kb: Base dissociation constant, higher Kb = stronger base
- Ka × Kb = Kw = 1.0 × 10^-14 at 25°C (for conjugate acid-base pairs)
- pH = -log[H+], pOH = -log[OH-], pH + pOH = 14 at 25°C

For weak acids: HA ⇌ H⁺ + A⁻
- Ka = [H⁺][A⁻]/[HA]
- Use quadratic formula when approximation x << C₀ is not valid

For weak bases: B + H₂O ⇌ BH⁺ + OH⁻
- Kb = [BH⁺][OH⁻]/[B]
- Similar mathematical treatment as weak acids
        """
        
        print(info)
    
    def _display_equilibrium_results(self, result, formula, concentration):
        """Helper method to display equilibrium calculation results."""
        print(f"Equilibrium Analysis for {formula} at {concentration} mol/L")
        print("-" * 50)
        
        # Display based on the type of calculation
        if result.get("acid_type") == "strong":
            print("Strong Acid Calculation:")
            print(f"pH: {result['ph']:.4f}")
            print(f"pOH: {result['poh']:.4f}")
            print(f"[H+]: {result['h_concentration']:.2e} mol/L")
            print(f"[OH-]: {result['oh_concentration']:.2e} mol/L")
            print("Complete dissociation assumed")
            
        elif result.get("acid_type") == "weak":
            print("Weak Acid Calculation:")
            print(f"Ka value: {result['ka']:.2e}")
            print(f"pH: {result['ph']:.4f}")
            print(f"pOH: {result['poh']:.4f}")
            print(f"[H+]: {result['h_concentration']:.2e} mol/L")
            print(f"[OH-]: {result['oh_concentration']:.2e} mol/L")
            print(f"[A⁻]: {result['conjugate_base_concentration']:.2e} mol/L")
            print(f"[HA]: {result['undissociated_concentration']:.2e} mol/L")
            print(f"Percent dissociation: {result['percent_dissociation']:.2f}%")
            
            if result['is_approximation_valid']:
                print("Note: The approximation x << C₀ is valid for this calculation.")
            else:
                print("Note: The approximation x << C₀ is NOT valid. Quadratic formula was used.")
                
        elif result.get("base_type") == "strong":
            print("Strong Base Calculation:")
            print(f"pH: {result['ph']:.4f}")
            print(f"pOH: {result['poh']:.4f}")
            print(f"[H+]: {result['h_concentration']:.2e} mol/L")
            print(f"[OH-]: {result['oh_concentration']:.2e} mol/L")
            print("Complete dissociation assumed")
            
        elif result.get("base_type") == "weak":
            print("Weak Base Calculation:")
            print(f"Kb value: {result['kb']:.2e}")
            print(f"pH: {result['ph']:.4f}")
            print(f"pOH: {result['poh']:.4f}")
            print(f"[H+]: {result['h_concentration']:.2e} mol/L")
            print(f"[OH-]: {result['oh_concentration']:.2e} mol/L")
            print(f"[BH⁺]: {result['conjugate_acid_concentration']:.2e} mol/L")
            print(f"[B]: {result['undissociated_concentration']:.2e} mol/L")
            print(f"Percent ionization: {result['percent_ionization']:.2f}%")
            
            if result['is_approximation_valid']:
                print("Note: The approximation x << C₀ is valid for this calculation.")
            else:
                print("Note: The approximation x << C₀ is NOT valid. Quadratic formula was used.")
        
        else:
            print("General calculation results:")
            for key, value in result.items():
                if key not in ["error", "formula", "initial_concentration"]:
                    print(f"{key}: {value}")