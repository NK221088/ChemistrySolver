"""
Terminal User Interface for Insoluble Salts & Qualitative Analysis
"""
from ui.terminal_ui import display_title, display_results_header, display_steps, wait_for_user
from chemistry_solver.insoluble_salts import (
    solve_qualitative_analysis_problem,
    analyze_specific_scenario,
    calculate_solubility_from_ksp,
    solve_ksp_problem,
    get_available_cations,
    get_available_reagents,
    get_available_compounds
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
            choice = input("\nEnter choice (0-6): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_qualitative_analysis()
            elif choice == "2":
                self._handle_ksp_calculations()
            elif choice == "3":
                self._handle_solubility_from_ksp()
            elif choice == "4":
                self._handle_predefined_scenarios()
            elif choice == "5":
                self._handle_compound_database()
            elif choice == "6":
                self._handle_solubility_rules()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the insoluble salts module menu."""
        menu = """
        [1] Qualitative Analysis Problem Solver
        [2] General Ksp Problem Solver
        [3] Calculate Solubility from Ksp
        [4] Analyze Predefined Scenarios
        [5] View Compound Database
        [6] Solubility Rules Reference
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
            precip_input = input("Enter reagents that cause precipitation (comma-separated, or press enter for none): ")
            if precip_input.strip():
                precipitates_with = [r.strip() for r in precip_input.split(",")]
            else:
                precipitates_with = []
            
            # Get reagents that DON'T cause precipitation
            no_precip_input = input("Enter reagents that DON'T cause precipitation (comma-separated, or press enter for none): ")
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
    
    def _handle_ksp_calculations(self):
        """Handle general Ksp problem solving."""
        print("\n===== GENERAL KSP PROBLEM SOLVER =====")
        
        try:
            print("\nAvailable problem types:")
            print("1. Calculate solubility from Ksp")
            print("2. Common ion effect (coming soon)")
            print("3. Precipitation calculations (coming soon)")
            
            prob_type = input("\nSelect problem type (1-3): ").strip()
            
            if prob_type == "1":
                self._handle_solubility_from_ksp()
            else:
                print("This problem type is not yet implemented.")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_solubility_from_ksp(self):
        """Handle solubility calculations from Ksp values."""
        print("\n===== CALCULATE SOLUBILITY FROM KSP =====")
        
        try:
            # Show available compounds
            print("\nAvailable compounds with Ksp data:")
            available_compounds = get_available_compounds()
            for i, compound in enumerate(available_compounds, 1):
                print(f"  {i:2d}. {compound}")
            
            # Get compound choice
            print("\nOptions:")
            print("1. Select from database")
            print("2. Enter custom compound and Ksp")
            
            choice = input("\nEnter choice (1 or 2): ").strip()
            
            if choice == "1":
                # Use database compound
                compound_input = input(f"\nEnter compound name from list above: ").strip()
                
                if compound_input not in available_compounds:
                    print(f"Compound '{compound_input}' not found in database.")
                    print(f"Available: {', '.join(available_compounds)}")
                    return
                
                result = calculate_solubility_from_ksp(compound_input)
                
            elif choice == "2":
                # Custom compound and Ksp
                compound_input = input("\nEnter compound formula (e.g., AgCl, PbCl2): ").strip()
                
                try:
                    ksp_input = input("Enter Ksp value (scientific notation like 1.6e-10): ").strip()
                    ksp_value = float(ksp_input)
                except ValueError:
                    print("Invalid Ksp value format. Please use scientific notation like 1.6e-10")
                    return
                
                temp_input = input("Enter temperature in °C (press enter for 25°C): ").strip()
                temperature = 25.0 if not temp_input else float(temp_input)
                
                result = calculate_solubility_from_ksp(compound_input, ksp_value, temperature)
                
            else:
                print("Invalid choice.")
                return
            
            # Display results
            display_results_header()
            
            if "error" in result:
                print(f"Error: {result['error']}")
                if "available_compounds" in result:
                    print("Available compounds:", ", ".join(result["available_compounds"]))
                return
            
            # Display calculation steps
            display_steps(result["steps"])
            
            # Display final results in a formatted table
            print(f"\n{'='*60}")
            print(f"FINAL RESULTS FOR {result['compound']}")
            print(f"{'='*60}")
            print(f"Ksp value:                {result['ksp']:.2e}")
            print(f"Molar mass:               {result['molar_mass']:.2f} g/mol")
            print(f"Temperature:              {result['temperature']}°C")
            print(f"")
            print(f"Molar solubility:         {result['molar_solubility']:.2e} mol/L")
            print(f"Mass solubility:          {result['mass_solubility_g_L']:.2e} g/L")
            print(f"Mass solubility:          {result['mass_solubility_mg_L']:.2e} mg/L")
            print(f"Mass solubility:          {result['mass_solubility_g_100mL']:.2e} g/100mL")
            print(f"{'='*60}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_predefined_scenarios(self):
        """Handle predefined analysis scenarios."""
        print("\n===== PREDEFINED ANALYSIS SCENARIOS =====")
        
        scenarios = {
            "1": {
                "id": "W20_8",
                "description": "Qualitative Analysis: Determine if unknown solution contains Ag+, Ba2+, or Pb2+"
            },
            "2": {
                "id": "S21_11",
                "description": "Ksp Calculation: AgCl solubility from Ksp = 1.6×10⁻¹⁰"
            }
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
            
            # Handle different types of results
            if "steps" in result:
                # Display analysis steps
                display_steps(result["steps"])
                
                if "conclusion" in result:
                    # Qualitative analysis result
                    print("\nConclusion:")
                    print(result["conclusion"])
                    
                    if "identified_cations" in result and result["identified_cations"]:
                        print("\nIdentified cation(s):", ", ".join(result["identified_cations"]))
                    elif "identified_cations" in result:
                        print("\nNo cation could be identified with the given constraints.")
                else:
                    # Ksp calculation result
                    print(f"\n{'='*60}")
                    print(f"FINAL RESULTS FOR {result['compound']}")
                    print(f"{'='*60}")
                    print(f"Ksp value:                {result['ksp']:.2e}")
                    print(f"Molar mass:               {result['molar_mass']:.2f} g/mol")
                    print(f"Temperature:              {result['temperature']}°C")
                    print(f"")
                    print(f"Molar solubility:         {result['molar_solubility']:.2e} mol/L")
                    print(f"Mass solubility:          {result['mass_solubility_g_L']:.2e} g/L")
                    print(f"Mass solubility:          {result['mass_solubility_mg_L']:.2e} mg/L")
                    print(f"Mass solubility:          {result['mass_solubility_g_100mL']:.2e} g/100mL")
                    print(f"{'='*60}")
        else:
            print("Invalid choice.")
    
    def _handle_compound_database(self):
        """Display information about compounds in the database."""
        print("\n===== COMPOUND DATABASE =====")
        
        compounds = get_available_compounds()
        
        if not compounds:
            print("No compounds found in database.")
            return
        
        print(f"\nDatabase contains {len(compounds)} compounds with Ksp data:")
        print(f"{'Compound':<15} {'Ksp Value':<12} {'Molar Mass':<12} {'Type':<8}")
        print("-" * 50)
        
        # Import the KSP_DATA to display detailed information
        from chemistry_solver.insoluble_salts import KSP_DATA
        
        for compound in sorted(compounds):
            data = KSP_DATA[compound]
            print(f"{compound:<15} {data['ksp']:<12.2e} {data['molar_mass']:<12.2f} {data['type']:<8}")
        
        print(f"\nCompound Types:")
        print(f"  AB   - Binary compound (e.g., AgCl → Ag⁺ + Cl⁻)")
        print(f"  AB2  - 1:2 compound (e.g., PbCl₂ → Pb²⁺ + 2Cl⁻)")
        print(f"  A2B  - 2:1 compound (e.g., Ag₂CrO₄ → 2Ag⁺ + CrO₄²⁻)")
        print(f"  AB3  - 1:3 compound (e.g., Fe(OH)₃ → Fe³⁺ + 3OH⁻)")
        
        # Show available cations and reagents too
        print(f"\nAvailable cations for qualitative analysis:")
        cations = get_available_cations()
        print(", ".join(cations))
        
        print(f"\nAvailable reagents for qualitative analysis:")
        reagents = get_available_reagents()
        print(", ".join(reagents))
    
    def _handle_solubility_rules(self):
        """Display information about solubility rules."""
        print("\n===== SOLUBILITY RULES REFERENCE =====")
        
        solubility_rules = """
General Solubility Rules:
------------------------
1. Most salts containing Group 1 elements (Li⁺, Na⁺, K⁺, etc.) are soluble.
2. Most salts containing ammonium ion (NH₄⁺) are soluble.
3. Most nitrates (NO₃⁻), acetates (CH₃COO⁻), and perchlorates (ClO₄⁻) are soluble.
4. Most chlorides (Cl⁻), bromides (Br⁻), and iodides (I⁻) are soluble.
   Exceptions: AgCl, PbCl₂, Hg₂Cl₂, HgCl₂ are insoluble or slightly soluble.
5. Most sulfates (SO₄²⁻) are soluble.
   Exceptions: CaSO₄, BaSO₄, PbSO₄, SrSO₄, Ag₂SO₄ are insoluble or slightly soluble.
6. Most hydroxides (OH⁻) are insoluble.
   Exceptions: Group 1 hydroxides and Ba(OH)₂ are soluble.
7. Most sulfides (S²⁻), carbonates (CO₃²⁻), phosphates (PO₄³⁻),
   and chromates (CrO₄²⁻) are insoluble.
   Exceptions: Group 1 and NH₄⁺ salts of these anions are soluble.

Common Precipitation Reagents:
----------------------------
- AgNO₃: Test for halides (Cl⁻, Br⁻, I⁻)
- BaCl₂: Test for sulfates (SO₄²⁻)
- NaOH: Test for metal hydroxides
- H₂SO₄: Test for Ba²⁺, Pb²⁺
- HCl: Test for Ag⁺, Pb²⁺
- (NH₄)₂CO₃: Test for carbonates

Ksp Expression Examples:
-----------------------
For AB type (e.g., AgCl):
  AgCl(s) ⇌ Ag⁺(aq) + Cl⁻(aq)
  Ksp = [Ag⁺][Cl⁻]

For AB₂ type (e.g., PbCl₂):
  PbCl₂(s) ⇌ Pb²⁺(aq) + 2Cl⁻(aq)
  Ksp = [Pb²⁺][Cl⁻]²

For A₂B type (e.g., Ag₂CrO₄):
  Ag₂CrO₄(s) ⇌ 2Ag⁺(aq) + CrO₄²⁻(aq)
  Ksp = [Ag⁺]²[CrO₄²⁻]

For AB₃ type (e.g., Fe(OH)₃):
  Fe(OH)₃(s) ⇌ Fe³⁺(aq) + 3OH⁻(aq)
  Ksp = [Fe³⁺][OH⁻]³

Solubility Calculations:
-----------------------
From Ksp to molar solubility (s):
- AB type:  s = √Ksp
- AB₂ type: s = ∛(Ksp/4)
- A₂B type: s = ∛(Ksp/4)
- AB₃ type: s = ⁴√(Ksp/27)

Mass solubility = molar solubility × molar mass
        """
        
        print(solubility_rules)