"""
Terminal User Interface for Colligative Properties Calculator
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user
from chemistry_solver.colligative_properties import (
    ColligativePropertyCalculator,
    calculate_molecular_weight,
    calculate_molality_and_freezing_point,
    calculate_molality_and_boiling_point,
    compare_colligative_properties,
    solve_multiple_choice_problem
)

class ColligativePropertiesUI:
    """UI class for colligative properties calculations."""
    
    def __init__(self):
        self.title = "COLLIGATIVE PROPERTIES CALCULATOR"
    
    def run(self):
        """Run the colligative properties UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-7): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_molality_freezing_point()
            elif choice == "2":
                self._handle_molality_boiling_point()
            elif choice == "3":
                self._handle_molecular_weight_calculation()
            elif choice == "4":
                self._handle_solution_comparison()
            elif choice == "5":
                self._handle_osmotic_pressure()
            elif choice == "6":
                self._handle_vapor_pressure()
            elif choice == "7":
                self._handle_about_colligative_properties()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the colligative properties module menu."""
        menu = """
        [1] Calculate molality and freezing point
        [2] Calculate molality and boiling point
        [3] Calculate molecular weight from colligative data
        [4] Compare solutions by colligative properties
        [5] Osmotic pressure calculations
        [6] Vapor pressure calculations
        [7] About colligative properties
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_molality_freezing_point(self):
        """Handle molality and freezing point calculations."""
        print("\n===== MOLALITY AND FREEZING POINT CALCULATOR =====")
        
        try:
            # Get input data
            solute_mass = float(input("Enter mass of solute (g): "))
            solute_molar_mass = float(input("Enter molar mass of solute (g/mol): "))
            solvent_mass = float(input("Enter mass of solvent (g): "))
            solvent = input("Enter solvent name (default: water): ").strip() or "water"
            
            # Get van 't Hoff factor
            ionization_input = input("Enter van 't Hoff factor (i) [default: 1]: ").strip()
            ionization_factor = float(ionization_input) if ionization_input else 1.0
            
            # Check for multiple choice options
            mc_input = input("Do you have multiple choice options? (y/n): ").strip().lower()
            answer_choices = None
            
            if mc_input == 'y':
                print("Enter answer choices as (molality, freezing_point) pairs.")
                print("Example: 0.08,5.6 0.06,5.4 0.04,5.4")
                choices_input = input("Enter choices (space-separated): ")
                try:
                    choices = []
                    for choice_str in choices_input.split():
                        m, fp = choice_str.split(',')
                        choices.append((float(m), float(fp)))
                    answer_choices = choices
                except:
                    print("Invalid format for choices. Proceeding without multiple choice.")
            
            # Calculate
            result = calculate_molality_and_freezing_point(
                solute_mass, solute_molar_mass, solvent_mass, 
                solvent, ionization_factor, answer_choices
            )
            
            # Display results
            display_results_header()
            if result['success']:
                for step in result['steps']:
                    print(step)
                
                print(f"\nFinal Results:")
                print(f"Molality: {result['molality']:.4f} mol/kg")
                print(f"Freezing point depression: {result['freezing_point_depression']:.2f}°C")
                print(f"New freezing point: {result['new_freezing_point']:.2f}°C")
                
                if 'closest_answer' in result:
                    print(f"\nClosest multiple choice answer: {result['closest_answer']}")
            else:
                print(f"Error: {result.get('error', 'Unknown error')}")
                
        except ValueError as e:
            print(f"Error: Invalid input - {str(e)}")
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_molality_boiling_point(self):
        """Handle molality and boiling point calculations."""
        print("\n===== MOLALITY AND BOILING POINT CALCULATOR =====")
        
        try:
            # Get input data
            solute_mass = float(input("Enter mass of solute (g): "))
            solute_molar_mass = float(input("Enter molar mass of solute (g/mol): "))
            solvent_mass = float(input("Enter mass of solvent (g): "))
            solvent = input("Enter solvent name (default: water): ").strip() or "water"
            
            # Get van 't Hoff factor
            ionization_input = input("Enter van 't Hoff factor (i) [default: 1]: ").strip()
            ionization_factor = float(ionization_input) if ionization_input else 1.0
            
            # Check for multiple choice options
            mc_input = input("Do you have multiple choice options? (y/n): ").strip().lower()
            answer_choices = None
            
            if mc_input == 'y':
                print("Enter answer choices as (molality, boiling_point) pairs.")
                print("Example: 0.08,100.4 0.06,100.3 0.04,100.2")
                choices_input = input("Enter choices (space-separated): ")
                try:
                    choices = []
                    for choice_str in choices_input.split():
                        m, bp = choice_str.split(',')
                        choices.append((float(m), float(bp)))
                    answer_choices = choices
                except:
                    print("Invalid format for choices. Proceeding without multiple choice.")
            
            # Calculate
            result = calculate_molality_and_boiling_point(
                solute_mass, solute_molar_mass, solvent_mass, 
                solvent, ionization_factor, answer_choices
            )
            
            # Display results
            display_results_header()
            if result['success']:
                for step in result['steps']:
                    print(step)
                
                print(f"\nFinal Results:")
                print(f"Molality: {result['molality']:.4f} mol/kg")
                print(f"Boiling point elevation: {result['boiling_point_elevation']:.2f}°C")
                print(f"New boiling point: {result['new_boiling_point']:.2f}°C")
                
                if 'closest_answer' in result:
                    print(f"\nClosest multiple choice answer: {result['closest_answer']}")
            else:
                print(f"Error: {result.get('error', 'Unknown error')}")
                
        except ValueError as e:
            print(f"Error: Invalid input - {str(e)}")
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_molecular_weight_calculation(self):
        """Handle molecular weight calculations from colligative data."""
        print("\n===== MOLECULAR WEIGHT CALCULATOR =====")
        
        try:
            # Choose method
            print("\nAvailable methods:")
            print("1. Freezing point depression")
            print("2. Boiling point elevation")
            print("3. Osmotic pressure")
            print("4. Vapor pressure lowering")
            
            method_choice = input("Choose method (1-4): ").strip()
            method_map = {
                "1": "freezing_point",
                "2": "boiling_point", 
                "3": "osmotic_pressure",
                "4": "vapor_pressure"
            }
            
            if method_choice not in method_map:
                print("Invalid method choice.")
                return
            
            method = method_map[method_choice]
            solvent = input("Enter solvent name (default: water): ").strip() or "water"
            
            # Get common parameters
            solute_mass = float(input("Enter mass of solute (g): "))
            
            # Get van 't Hoff factor
            ionization_input = input("Enter van 't Hoff factor (i) [default: 1]: ").strip()
            ionization_factor = float(ionization_input) if ionization_input else 1.0
            
            # Method-specific parameters
            kwargs = {'ionization_factor': ionization_factor}
            
            if method in ["freezing_point", "boiling_point"]:
                delta_T = float(input(f"Enter temperature change (°C): "))
                solvent_mass = float(input("Enter mass of solvent (g): "))
                kwargs.update({
                    'delta_T': delta_T,
                    'solute_mass': solute_mass,
                    'solvent_mass': solvent_mass
                })
                
            elif method == "osmotic_pressure":
                osmotic_pressure = float(input("Enter osmotic pressure (atm): "))
                temperature_c = float(input("Enter temperature (°C): "))
                solution_volume = float(input("Enter solution volume (L): "))
                kwargs.update({
                    'delta_T': osmotic_pressure,  # Using delta_T parameter name
                    'solute_mass': solute_mass,
                    'temperature_c': temperature_c,
                    'solution_volume_L': solution_volume
                })
                
            elif method == "vapor_pressure":
                P_pure = float(input("Enter pure solvent vapor pressure: "))
                P_solution = float(input("Enter solution vapor pressure: "))
                solvent_formula = input("Enter solvent formula: ")
                solvent_mass = float(input("Enter mass of solvent (g): "))
                kwargs.update({
                    'P_pure': P_pure,
                    'P_solution': P_solution,
                    'solute_mass': solute_mass,
                    'solvent_formula': solvent_formula,
                    'solvent_mass': solvent_mass
                })
            
            # Check for multiple choice
            mc_input = input("Do you have multiple choice options? (y/n): ").strip().lower()
            if mc_input == 'y':
                choices_input = input("Enter MW choices (space-separated): ")
                try:
                    answer_choices = [float(x) for x in choices_input.split()]
                    kwargs['answer_choices'] = answer_choices
                except:
                    print("Invalid format for choices.")
            
            # Calculate
            result = calculate_molecular_weight(method, solvent, **kwargs)
            
            # Display results
            display_results_header()
            if result['success']:
                for step in result['steps']:
                    print(step)
                
                print(f"\nCalculated Molecular Weight: {result['molecular_weight']:.2f} g/mol")
                
                if 'closest_answer' in result:
                    print(f"Closest multiple choice answer: {result['closest_answer']} g/mol")
            else:
                print(f"Error: {result.get('error', 'Unknown error')}")
                
        except ValueError as e:
            print(f"Error: Invalid input - {str(e)}")
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_solution_comparison(self):
        """Handle comparison of solutions by colligative properties."""
        print("\n===== SOLUTION COMPARISON =====")
        
        try:
            # Choose property type
            print("\nProperty types:")
            print("1. Freezing point depression")
            print("2. Boiling point elevation")
            
            prop_choice = input("Choose property type (1-2): ").strip()
            property_map = {
                "1": "freezing_point_depression",
                "2": "boiling_point_elevation"
            }
            
            if prop_choice not in property_map:
                print("Invalid property choice.")
                return
            
            property_type = property_map[prop_choice]
            solvent = input("Enter solvent name (default: water): ").strip() or "water"
            
            # Get solutions data
            num_solutions = int(input("Enter number of solutions to compare: "))
            solutions = []
            
            for i in range(num_solutions):
                print(f"\nSolution {i+1}:")
                name = input(f"Enter name for solution {i+1} (optional): ").strip()
                molality = float(input("Enter molality (mol/kg): "))
                
                ionization_input = input("Enter van 't Hoff factor (i) [default: 1]: ").strip()
                ionization_factor = float(ionization_input) if ionization_input else 1.0
                
                solution = {
                    'molality': molality,
                    'ionization_factor': ionization_factor
                }
                
                if name:
                    solution['name'] = name
                
                solutions.append(solution)
            
            # Compare solutions
            result = compare_colligative_properties(solutions, property_type, solvent)
            
            # Display results
            display_results_header()
            if result['success']:
                for step in result['steps']:
                    print(step)
            else:
                print(f"Error: {result.get('error', 'Unknown error')}")
                
        except ValueError as e:
            print(f"Error: Invalid input - {str(e)}")
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_osmotic_pressure(self):
        """Handle osmotic pressure calculations."""
        print("\n===== OSMOTIC PRESSURE CALCULATOR =====")
        print("This feature calculates molecular weight from osmotic pressure data.")
        print("Use option 3 (Molecular Weight Calculator) and select osmotic pressure method.")
    
    def _handle_vapor_pressure(self):
        """Handle vapor pressure calculations."""
        print("\n===== VAPOR PRESSURE CALCULATOR =====")
        print("This feature calculates molecular weight from vapor pressure lowering data.")
        print("Use option 3 (Molecular Weight Calculator) and select vapor pressure method.")
    
    def _handle_about_colligative_properties(self):
        """Display information about colligative properties."""
        print("\n===== ABOUT COLLIGATIVE PROPERTIES =====")
        
        info = """
Colligative properties are properties of solutions that depend only on the 
number of solute particles present, not on the identity of the solute particles.

The four main colligative properties are:

1. FREEZING POINT DEPRESSION (ΔTf)
   - Formula: ΔTf = Kf × m × i
   - Kf = freezing point depression constant (depends on solvent)
   - m = molality of solution
   - i = van 't Hoff factor (number of particles per formula unit)

2. BOILING POINT ELEVATION (ΔTb)
   - Formula: ΔTb = Kb × m × i
   - Kb = boiling point elevation constant (depends on solvent)
   - m = molality of solution
   - i = van 't Hoff factor

3. OSMOTIC PRESSURE (π)
   - Formula: π = iMRT
   - i = van 't Hoff factor
   - M = molarity of solution
   - R = gas constant (0.08206 L·atm/(mol·K))
   - T = temperature in Kelvin

4. VAPOR PRESSURE LOWERING (ΔP)
   - Formula: ΔP/P° = X_solute
   - ΔP = vapor pressure lowering
   - P° = vapor pressure of pure solvent
   - X_solute = mole fraction of solute

Applications:
- Determining molecular weights of unknown compounds
- Comparing the relative effects of different solutes
- Calculating concentrations from property measurements
- Understanding solution behavior in various conditions

Note: The van 't Hoff factor (i) accounts for dissociation:
- i = 1 for non-electrolytes (glucose, sucrose)
- i = 2 for binary electrolytes (NaCl, KBr)
- i = 3 for ternary electrolytes (CaCl₂, Na₂SO₄)

Decision tree for determining Van't Hoff factor:

Is it an ionic compound (salt)? 
├─ YES → Count the ions (NaCl = 2, CaCl₂ = 3)
└─ NO → Is it a strong acid/base in water?
   ├─ YES → Count H⁺/OH⁻ produced  
   └─ NO → i = 1
        """
        
        print(info)