"""
Terminal User Interface for Gas Concentration Calculations
Extended with gas density calculations.
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user
from chemistry_solver.gas_concentration import (
    ppm_to_mass, mass_to_ppm, calculate_gas_density, 
    identify_gas_by_density, solve_density_problem
)

class GasConcentrationUI:
    """UI class for gas concentration calculations."""
    
    def __init__(self):
        self.title = "GAS CONCENTRATION & DENSITY CALCULATOR"
    
    def run(self):
        """Run the gas concentration UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-5): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_ppm_to_mass()
            elif choice == "2":
                self._handle_mass_to_ppm()
            elif choice == "3":
                self._handle_gas_density()
            elif choice == "4":
                self._handle_identify_by_density()
            elif choice == "5":
                self._handle_solve_chemistry_problem()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the gas concentration module menu."""
        menu = """
        [1] Calculate mass from concentration (ppm to grams)
        [2] Calculate concentration from mass (grams to ppm)
        [3] Calculate gas density at given conditions
        [4] Identify gas by density (multiple candidates)
        [5] Solve chemistry density problem (like exam questions)
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_ppm_to_mass(self):
        """Handle calculation of mass from concentration in ppm."""
        print("\n===== PPM TO MASS CALCULATOR =====")
        
        try:
            formula = input("Enter chemical formula of the gas: ")
            ppm = float(input("Enter concentration in ppm: "))
            volume = float(input("Enter room volume in m³: "))
            
            temp_input = input("Enter temperature in °C (default 25): ").strip()
            temperature = float(temp_input) if temp_input else 25
            
            pressure_input = input("Enter pressure in atm (default 1): ").strip()
            pressure = float(pressure_input) if pressure_input else 1
            
            result = ppm_to_mass(formula, ppm, volume, temperature, pressure)
            
            display_results_header()
            if result['success']:
                print(f"Formula: {result['formula']}")
                print(f"Molar Mass: {result['molar_mass']:.4f} g/mol")
                print(f"Concentration: {result['ppm']:.2f} ppm")
                print(f"Room Volume: {result['volume_m3']:.2f} m³ ({result['volume_liters']:.2f} liters)")
                print(f"Temperature: {result['temperature_celsius']:.2f}°C ({result['temperature_kelvin']:.2f} K)")
                print(f"Pressure: {result['pressure_atm']:.2f} atm")
                print(f"\nCalculated Mass: {result['mass_grams']:.4f} g")
                print(f"\nDetails:")
                print(f"- Moles of air: {result['moles_air']:.4f} mol")
                print(f"- Moles of gas: {result['moles_gas']:.8f} mol")
            else:
                print(f"Error: {result['error']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_mass_to_ppm(self):
        """Handle calculation of concentration in ppm from mass."""
        print("\n===== MASS TO PPM CALCULATOR =====")
        
        try:
            formula = input("Enter chemical formula of the gas: ")
            mass = float(input("Enter mass in grams: "))
            volume = float(input("Enter room volume in m³: "))
            
            temp_input = input("Enter temperature in °C (default 25): ").strip()
            temperature = float(temp_input) if temp_input else 25
            
            pressure_input = input("Enter pressure in atm (default 1): ").strip()
            pressure = float(pressure_input) if pressure_input else 1
            
            result = mass_to_ppm(formula, mass, volume, temperature, pressure)
            
            display_results_header()
            if result['success']:
                print(f"Formula: {result['formula']}")
                print(f"Molar Mass: {result['molar_mass']:.4f} g/mol")
                print(f"Mass: {result['mass_grams']:.4f} g")
                print(f"Room Volume: {result['volume_m3']:.2f} m³ ({result['volume_liters']:.2f} liters)")
                print(f"Temperature: {result['temperature_celsius']:.2f}°C ({result['temperature_kelvin']:.2f} K)")
                print(f"Pressure: {result['pressure_atm']:.2f} atm")
                print(f"\nCalculated Concentration: {result['ppm']:.4f} ppm")
                print(f"\nDetails:")
                print(f"- Moles of air: {result['moles_air']:.4f} mol")
                print(f"- Moles of gas: {result['moles_gas']:.8f} mol")
            else:
                print(f"Error: {result['error']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_gas_density(self):
        """Handle calculation of gas density at given conditions."""
        print("\n===== GAS DENSITY CALCULATOR =====")
        
        try:
            formula = input("Enter chemical formula of the gas: ")
            
            temp_input = input("Enter temperature in °C (default 25): ").strip()
            temperature = float(temp_input) if temp_input else 25
            
            pressure_input = input("Enter pressure in atm (default 1): ").strip()
            pressure = float(pressure_input) if pressure_input else 1
            
            result = calculate_gas_density(formula, temperature, pressure)
            
            display_results_header()
            if result['success']:
                print(f"Formula: {result['formula']}")
                print(f"Molar Mass: {result['molar_mass']:.4f} g/mol")
                print(f"Temperature: {result['temperature_celsius']:.2f}°C ({result['temperature_kelvin']:.2f} K)")
                print(f"Pressure: {result['pressure_atm']:.2f} atm")
                print(f"\nCalculated Density: {result['density_g_per_L']:.4f} g/L")
                print(f"\nCalculation: density = (P × M) / (R × T)")
                print(f"            density = ({pressure} × {result['molar_mass']:.4f}) / (0.08206 × {result['temperature_kelvin']:.2f})")
                print(f"            density = {result['density_g_per_L']:.4f} g/L")
            else:
                print(f"Error: {result['error']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_identify_by_density(self):
        """Handle identification of gas by comparing densities of multiple candidates."""
        print("\n===== IDENTIFY GAS BY DENSITY =====")
        
        try:
            print("Enter candidate chemical formulas (one per line, press Enter twice when done):")
            formulas = []
            while True:
                formula = input().strip()
                if not formula:
                    break
                formulas.append(formula)
            
            if not formulas:
                print("No formulas entered.")
                return
            
            target_density = float(input("Enter target density in g/L: "))
            
            temp_input = input("Enter temperature in °C (default 25): ").strip()
            temperature = float(temp_input) if temp_input else 25
            
            pressure_input = input("Enter pressure in atm (default 1): ").strip()
            pressure = float(pressure_input) if pressure_input else 1
            
            tolerance_input = input("Enter tolerance in g/L (default 0.01): ").strip()
            tolerance = float(tolerance_input) if tolerance_input else 0.01
            
            result = identify_gas_by_density(formulas, target_density, temperature, pressure, tolerance)
            
            display_results_header()
            if result['success']:
                print(f"Target density: {result['target_density']:.4f} g/L")
                print(f"Conditions: {result['temperature_celsius']:.2f}°C, {result['pressure_atm']:.2f} atm")
                print(f"Tolerance: ±{result['tolerance']:.4f} g/L")
                
                print(f"\nCalculated densities for all candidates:")
                for entry in result['all_results']:
                    if 'calculated_density' in entry:
                        print(f"  {entry['formula']}: {entry['calculated_density']:.4f} g/L (difference: {entry['difference']:.4f})")
                    else:
                        print(f"  {entry['formula']}: Error - {entry['error']}")
                
                if result['matches']:
                    print(f"\nMatches within tolerance:")
                    for match in result['matches']:
                        print(f"  ✓ {match['formula']} - {match['calculated_density']:.4f} g/L")
                else:
                    print(f"\nNo exact matches within tolerance.")
                    if result['best_match']:
                        best = result['best_match']
                        print(f"Closest match: {best['formula']} - {best['calculated_density']:.4f} g/L")
            else:
                print(f"Error: {result['error']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_solve_chemistry_problem(self):
        """Handle solving chemistry problems like the exam question."""
        print("\n===== CHEMISTRY PROBLEM SOLVER =====")
        print("This mode is designed for problems like:")
        print("'Which compound has a density of X g/L at Y°C and Z atm?'")
        
        try:
            print("\nEnter candidate chemical formulas (one per line, press Enter twice when done):")
            formulas = []
            while True:
                formula = input().strip()
                if not formula:
                    break
                formulas.append(formula)
            
            if not formulas:
                print("No formulas entered.")
                return
            
            target_density = float(input("Enter target density in g/L: "))
            temperature = float(input("Enter temperature in °C: "))
            pressure = float(input("Enter pressure in atm: "))
            
            # This function will print detailed step-by-step solution
            result = solve_density_problem(formulas, target_density, temperature, pressure)
            
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def solve_example_problem(self):
        """Solve the specific problem from your chemistry question."""
        print("=== SOLVING YOUR CHEMISTRY QUESTION ===")
        
        # The specific problem: density of 2.61 g/L at 100°C and 1 atm
        # Candidates: SO2, SO3, SO, S6O, S2O2
        candidate_formulas = ["SO2", "SO3", "SO", "S6O", "S2O2"]
        target_density = 2.61
        temperature_celsius = 100
        pressure_atm = 1
        
        result = solve_density_problem(candidate_formulas, target_density, temperature_celsius, pressure_atm)
        
        return result