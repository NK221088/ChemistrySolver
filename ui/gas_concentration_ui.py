"""
Terminal User Interface for Gas Concentration Calculations
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user
from chemistry_solver.gas_concentration import ppm_to_mass, mass_to_ppm

class GasConcentrationUI:
    """UI class for gas concentration calculations."""
    
    def __init__(self):
        self.title = "GAS CONCENTRATION CALCULATOR"
    
    def run(self):
        """Run the gas concentration UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-2): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_ppm_to_mass()
            elif choice == "2":
                self._handle_mass_to_ppm()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the gas concentration module menu."""
        menu = """
        [1] Calculate mass from concentration (ppm to grams)
        [2] Calculate concentration from mass (grams to ppm)
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