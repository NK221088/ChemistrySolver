"""
Terminal User Interface for Chemical Kinetics Problems
"""
from ui.terminal_ui import display_title, display_results_header, display_steps, wait_for_user
from chemistry_solver.kinetics import (
    solve_first_order_kinetics, 
    calculate_rate_constant_from_half_life,
    calculate_half_life_from_rate_constant, 
    calculate_concentration_after_time,
    calculate_fraction_remaining, 
    as_simplified_fraction,
    solve_arrhenius_problem,
    celsius_to_kelvin,
    kelvin_to_celsius,
    calculate_rate_constant_arrhenius,
    calculate_temperature_for_rate_constant,
    calculate_activation_energy
)
import math

class KineticsUI:
    """UI class for chemical kinetics problems."""
    
    def __init__(self):
        self.title = "CHEMICAL KINETICS PROBLEM SOLVER"
    
    def run(self):
        """Run the chemical kinetics UI."""
        display_title(self.title)
        
        while True:
            self._display_main_menu()
            choice = input("\nEnter choice (0-2): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._run_first_order_kinetics()
            elif choice == "2":
                self._run_arrhenius_equation()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_main_menu(self):
        """Display the main chemical kinetics menu."""
        menu = """
        ===== CHEMICAL KINETICS MODULES =====
        [1] First-Order Kinetics Problems
        [2] Arrhenius Equation Problems
        [0] Return to main menu
        """
        print(menu)
    
    def _run_first_order_kinetics(self):
        """Run the first-order kinetics submenu."""
        while True:
            print("\n" + "="*50)
            print("FIRST-ORDER KINETICS PROBLEMS")
            print("="*50)
            self._display_first_order_menu()
            choice = input("\nEnter choice (0-6): ").strip()
            
            if choice == "0":
                return
            elif choice == "1":
                self._handle_rate_constant_calculation()
            elif choice == "2":
                self._handle_half_life_calculation()
            elif choice == "3":
                self._handle_concentration_calculation()
            elif choice == "4":
                self._handle_time_calculation()
            elif choice == "5":
                self._handle_fraction_calculation()
            elif choice == "6":
                self._handle_complete_analysis()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_first_order_menu(self):
        """Display the first-order kinetics submenu."""
        menu = """
        [1] Calculate rate constant from half-life
        [2] Calculate half-life from rate constant
        [3] Calculate concentration after time
        [4] Calculate time to reach concentration
        [5] Calculate fraction remaining after time
        [6] Complete first-order kinetics analysis
        [0] Back to kinetics menu
        """
        print(menu)
    
    def _run_arrhenius_equation(self):
        """Run the Arrhenius equation submenu."""
        while True:
            print("\n" + "="*50)
            print("ARRHENIUS EQUATION PROBLEMS")
            print("="*50)
            self._display_arrhenius_menu()
            choice = input("\nEnter choice (0-5): ").strip()
            
            if choice == "0":
                return
            elif choice == "1":
                self._handle_temperature_conversion()
            elif choice == "2":
                self._handle_rate_constant_arrhenius()
            elif choice == "3":
                self._handle_activation_energy_calculation()
            elif choice == "4":
                self._handle_temperature_for_rate_calculation()
            elif choice == "5":
                self._handle_complete_arrhenius_analysis()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_arrhenius_menu(self):
        """Display the Arrhenius equation submenu."""
        menu = """
        [1] Temperature conversion (°C ↔ K)
        [2] Calculate rate constant using Arrhenius equation
        [3] Calculate activation energy from two rate constants
        [4] Calculate temperature for desired rate constant
        [5] Complete Arrhenius equation analysis
        [0] Back to kinetics menu
        """
        print(menu)
    
    def _handle_rate_constant_calculation(self):
        """Handle calculation of rate constant from half-life."""
        print("\n===== RATE CONSTANT CALCULATOR =====")
        
        try:
            half_life = float(input("Enter half-life value: "))
            unit = input("Enter the unit of time (e.g., s, min, h): ")
            
            rate_constant = calculate_rate_constant_from_half_life(half_life)
            
            display_results_header()
            print(f"Half-life: {half_life} {unit}")
            print(f"Rate constant (k): {rate_constant:.6f} {unit}⁻¹")
            print("\nFormula used: k = ln(2) / t₁/₂")
            print(f"k = 0.693147 / {half_life} = {rate_constant:.6f} {unit}⁻¹")
        except ValueError:
            print("Error: Please enter a valid number for half-life.")
    
    def _handle_half_life_calculation(self):
        """Handle calculation of half-life from rate constant."""
        print("\n===== HALF-LIFE CALCULATOR =====")
        
        try:
            rate_constant = float(input("Enter rate constant (k) value: "))
            unit = input("Enter the unit of time (e.g., s, min, h): ")
            
            half_life = calculate_half_life_from_rate_constant(rate_constant)
            
            display_results_header()
            print(f"Rate constant (k): {rate_constant} {unit}⁻¹")
            print(f"Half-life: {half_life:.6f} {unit}")
            print("\nFormula used: t₁/₂ = ln(2) / k")
            print(f"t₁/₂ = 0.693147 / {rate_constant} = {half_life:.6f} {unit}")
        except ValueError:
            print("Error: Please enter a valid number for rate constant.")
    
    def _handle_concentration_calculation(self):
        """Handle calculation of concentration after time."""
        print("\n===== CONCENTRATION AFTER TIME CALCULATOR =====")
        
        try:
            initial_conc = float(input("Enter initial concentration [A]₀: "))
            rate_constant = float(input("Enter rate constant (k): "))
            time = float(input("Enter time elapsed: "))
            
            conc_unit = input("Enter concentration unit (e.g., M, mol/L): ")
            time_unit = input("Enter time unit (e.g., s, min, h): ")
            
            final_conc = calculate_concentration_after_time(initial_conc, rate_constant, time)
            fraction = calculate_fraction_remaining(rate_constant, time)
            
            display_results_header()
            print(f"Initial concentration [A]₀: {initial_conc} {conc_unit}")
            print(f"Rate constant (k): {rate_constant} {time_unit}⁻¹")
            print(f"Time elapsed: {time} {time_unit}")
            print(f"Final concentration [A]t: {final_conc:.6f} {conc_unit}")
            print(f"Fraction remaining: {fraction:.6f}")
            
            # Try to represent as simplified fraction if possible
            fraction_as_tuple = as_simplified_fraction(fraction)
            if fraction_as_tuple:
                num, den = fraction_as_tuple
                if den <= 100:  # Only show if denominator is reasonable
                    print(f"This fraction can be expressed as {num}/{den}")
            
            print("\nFormula used: [A]t = [A]₀ × e⁻ᵏᵗ")
            print(f"[A]t = {initial_conc} × e^(-{rate_constant} × {time})")
            print(f"[A]t = {initial_conc} × {fraction:.6f}")
            print(f"[A]t = {final_conc:.6f} {conc_unit}")
        except ValueError:
            print("Error: Please enter valid numbers for all values.")
    
    def _handle_time_calculation(self):
        """Handle calculation of time to reach a concentration."""
        print("\n===== TIME TO REACH CONCENTRATION CALCULATOR =====")
        
        try:
            initial_conc = float(input("Enter initial concentration [A]₀: "))
            final_conc = float(input("Enter final concentration [A]t: "))
            rate_constant = float(input("Enter rate constant (k): "))
            
            conc_unit = input("Enter concentration unit (e.g., M, mol/L): ")
            time_unit = input("Enter time unit (e.g., s, min, h): ")
            
            # Calculate time using ln([A]₀/[A]t) / k
            if final_conc <= 0 or final_conc >= initial_conc:
                raise ValueError("Final concentration must be positive and less than initial concentration")
            
            time = math.log(initial_conc/final_conc) / rate_constant
            
            display_results_header()
            print(f"Initial concentration [A]₀: {initial_conc} {conc_unit}")
            print(f"Final concentration [A]t: {final_conc} {conc_unit}")
            print(f"Rate constant (k): {rate_constant} {time_unit}⁻¹")
            print(f"Time required: {time:.6f} {time_unit}")
            
            print("\nFormula used: t = ln([A]₀/[A]t) / k")
            print(f"t = ln({initial_conc}/{final_conc}) / {rate_constant}")
            print(f"t = ln({initial_conc/final_conc:.6f}) / {rate_constant}")
            print(f"t = {time:.6f} {time_unit}")
        except ValueError as e:
            print(f"Error: {e}")
    
    def _handle_fraction_calculation(self):
        """Handle calculation of fraction remaining after time."""
        print("\n===== FRACTION REMAINING CALCULATOR =====")
        
        try:
            rate_constant = float(input("Enter rate constant (k): "))
            time = float(input("Enter time elapsed: "))
            time_unit = input("Enter time unit (e.g., s, min, h): ")
            
            fraction = calculate_fraction_remaining(rate_constant, time)
            
            display_results_header()
            print(f"Rate constant (k): {rate_constant} {time_unit}⁻¹")
            print(f"Time elapsed: {time} {time_unit}")
            print(f"Fraction remaining: {fraction:.6f}")
            
            # Try to represent as simplified fraction
            fraction_as_tuple = as_simplified_fraction(fraction)
            if fraction_as_tuple:
                num, den = fraction_as_tuple
                if den <= 100:
                    print(f"This can be expressed as {num}/{den}")
            
            print(f"Percentage remaining: {fraction * 100:.2f}%")
            print(f"Percentage decomposed: {(1 - fraction) * 100:.2f}%")
            
            print("\nFormula used: Fraction = e⁻ᵏᵗ")
            print(f"Fraction = e^(-{rate_constant} × {time}) = {fraction:.6f}")
            
        except ValueError:
            print("Error: Please enter valid numbers for all values.")
    
    def _handle_complete_analysis(self):
        """Handle complete first-order kinetics analysis."""
        print("\n=== Complete First-Order Kinetics Analysis ===")
        print("This performs a complete analysis given any two parameters.")
        
        print("\nWhat information do you have? Enter values for any two of the following:")
        print("(Leave blank if unknown)")
        
        half_life_input = input("Half-life (leave blank if unknown): ")
        half_life = float(half_life_input) if half_life_input else None
        
        rate_constant_input = input("Rate constant k (leave blank if unknown): ")
        rate_constant = float(rate_constant_input) if rate_constant_input else None
        
        time_input = input("Time elapsed (leave blank if unknown): ")
        time = float(time_input) if time_input else None
        
        fraction_input = input("Fraction remaining (e.g., 0.5 for half) (leave blank if unknown): ")
        fraction_remaining = float(fraction_input) if fraction_input else None
        
        initial_conc_input = input("Initial concentration [A]₀ (leave blank if unknown): ")
        initial_concentration = float(initial_conc_input) if initial_conc_input else None
        
        final_conc_input = input("Final concentration [A]t (leave blank if unknown): ")
        final_concentration = float(final_conc_input) if final_conc_input else None
        
        unit = input("Enter the unit of time (e.g., s, min, h): ")
        conc_unit = input("Enter the concentration unit (if applicable, e.g., M, mol/L): ") if initial_concentration or final_concentration else None
        
        try:
            result = solve_first_order_kinetics(
                half_life=half_life,
                rate_constant=rate_constant,
                initial_concentration=initial_concentration,
                final_concentration=final_concentration,
                time=time,
                fraction_remaining=fraction_remaining
            )
            
            display_results_header()
            display_steps(result["steps"])
            
            print(f"\nFinal Results:")
            if result["half_life"] is not None:
                print(f"Half-life: {result['half_life']:.6f} {unit}")
            if result["rate_constant"] is not None:
                print(f"Rate constant (k): {result['rate_constant']:.6f} {unit}⁻¹")
            if result["time"] is not None:
                print(f"Time elapsed: {result['time']:.6f} {unit}")
            if result["fraction_remaining"] is not None:
                print(f"Fraction remaining: {result['fraction_remaining']:.6f}")
                print(f"Percentage remaining: {result['fraction_remaining'] * 100:.2f}%")
                if result["fraction_as_tuple"]:
                    num, den = result["fraction_as_tuple"]
                    if den <= 100:  # Only show if denominator is reasonable
                        print(f"This fraction can be expressed as {num}/{den}")
            if result["initial_concentration"] is not None:
                print(f"Initial concentration [A]₀: {result['initial_concentration']:.6f} {conc_unit}")
            if result["final_concentration"] is not None:
                print(f"Final concentration [A]t: {result['final_concentration']:.6f} {conc_unit}")
                
        except Exception as e:
            print(f"\nError: {e}")
            print("Please provide sufficient information to solve the problem.")
    
    def _handle_temperature_conversion(self):
        """Handle temperature conversion between Celsius and Kelvin."""
        print("\n===== TEMPERATURE CONVERSION =====")
        
        print("Convert from:")
        print("[1] Celsius to Kelvin")
        print("[2] Kelvin to Celsius")
        
        choice = input("Enter choice (1 or 2): ").strip()
        
        try:
            if choice == "1":
                celsius = float(input("Enter temperature in Celsius: "))
                kelvin = celsius_to_kelvin(celsius)
                display_results_header()
                print(f"Temperature: {celsius}°C = {kelvin:.2f} K")
            elif choice == "2":
                kelvin = float(input("Enter temperature in Kelvin: "))
                celsius = kelvin_to_celsius(kelvin)
                display_results_header()
                print(f"Temperature: {kelvin} K = {celsius:.2f}°C")
            else:
                print("Invalid choice.")
        except ValueError:
            print("Error: Please enter a valid number for temperature.")
    
    def _handle_rate_constant_arrhenius(self):
        """Handle calculation of rate constant using Arrhenius equation."""
        print("\n===== RATE CONSTANT FROM ARRHENIUS EQUATION =====")
        
        try:
            activation_energy = float(input("Enter activation energy (J/mol): "))
            
            print("Enter temperature:")
            temp_choice = input("In Celsius (C) or Kelvin (K)? Enter C or K: ").strip().upper()
            
            if temp_choice == "C":
                temp_celsius = float(input("Enter temperature in Celsius: "))
                temp_kelvin = celsius_to_kelvin(temp_celsius)
            elif temp_choice == "K":
                temp_kelvin = float(input("Enter temperature in Kelvin: "))
                temp_celsius = kelvin_to_celsius(temp_kelvin)
            else:
                print("Invalid choice. Please enter C or K.")
                return
            
            pre_exponential = float(input("Enter pre-exponential factor (A): "))
            
            rate_constant = calculate_rate_constant_arrhenius(activation_energy, temp_kelvin, pre_exponential)
            
            display_results_header()
            print(f"Activation energy (Ea): {activation_energy:,.0f} J/mol")
            print(f"Temperature: {temp_celsius:.1f}°C = {temp_kelvin:.2f} K")
            print(f"Pre-exponential factor (A): {pre_exponential}")
            print(f"Rate constant (k): {rate_constant:.6e}")
            
            print("\nFormula used: k = A × e^(-Ea/RT)")
            print(f"k = {pre_exponential} × e^(-{activation_energy}/(8.314 × {temp_kelvin}))")
            print(f"k = {pre_exponential} × e^({-activation_energy/(8.314 * temp_kelvin):.6f})")
            print(f"k = {rate_constant:.6e}")
            
        except ValueError:
            print("Error: Please enter valid numbers for all values.")
    
    def _handle_activation_energy_calculation(self):
        """Handle calculation of activation energy from two rate constants."""
        print("\n===== ACTIVATION ENERGY CALCULATOR =====")
        
        try:
            print("Enter data for the first temperature:")
            k1 = float(input("Rate constant k1: "))
            
            temp1_choice = input("Temperature in Celsius (C) or Kelvin (K)? Enter C or K: ").strip().upper()
            if temp1_choice == "C":
                T1_celsius = float(input("Temperature T1 in Celsius: "))
                T1_kelvin = celsius_to_kelvin(T1_celsius)
            elif temp1_choice == "K":
                T1_kelvin = float(input("Temperature T1 in Kelvin: "))
                T1_celsius = kelvin_to_celsius(T1_kelvin)
            else:
                print("Invalid choice.")
                return
            
            print("\nEnter data for the second temperature:")
            k2 = float(input("Rate constant k2: "))
            
            temp2_choice = input("Temperature in Celsius (C) or Kelvin (K)? Enter C or K: ").strip().upper()
            if temp2_choice == "C":
                T2_celsius = float(input("Temperature T2 in Celsius: "))
                T2_kelvin = celsius_to_kelvin(T2_celsius)
            elif temp2_choice == "K":
                T2_kelvin = float(input("Temperature T2 in Kelvin: "))
                T2_celsius = kelvin_to_celsius(T2_kelvin)
            else:
                print("Invalid choice.")
                return
            
            activation_energy = calculate_activation_energy(k1, T1_kelvin, k2, T2_kelvin)
            
            display_results_header()
            print(f"Rate constant k1: {k1} at T1 = {T1_celsius:.1f}°C ({T1_kelvin:.2f} K)")
            print(f"Rate constant k2: {k2} at T2 = {T2_celsius:.1f}°C ({T2_kelvin:.2f} K)")
            print(f"Activation energy (Ea): {activation_energy:,.0f} J/mol")
            print(f"Activation energy (Ea): {activation_energy/1000:.1f} kJ/mol")
            
            print("\nFormula used: Ea = -R × ln(k2/k1) / (1/T2 - 1/T1)")
            ln_ratio = math.log(k2/k1)
            temp_term = (1/T2_kelvin) - (1/T1_kelvin)
            print(f"ln(k2/k1) = ln({k2}/{k1}) = {ln_ratio:.6f}")
            print(f"1/T2 - 1/T1 = 1/{T2_kelvin:.2f} - 1/{T1_kelvin:.2f} = {temp_term:.8f}")
            print(f"Ea = -8.314 × {ln_ratio:.6f} / {temp_term:.8f}")
            print(f"Ea = {activation_energy:.0f} J/mol")
            
        except ValueError:
            print("Error: Please enter valid numbers for all values.")
    
    def _handle_temperature_for_rate_calculation(self):
        """Handle calculation of temperature for desired rate constant."""
        print("\n===== TEMPERATURE FOR DESIRED RATE CONSTANT =====")
        
        try:
            # Get activation energy
            activation_energy = float(input("Enter activation energy (J/mol): "))
            if activation_energy <= 0:
                print("Error: Activation energy must be positive.")
                return
            
            print("\nEnter initial conditions:")
            
            # Get initial rate constant
            k1 = float(input("Initial rate constant k1: "))
            if k1 <= 0:
                print("Error: Rate constant must be positive.")
                return
            
            # Get temperature unit choice
            temp1_choice = input("Initial temperature in Celsius (C) or Kelvin (K)? Enter C or K: ").strip().upper()
             
            if temp1_choice == "C":
                T1_celsius = float(input("Initial temperature T1 in Celsius: "))
                if T1_celsius < -273.15:
                    print("Error: Temperature cannot be below absolute zero (-273.15°C).")
                    return
                T1_kelvin = celsius_to_kelvin(T1_celsius)
            elif temp1_choice == "K":
                T1_kelvin = float(input("Initial temperature T1 in Kelvin: "))
                if T1_kelvin <= 0:
                    print("Error: Temperature in Kelvin must be positive.")
                    return
                T1_celsius = kelvin_to_celsius(T1_kelvin)
            else:
                print("Error: Please enter 'C' for Celsius or 'K' for Kelvin.")
                return
            
            # Get desired rate constant
            k2 = float(input("\nEnter desired rate constant k2: "))
            if k2 <= 0:
                print("Error: Rate constant must be positive.")
                return
            
            # Calculate the required temperature
            T2_kelvin = calculate_temperature_for_rate_constant(activation_energy, k1, T1_kelvin, k2)
            T2_celsius = kelvin_to_celsius(T2_kelvin)
            
            # Display results
            display_results_header()
            print(f"Activation energy (Ea): {activation_energy:,.0f} J/mol")
            print(f"Initial conditions: k1 = {k1} at T1 = {T1_celsius:.1f}°C ({T1_kelvin:.2f} K)")
            print(f"Desired rate constant: k2 = {k2}")
            print(f"Required temperature: T2 = {T2_celsius:.1f}°C ({T2_kelvin:.2f} K)")
            
            temp_change = T2_celsius - T1_celsius
            rate_change = k2 / k1
            print(f"\nTemperature change: {temp_change:+.1f}°C")
            print(f"Rate constant change: {rate_change:.2f}× {'increase' if rate_change > 1 else 'decrease'}")
            
            print("\nFormula used: 1/T2 = 1/T1 - (R/Ea) × ln(k2/k1)")
            ln_ratio = math.log(k2/k1)
            print(f"ln(k2/k1) = ln({k2}/{k1}) = {ln_ratio:.6f}")
            print(f"1/T2 = 1/{T1_kelvin:.2f} - (8.314/{activation_energy}) × {ln_ratio:.6f}")
            print(f"1/T2 = {1/T1_kelvin:.8f} - {8.314/activation_energy * ln_ratio:.8f}")
            print(f"1/T2 = {1/T2_kelvin:.8f}")
            print(f"T2 = {T2_kelvin:.2f} K = {T2_celsius:.1f}°C")
            
        except ValueError:
            print("Error: Please enter valid numbers for all values.")
        except Exception as e:
            print(f"Error: An unexpected error occurred: {str(e)}")
        
        input("\nPress Enter to continue...")
    
    def _handle_complete_arrhenius_analysis(self):
        """Handle complete Arrhenius equation analysis."""
        print("\n=== Complete Arrhenius Equation Analysis ===")
        print("This performs a complete analysis given sufficient parameters.")
        
        print("\nWhat information do you have? (Leave blank if unknown)")
        
        ea_input = input("Activation energy in J/mol (leave blank if unknown): ")
        activation_energy = float(ea_input) if ea_input else None
        
        k1_input = input("First rate constant k1 (leave blank if unknown): ")
        k1 = float(k1_input) if k1_input else None
        
        print("First temperature T1:")
        t1c_input = input("  In Celsius (leave blank if unknown): ")
        T1_celsius = float(t1c_input) if t1c_input else None
        
        t1k_input = input("  In Kelvin (leave blank if unknown): ")
        T1_kelvin = float(t1k_input) if t1k_input else None
        
        k2_input = input("Second rate constant k2 (leave blank if unknown): ")
        k2 = float(k2_input) if k2_input else None
        
        print("Second temperature T2:")
        t2c_input = input("  In Celsius (leave blank if unknown): ")
        T2_celsius = float(t2c_input) if t2c_input else None
        
        t2k_input = input("  In Kelvin (leave blank if unknown): ")
        T2_kelvin = float(t2k_input) if t2k_input else None
        
        pef_input = input("Pre-exponential factor A (leave blank if unknown): ")
        pre_exponential_factor = float(pef_input) if pef_input else None
        
        try:
            result = solve_arrhenius_problem(
                activation_energy=activation_energy,
                k1=k1,
                T1_celsius=T1_celsius,
                T1_kelvin=T1_kelvin,
                k2=k2,
                T2_celsius=T2_celsius,
                T2_kelvin=T2_kelvin,
                pre_exponential_factor=pre_exponential_factor
            )
            
            display_results_header()
            display_steps(result["steps"])
            
            print(f"\nFinal Results:")
            if result["activation_energy"] is not None:
                print(f"Activation energy (Ea): {result['activation_energy']:,.0f} J/mol ({result['activation_energy']/1000:.1f} kJ/mol)")
            if result["k1"] is not None:
                print(f"Rate constant k1: {result['k1']}")
            if result["T1_celsius"] is not None and result["T1_kelvin"] is not None:
                print(f"Temperature T1: {result['T1_celsius']:.1f}°C ({result['T1_kelvin']:.2f} K)")
            if result["k2"] is not None:
                print(f"Rate constant k2: {result['k2']}")
            if result["T2_celsius"] is not None and result["T2_kelvin"] is not None:
                print(f"Temperature T2: {result['T2_celsius']:.1f}°C ({result['T2_kelvin']:.2f} K)")
            if result["pre_exponential_factor"] is not None:
                print(f"Pre-exponential factor (A): {result['pre_exponential_factor']}")
            
            # Additional insights
            if result["T1_celsius"] is not None and result["T2_celsius"] is not None:
                temp_change = result["T2_celsius"] - result["T1_celsius"]
                print(f"\nTemperature change: {temp_change:+.1f}°C")
            
            if result["k1"] is not None and result["k2"] is not None:
                rate_change = result["k2"] / result["k1"]
                print(f"Rate constant change: {rate_change:.2f}× {'increase' if rate_change > 1 else 'decrease'}")
                
        except Exception as e:
            print(f"\nError: {e}")
            print("Please provide sufficient information to solve the problem.")
            print("Common combinations that work:")
            print("- Ea, k1, T1, k2 → calculates T2")
            print("- k1, T1, k2, T2 → calculates Ea")
            print("- Ea, k1, T1, T2 → calculates k2")
            print("- Ea, T, A → calculates k using basic Arrhenius equation")