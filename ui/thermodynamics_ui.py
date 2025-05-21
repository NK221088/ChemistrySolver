"""
Terminal User Interface for Thermodynamics Calculations
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user
from chemistry_solver.thermodynamics import (calculate_heat, calculate_temperature_change, calculate_final_temperature,
                          calculate_molar_heat, solve_thermal_equilibrium, solve_thermal_equilibrium_with_molar_heat,
                          handle_heat_transfer_problem, handle_heat_transfer_with_molar_heat, solve_mixture_problem,
                          calculate_boiling_point_with_pressure, calculate_pressure_with_temperature, 
                          calculate_heat_of_vaporization)

try:
    from chemistry_solver.molar_mass import calculate_molar_mass
except ImportError:
    # Define a fallback function if the module is not available
    def calculate_molar_mass(formula):
        return {'success': False, 'error': 'molmass module not installed'}

class ThermodynamicsUI:
    """UI class for thermodynamics calculations."""
    
    def __init__(self):
        self.title = "THERMODYNAMICS CALCULATOR"
    
    def run(self):
        """Run the thermodynamics UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-7): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_heat_calculation()
            elif choice == "2":
                self._handle_temperature_change_calculation()
            elif choice == "3":
                self._handle_final_temperature_calculation()
            elif choice == "4":
                self._handle_thermal_equilibrium()
            elif choice == "5":
                self._handle_thermal_equilibrium_molar()
            elif choice == "6":
                self._handle_mixture_thermal_equilibrium()
            elif choice == "7":
                self._handle_clausius_clapeyron_calculations()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the thermodynamics module menu."""
        menu = """
        [1] Calculate heat transfer
        [2] Calculate temperature change
        [3] Calculate final temperature
        [4] Solve thermal equilibrium problem
        [5] Solve thermal equilibrium with molar heat capacities
        [6] Solve multi-substance thermal equilibrium problem
        [7] Clausius-Clapeyron equation calculations
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_heat_calculation(self):
        """
        Handler function for calculating heat using q = m × c × ΔT.
        """
        print("\n=== Heat Transfer Calculation ===")
        mass = float(input("Enter mass (g): "))
        specific_heat = float(input("Enter specific heat capacity (J/(g·K)): "))
        initial_temp = float(input("Enter initial temperature (°C): "))
        final_temp = float(input("Enter final temperature (°C): "))
        
        delta_t = final_temp - initial_temp
        heat = calculate_heat(mass, specific_heat, delta_t)
        
        display_results_header()
        print(f"Temperature change (ΔT): {delta_t:.2f} °C")
        print(f"Heat transferred (q): {heat:.2f} J")
        print(f"Heat transferred: {heat/1000:.2f} kJ")

    def _handle_temperature_change_calculation(self):
        """
        Handler function for calculating temperature change using ΔT = q / (m × c).
        """
        print("\n=== Temperature Change Calculation ===")
        heat = float(input("Enter heat energy (J): "))
        mass = float(input("Enter mass (g): "))
        specific_heat = float(input("Enter specific heat capacity (J/(g·K)): "))
        
        delta_t = calculate_temperature_change(heat, mass, specific_heat)
        
        display_results_header()
        print(f"Temperature change (ΔT): {delta_t:.2f} °C")

    def _handle_final_temperature_calculation(self):
        """
        Handler function for calculating final temperature.
        """
        print("\n=== Final Temperature Calculation ===")
        initial_temp = float(input("Enter initial temperature (°C): "))
        delta_t = float(input("Enter temperature change (°C): "))
        
        final_temp = calculate_final_temperature(initial_temp, delta_t)
        
        display_results_header()
        print(f"Final temperature: {final_temp:.2f} °C")

    def _handle_thermal_equilibrium(self):
        """
        Handler function for solving thermal equilibrium problems.
        """
        print("\n=== Thermal Equilibrium Problem ===")
        print("This will calculate the final temperature when two substances reach thermal equilibrium.")
        
        # Substance 1
        print("\nSubstance 1:")
        mass1 = float(input("Enter mass (g): "))
        specific_heat1 = float(input("Enter specific heat capacity (J/(g·K)): "))
        initial_temp1 = float(input("Enter initial temperature (°C): "))
        
        # Substance 2
        print("\nSubstance 2:")
        mass2 = float(input("Enter mass (g): "))
        specific_heat2 = float(input("Enter specific heat capacity (J/(g·K)): "))
        initial_temp2 = float(input("Enter initial temperature (°C): "))
        
        # Solve problem
        result = handle_heat_transfer_problem(
            mass1, specific_heat1, initial_temp1,
            mass2, specific_heat2, initial_temp2
        )
        
        display_results_header()
        for step in result["steps"]:
            print(step)
        
        print(f"\nFinal equilibrium temperature: {result['final_temp']:.2f} °C")
        print(f"Heat transferred by substance 1: {result['heat_transferred_1']:.2f} J")
        print(f"Heat transferred by substance 2: {result['heat_transferred_2']:.2f} J")

    def _handle_thermal_equilibrium_molar(self):
        """
        Handler function for solving thermal equilibrium problems using molar heat capacities.
        """
        print("\n=== Thermal Equilibrium Problem (Molar Heat Capacities) ===")
        print("This will calculate the final temperature when two substances reach thermal equilibrium.")
        
        # Substance 1
        print("\nSubstance 1:")
        name1 = input("Enter substance name or chemical formula: ")
        mass1 = float(input("Enter mass (g): "))
        
        # Try to calculate molar mass if a chemical formula is given
        try:
            molar_mass_result = calculate_molar_mass(name1)
            if molar_mass_result['success']:
                molar_mass1 = molar_mass_result['molar_mass']
                print(f"Calculated molar mass: {molar_mass1:.4f} g/mol")
            else:
                print("Could not automatically calculate molar mass.")
                molar_mass1 = float(input("Enter molar mass (g/mol): "))
        except Exception:
            print("Could not automatically calculate molar mass.")
            molar_mass1 = float(input("Enter molar mass (g/mol): "))
        
        molar_heat_capacity1 = float(input("Enter molar heat capacity (J/(mol·K)): "))
        initial_temp1 = float(input("Enter initial temperature (°C): "))
        
        # Substance 2
        print("\nSubstance 2:")
        name2 = input("Enter substance name or chemical formula: ")
        mass2 = float(input("Enter mass (g): "))
        
        # Try to calculate molar mass if a chemical formula is given
        try:
            molar_mass_result = calculate_molar_mass(name2)
            if molar_mass_result['success']:
                molar_mass2 = molar_mass_result['molar_mass']
                print(f"Calculated molar mass: {molar_mass2:.4f} g/mol")
            else:
                print("Could not automatically calculate molar mass.")
                molar_mass2 = float(input("Enter molar mass (g/mol): "))
        except Exception:
            print("Could not automatically calculate molar mass.")
            molar_mass2 = float(input("Enter molar mass (g/mol): "))
        
        molar_heat_capacity2 = float(input("Enter molar heat capacity (J/(mol·K)): "))
        initial_temp2 = float(input("Enter initial temperature (°C): "))
        
        # Solve problem
        result = handle_heat_transfer_with_molar_heat(
            mass1, molar_mass1, molar_heat_capacity1, initial_temp1,
            mass2, molar_mass2, molar_heat_capacity2, initial_temp2
        )
        
        display_results_header()
        for step in result["steps"]:
            print(step)
        
        print(f"\nFinal equilibrium temperature: {result['final_temp']:.2f} °C")
        print(f"Heat transferred by substance 1: {result['heat_transferred_1']:.2f} J")
        print(f"Heat transferred by substance 2: {result['heat_transferred_2']:.2f} J")

    def _handle_mixture_thermal_equilibrium(self):
        """
        Handler function for solving thermal equilibrium problems with multiple substances.
        """
        print("\n=== Multiple Substance Thermal Equilibrium Problem ===")
        print("This will calculate the final temperature when multiple substances reach thermal equilibrium.")
        
        substance_count = int(input("Enter number of substances: "))
        substances = []
        
        for i in range(substance_count):
            print(f"\nSubstance {i+1}:")
            name = input("Enter substance name or chemical formula: ")
            mass = float(input("Enter mass (g): "))
            
            # Try to automatically calculate molar mass
            use_molar_heat = input("Use molar heat capacity? (y/n): ").lower() == 'y'
            
            if use_molar_heat:
                # Try to calculate molar mass if a chemical formula is given
                try:
                    molar_mass_result = calculate_molar_mass(name)
                    if molar_mass_result['success']:
                        molar_mass = molar_mass_result['molar_mass']
                        print(f"Calculated molar mass: {molar_mass:.4f} g/mol")
                    else:
                        print("Could not automatically calculate molar mass.")
                        molar_mass = float(input("Enter molar mass (g/mol): "))
                except Exception:
                    print("Could not automatically calculate molar mass.")
                    molar_mass = float(input("Enter molar mass (g/mol): "))
                    
                molar_heat_capacity = float(input("Enter molar heat capacity (J/(mol·K)): "))
                # Calculate specific heat
                specific_heat = molar_heat_capacity / molar_mass
                print(f"Calculated specific heat: {specific_heat:.4f} J/(g·K)")
            else:
                specific_heat = float(input("Enter specific heat capacity (J/(g·K)): "))
            
            initial_temp = float(input("Enter initial temperature (°C): "))
            
            substance = {
                'name': name,
                'mass': mass,
                'specific_heat': specific_heat,
                'initial_temp': initial_temp
            }
            
            if use_molar_heat:
                substance['molar_mass'] = molar_mass
                substance['molar_heat_capacity'] = molar_heat_capacity
            
            substances.append(substance)
        
        # Solve problem
        result = solve_mixture_problem(substances)
        
        display_results_header()
        for step in result["steps"]:
            print(step)
        
        print(f"\nFinal equilibrium temperature: {result['final_temp']:.2f} °C")
        print("\nHeat transferred by each substance:")
        for transfer in result["heat_transfers"]:
            print(f"  {transfer['name']}: {transfer['heat']:.2f} J (ΔT = {transfer['delta_t']:.2f} °C)")

    def _handle_clausius_clapeyron_calculations(self):
        """
        Handler function for Clausius-Clapeyron equation calculations.
        """
        print("\n=== Clausius-Clapeyron Equation Calculations ===")
        print("Select the type of calculation:")
        print("[1] Calculate boiling point at a given pressure")
        print("[2] Calculate vapor pressure at a given temperature")
        print("[3] Calculate heat of vaporization from two data points")
        
        choice = input("\nEnter choice (1-3): ").strip()
        
        if choice == "1":
            print("\n--- Calculate Boiling Point at a Given Pressure ---")
            normal_bp = float(input("Enter normal boiling point (°C): "))
            heat_vap = float(input("Enter heat of vaporization (kJ/mol): "))
            init_pressure = float(input("Enter initial pressure (atm, default=1): ") or "1")
            final_pressure = float(input("Enter final pressure (atm): "))
            
            result = calculate_boiling_point_with_pressure(
                normal_bp, heat_vap, init_pressure, final_pressure
            )
            
            display_results_header()
            for step in result["steps"]:
                print(step)
            
            print(f"\nThe boiling point at {final_pressure} atm is {result['new_boiling_point_c']:.2f} °C")
            
        elif choice == "2":
            print("\n--- Calculate Vapor Pressure at a Given Temperature ---")
            normal_bp = float(input("Enter normal boiling point (°C): "))
            heat_vap = float(input("Enter heat of vaporization (kJ/mol): "))
            init_pressure = float(input("Enter initial pressure (atm, default=1): ") or "1")
            final_temp = float(input("Enter temperature (°C): "))
            
            result = calculate_pressure_with_temperature(
                normal_bp, heat_vap, init_pressure, final_temp
            )
            
            display_results_header()
            for step in result["steps"]:
                print(step)
            
            print(f"\nThe vapor pressure at {final_temp} °C is {result['final_pressure']:.6f} atm")
            
        elif choice == "3":
            print("\n--- Calculate Heat of Vaporization from Two Data Points ---")
            temp1 = float(input("Enter first temperature (°C): "))
            pressure1 = float(input("Enter pressure at first temperature (atm): "))
            temp2 = float(input("Enter second temperature (°C): "))
            pressure2 = float(input("Enter pressure at second temperature (atm): "))
            
            result = calculate_heat_of_vaporization(
                temp1, pressure1, temp2, pressure2
            )
            
            display_results_header()
            for step in result["steps"]:
                print(step)
            
            print(f"\nThe heat of vaporization is {result['heat_of_vaporization']:.2f} kJ/mol")
            
        else:
            print("Invalid choice.")