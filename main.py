#!/usr/bin/env python3
"""
Chemistry Problem Solver - Main Entry Point
"""
import sys
from ui.terminal_ui import display_main_menu, get_user_choice
from ui.electron_ui import ElectronConfigUI
from ui.kinetics_ui import KineticsUI
from ui.molar_mass_ui import MolarMassUI
from ui.balancer_ui import BalancerUI
from ui.thermodynamics_ui import ThermodynamicsUI
from ui.acid_base_ui import AcidBaseUI
from ui.name_to_formula_ui import NameToFormulaUI
from ui.stoichiometry_ui import StoichiometryUI
from ui.functional_groups_ui import FunctionalGroupsUI
from ui.colligative_properties_ui import ColligativePropertiesUI
from ui.insoluble_salts_ui import InsolubleSaltsUI
from ui.oxidation_state_ui import OxidationStateUI
from ui.redox_ui import RedoxUI
from ui.gas_concentration_ui import GasConcentrationUI
from ui.ionic_character_ui import IonicCharacterUI

def main():
    """Main application entry point."""
    # Display welcome message
    print("\n" + "="*50)
    print("CHEMISTRY PROBLEM SOLVER".center(50))
    print("="*50 + "\n")
    
    while True:
        # Display the main menu with all available modules
        display_main_menu()
        
        # Get user choice
        choice = get_user_choice()
        
        if choice == "0":
            print("\nExiting program. Goodbye!")
            sys.exit(0)
        elif choice == "1":
            # Launch Electron Configuration module
            electron_ui = ElectronConfigUI()
            electron_ui.run()
        elif choice == "2":
            # Launch Chemical Kinetics module
            kinetics_ui = KineticsUI()
            kinetics_ui.run()
        elif choice == "3":
            # Launch Molar Mass Calculator module
            molar_mass_ui = MolarMassUI()
            molar_mass_ui.run()
        elif choice == "4":
            # Launch Chemical Equation Balancer module
            balancer_ui = BalancerUI()
            balancer_ui.run()
        elif choice == "5":
            # Launch Thermodynamics Calculator module
            thermo_ui = ThermodynamicsUI()
            thermo_ui.run()
        elif choice == "6":
            # Launch Acid-Base Analyzer module
            acid_base_ui = AcidBaseUI()
            acid_base_ui.run()
        elif choice == "7":
            # Launch Name to Formula Converter module
            name_formula_ui = NameToFormulaUI()
            name_formula_ui.run()
        elif choice == "8":
            # Launch Stoichiometry Calculator module
            stoichiometry_ui = StoichiometryUI()
            stoichiometry_ui.run()
        elif choice == "9":
            # Launch Functional Groups Analyzer module
            functional_groups_ui = FunctionalGroupsUI()
            functional_groups_ui.run()
        elif choice == "10":
            # Launch Colligative Properties Calculator module
            colligative_ui = ColligativePropertiesUI()
            colligative_ui.run()
        elif choice == "11":
            # Launch Insoluble Salts & Qualitative Analysis module
            insoluble_salts_ui = InsolubleSaltsUI()
            insoluble_salts_ui.run()
        elif choice == "12":
            # Launch Oxidation State Calculator module
            oxidation_ui = OxidationStateUI()
            oxidation_ui.run()
        elif choice == "13":
            # Launch Redox Reactions Analyzer module
            redox_ui = RedoxUI()
            redox_ui.run()
        elif choice == "14":
            # Launch Gas Concentration Calculator module
            gas_concentration_ui = GasConcentrationUI()
            gas_concentration_ui.run()
        elif choice == "15":
            # Launch Ionic Character Calculator module
            ionic_character_ui = IonicCharacterUI()
            ionic_character_ui.run()
        # Add more modules as needed
        else:
            print("Invalid choice. Please try again.")

if __name__ == "__main__":
    main()