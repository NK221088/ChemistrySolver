"""
Terminal User Interface for Functional Groups Analysis
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user
from chemistry_solver.functional_groups import (
    identify_functional_groups,
    explain_functional_groups_in_compound,
    solve_functional_group_problem,
    check_functional_groups
)

class FunctionalGroupsUI:
    """UI class for functional groups analysis."""
    
    def __init__(self):
        self.title = "FUNCTIONAL GROUPS ANALYZER"
    
    def run(self):
        """Run the functional groups UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-3): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_identify_functional_groups()
            elif choice == "2":
                self._handle_check_functional_groups()
            elif choice == "3":
                self._handle_solve_problem()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the functional groups module menu."""
        menu = """
        [1] Identify all functional groups in a compound
        [2] Check for specific functional groups
        [3] Find which functional group is NOT present
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_identify_functional_groups(self):
        """Handle identification of all functional groups in a compound."""
        print("\n===== FUNCTIONAL GROUPS IDENTIFICATION =====")
        
        try:
            input_type = input("Do you want to enter a compound name or SMILES notation? (name/smiles): ").lower()
            
            compound_name = None
            smiles = None
            
            if input_type == 'smiles':
                smiles = input("Enter SMILES notation: ")
            else:
                compound_name = input("Enter compound name: ")
            
            result = explain_functional_groups_in_compound(compound_name=compound_name, smiles=smiles)
            
            display_results_header()
            print(f"Compound: {result['compound']}")
            
            if result['functional_groups']:
                print("\nIdentified functional groups:")
                print("-" * 50)
                for i, (group, explanation) in enumerate(zip(result['functional_groups'], result['explanations'])):
                    print(f"{i+1}. {group.capitalize()}: {explanation}")
            else:
                print("\nNo functional groups identified in this compound.")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_check_functional_groups(self):
        """Handle checking for specific functional groups in a compound."""
        print("\n===== CHECK SPECIFIC FUNCTIONAL GROUPS =====")
        
        try:
            input_type = input("Do you want to enter a compound name or SMILES notation? (name/smiles): ").lower()
            
            compound_name = None
            smiles = None
            
            if input_type == 'smiles':
                smiles = input("Enter SMILES notation: ")
            else:
                compound_name = input("Enter compound name: ")
            
            print("\nEnter functional groups to check (comma-separated):")
            groups_input = input("e.g., methyl, carboxyl, hydroxyl, amine, halogen: ")
            groups_to_check = [group.strip().lower() for group in groups_input.split(",")]
            
            result = check_functional_groups(compound_name=compound_name, smiles=smiles, groups_to_check=groups_to_check)
            
            display_results_header()
            print(f"Compound: {compound_name or smiles}")
            print("\nFunctional groups present:")
            print("-" * 40)
            
            present_groups = []
            absent_groups = []
            
            for group, is_present in result.items():
                if is_present:
                    present_groups.append(group)
                else:
                    absent_groups.append(group)
            
            if present_groups:
                for group in present_groups:
                    print(f"✓ {group.capitalize()}")
            else:
                print("None of the specified groups are present.")
            
            print("\nFunctional groups absent:")
            print("-" * 40)
            
            if absent_groups:
                for group in absent_groups:
                    print(f"✗ {group.capitalize()}")
            else:
                print("All of the specified groups are present.")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_solve_problem(self):
        """Handle solving which functional group is NOT present in a compound."""
        print("\n===== FUNCTIONAL GROUP PROBLEM SOLVER =====")
        
        try:
            compound_name = input("Enter compound name: ")
            print("\nEnter functional groups to check (comma-separated):")
            options_input = input("e.g., methyl, carboxyl, hydroxyl, amine, halogen: ")
            options = [opt.strip() for opt in options_input.split(",")]
            
            result = solve_functional_group_problem(compound_name, options)
            
            display_results_header()
            print(f"Compound: {result['compound']}")
            
            print("\nPresent functional groups:")
            print("-" * 40)
            for group in result['present_groups']:
                print(f"✓ {group}")
            
            print("\nMissing functional groups:")
            print("-" * 40)
            for group in result['missing_groups']:
                print(f"✗ {group}")
            
            print("\nAnalysis:")
            if len(result['missing_groups']) == 1:
                print(f"The functional group {result['missing_groups'][0]} is NOT present in {result['compound']}.")
            elif len(result['missing_groups']) > 1:
                print(f"Multiple functional groups are not present in {result['compound']}.")
            else:
                print(f"All specified functional groups are present in {result['compound']}.")
                
        except Exception as e:
            print(f"Error: {str(e)}")