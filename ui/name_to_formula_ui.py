"""
Terminal User Interface for Chemical Name to Formula Conversion
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user
from chemistry_solver.name_to_formula import get_formula_from_name, search_compound_properties

class NameToFormulaUI:
    """UI class for chemical name to formula conversion."""
    
    def __init__(self):
        self.title = "CHEMICAL NAME TO FORMULA CONVERTER"
    
    def run(self):
        """Run the name to formula UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-3): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_name_to_formula()
            elif choice == "2":
                self._handle_detailed_search()
            elif choice == "3":
                self._handle_about_pubchem()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the name to formula module menu."""
        menu = """
        [1] Convert chemical name to formula
        [2] Detailed compound search
        [3] About PubChem database
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_name_to_formula(self):
        """Handle conversion of chemical name to formula."""
        print("\n===== CHEMICAL NAME TO FORMULA CONVERTER =====")
        
        try:
            compound_name = input("Enter chemical name: ")
            
            result = get_formula_from_name(compound_name)
            
            display_results_header()
            if result['success']:
                print(f"Name: {result['name']}")
                print(f"Formula: {result['formula']}")
                print(f"IUPAC Name: {result['iupac_name']}")
                print(f"Molecular Weight: {float(result['weight']):.4f} g/mol")
            else:
                print(f"Error: {result['error']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_detailed_search(self):
        """Handle detailed compound search."""
        print("\n===== DETAILED COMPOUND SEARCH =====")
        
        try:
            compound_name = input("Enter chemical name: ")
            
            result = search_compound_properties(compound_name)
            
            display_results_header()
            if result['success']:
                print(f"Name: {result['name']}")
                print(f"Formula: {result['formula']}")
                print(f"IUPAC Name: {result['iupac_name']}")
                print(f"Molecular Weight: {float(result['weight']):.4f} g/mol")
                print(f"PubChem CID: {result['cid']}")
                print(f"Charge: {result['charge']}")
                
                if 'smiles' in result:
                    print(f"SMILES: {result['smiles']}")
                
                if 'synonyms' in result:
                    print("\nCommon Synonyms:")
                    for synonym in result['synonyms']:
                        print(f"- {synonym}")
            else:
                print(f"Error: {result['error']}")
                
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_about_pubchem(self):
        """Display information about PubChem database."""
        print("\n===== ABOUT PUBCHEM DATABASE =====")
        
        info = """
PubChem is an open chemistry database maintained by the National Center for 
Biotechnology Information (NCBI), which is part of the United States National 
Library of Medicine (NLM), an agency of the United States National Institutes 
of Health (NIH).

Key Features:
- Contains information on chemical structures, identifiers, chemical and 
  physical properties, biological activities, and more
- Over 100 million compounds in the database
- Free and open for public access
- Updated daily with new information

This module uses PubChemPy, a Python library that provides a way to interact 
with PubChem's PUG REST API.

Note: Internet connection is required for this functionality to work properly.
        """
        
        print(info)