"""
Chemical Name to Formula Conversion Module for Chemistry Problem Solver
"""
import pubchempy as pcp

def get_formula_from_name(compound_name):
    """
    Retrieves chemical formula for a given compound name using PubChemPy
    
    Args:
        compound_name (str): Name of the chemical compound
        
    Returns:
        dict: Dictionary containing formula and other information if successful,
              or error message if not successful
    """
    try:
        compounds = pcp.get_compounds(compound_name, 'name')
        if compounds:
            compound = compounds[0]  # Take the first match
            return {
                'success': True,
                'name': compound_name,
                'formula': compound.molecular_formula,
                'iupac_name': compound.iupac_name,
                'weight': compound.molecular_weight
            }
        else:
            return {
                'success': False,
                'error': f"Could not find compound: {compound_name}"
            }
    except Exception as e:
        return {
            'success': False,
            'error': f"Error: {str(e)}"
        }

def search_compound_properties(compound_name):
    """
    Retrieves detailed properties for a given compound name using PubChemPy
    
    Args:
        compound_name (str): Name of the chemical compound
        
    Returns:
        dict: Dictionary containing extended properties if successful,
              or error message if not successful
    """
    try:
        compounds = pcp.get_compounds(compound_name, 'name')
        if compounds:
            compound = compounds[0]  # Take the first match
            
            # Get basic properties
            properties = {
                'success': True,
                'name': compound_name,
                'formula': compound.molecular_formula,
                'iupac_name': compound.iupac_name,
                'weight': compound.molecular_weight,
                'cid': compound.cid,
                'charge': compound.charge
            }
            
            # Add synonyms if available
            if hasattr(compound, 'synonyms') and compound.synonyms:
                # Limit to first 5 synonyms to keep it manageable
                properties['synonyms'] = compound.synonyms[:5]
            
            # Add canonical SMILES if available
            if hasattr(compound, 'canonical_smiles') and compound.canonical_smiles:
                properties['smiles'] = compound.canonical_smiles
                
            return properties
        else:
            return {
                'success': False,
                'error': f"Could not find compound: {compound_name}"
            }
    except Exception as e:
        return {
            'success': False,
            'error': f"Error: {str(e)}"
        }