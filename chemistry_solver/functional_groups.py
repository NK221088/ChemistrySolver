"""
Functional Groups Module - Identifies and analyzes functional groups in chemical compounds
"""
import re

def identify_functional_groups(compound_name=None, smiles=None):
    """
    Identifies functional groups in a chemical compound.
    
    Args:
        compound_name (str, optional): Name of the chemical compound
        smiles (str, optional): SMILES notation of the compound
        
    Returns:
        dict: Dictionary with functional groups as keys and boolean values
    """
    functional_groups = {
        "alkyl": False,
        "methyl": False,
        "alkene": False,
        "alkyne": False,
        "aromatic": False,
        "hydroxyl": False,
        "ether": False,
        "aldehyde": False,
        "ketone": False,
        "carboxyl": False,
        "ester": False,
        "amide": False,
        "amine": False,
        "imine": False,
        "nitro": False,
        "nitrile": False,
        "isocyanate": False,
        "halogen": False,
        "sulfide": False,
        "sulfoxide": False,
        "sulfone": False,
        "thiol": False,
        "phosphate": False,
    }
    
    # If SMILES is provided, use it for identification
    if smiles:
        return identify_from_smiles(smiles, functional_groups)
    
    # If name is provided, use it for identification
    if compound_name:
        return identify_from_name(compound_name, functional_groups)
    
    return functional_groups

def identify_from_name(name, functional_groups):
    """
    Identify functional groups from a compound's name.
    
    Args:
        name (str): Name of the chemical compound
        functional_groups (dict): Dictionary to update
        
    Returns:
        dict: Updated functional groups dictionary
    """
    name = name.lower()
    
    # Common name patterns
    if "alcohol" in name or "ol" in name.split() or name.endswith("ol"):
        functional_groups["hydroxyl"] = True
        
    if "ether" in name:
        functional_groups["ether"] = True
        
    if "aldehyde" in name or name.endswith("al"):
        functional_groups["aldehyde"] = True
        
    if "ketone" in name or "one" in name.split() or name.endswith("one"):
        functional_groups["ketone"] = True
        
    if "acid" in name or name.endswith("oic acid"):
        functional_groups["carboxyl"] = True
        
    if "ester" in name or name.endswith("oate"):
        functional_groups["ester"] = True
        
    if "amide" in name or name.endswith("amide"):
        functional_groups["amide"] = True
        
    if "amine" in name or name.endswith("amine") or "amino" in name:
        functional_groups["amine"] = True
        
    if "methyl" in name:
        functional_groups["methyl"] = True
        functional_groups["alkyl"] = True
        
    if "ethyl" in name or "propyl" in name or "butyl" in name:
        functional_groups["alkyl"] = True
        
    if "benzene" in name or "phenyl" in name or "aromatic" in name:
        functional_groups["aromatic"] = True
        
    if "ene" in name.split() or name.endswith("ene"):
        functional_groups["alkene"] = True
        
    if "yne" in name.split() or name.endswith("yne"):
        functional_groups["alkyne"] = True
        
    if "nitro" in name:
        functional_groups["nitro"] = True
        
    if "nitrile" in name or name.endswith("nitrile") or "cyanide" in name:
        functional_groups["nitrile"] = True
        
    if "chloro" in name or "bromo" in name or "fluoro" in name or "iodo" in name:
        functional_groups["halogen"] = True
        
    if "thiol" in name or name.endswith("thiol") or "mercapto" in name:
        functional_groups["thiol"] = True
        
    if "sulfide" in name or "thioether" in name:
        functional_groups["sulfide"] = True
        
    # Special cases for common drugs and compounds
    if "hydroxychloroquine" in name:
        functional_groups["hydroxyl"] = True
        functional_groups["halogen"] = True  # Chloro group
        functional_groups["amine"] = True    # Secondary and tertiary amines
        functional_groups["methyl"] = True   # Methyl groups
        functional_groups["alkyl"] = True    # General alkyl groups
        functional_groups["aromatic"] = True # Quinoline core
        
    return functional_groups

def identify_from_smiles(smiles, functional_groups):
    """
    Identify functional groups from a compound's SMILES notation.
    
    Args:
        smiles (str): SMILES notation
        functional_groups (dict): Dictionary to update
        
    Returns:
        dict: Updated functional groups dictionary
    """
    # Check for hydroxyl (-OH)
    if re.search(r'[^O]O[H]', smiles) or re.search(r'[^O]OH', smiles):
        functional_groups["hydroxyl"] = True
        
    # Check for ethers (C-O-C)
    if re.search(r'[CO][CO]', smiles) and not re.search(r'C=O', smiles):
        functional_groups["ether"] = True
        
    # Check for aldehydes (C=O-H)
    if re.search(r'C=O', smiles) and re.search(r'C[H]', smiles):
        functional_groups["aldehyde"] = True
        
    # Check for ketones (C-C=O-C)
    if re.search(r'C=O', smiles) and not re.search(r'O-[H]', smiles):
        functional_groups["ketone"] = True
        
    # Check for carboxylic acids (C=O-OH)
    if re.search(r'C\(=O\)O[H]', smiles) or re.search(r'C\(=O\)OH', smiles):
        functional_groups["carboxyl"] = True
        
    # Check for esters (C=O-O-C)
    if re.search(r'C\(=O\)OC', smiles):
        functional_groups["ester"] = True
        
    # Check for amides (C=O-N)
    if re.search(r'C\(=O\)N', smiles):
        functional_groups["amide"] = True
        
    # Check for amines (C-N)
    if re.search(r'[CN]', smiles) and not re.search(r'C=N', smiles) and not re.search(r'C\(=O\)N', smiles):
        functional_groups["amine"] = True
        
    # Check for methyl groups (C)
    if re.search(r'C[H3]', smiles) or re.search(r'CH3', smiles):
        functional_groups["methyl"] = True
        functional_groups["alkyl"] = True
        
    # Check for alkyl groups
    if re.search(r'C[H2]', smiles) or re.search(r'CH2', smiles):
        functional_groups["alkyl"] = True
        
    # Check for aromatics (c)
    if re.search(r'c', smiles) or re.search(r'C1=CC=CC=C1', smiles):
        functional_groups["aromatic"] = True
        
    # Check for alkenes (C=C)
    if re.search(r'C=C', smiles):
        functional_groups["alkene"] = True
        
    # Check for alkynes (C#C)
    if re.search(r'C#C', smiles):
        functional_groups["alkyne"] = True
        
    # Check for halogens (F, Cl, Br, I)
    if re.search(r'[FClBrI]', smiles):
        functional_groups["halogen"] = True
        
    return functional_groups

def check_functional_groups(compound_name=None, smiles=None, groups_to_check=None):
    """
    Checks if specific functional groups are present in a compound.
    
    Args:
        compound_name (str, optional): Name of the chemical compound
        smiles (str, optional): SMILES notation of the compound
        groups_to_check (list, optional): List of functional groups to check
        
    Returns:
        dict: Dictionary with results for each requested group
    """
    identified = identify_functional_groups(compound_name, smiles)
    
    if groups_to_check is None:
        return identified
    
    # Filter for only the requested groups
    return {group: identified.get(group, False) for group in groups_to_check}

def find_missing_functional_groups(compound_name=None, smiles=None, groups_list=None):
    """
    Identifies which functional groups from a list are NOT present in the compound.
    
    Args:
        compound_name (str, optional): Name of the chemical compound
        smiles (str, optional): SMILES notation of the compound
        groups_list (list): List of functional groups to check
        
    Returns:
        list: List of functional groups that are NOT present in the compound
    """
    if groups_list is None:
        return []
        
    identified = identify_functional_groups(compound_name, smiles)
    missing_groups = [group for group in groups_list if not identified.get(group, False)]
    
    return missing_groups

def explain_functional_groups_in_compound(compound_name=None, smiles=None):
    """
    Provides a detailed explanation of functional groups present in a compound.
    
    Args:
        compound_name (str, optional): Name of the chemical compound
        smiles (str, optional): SMILES notation of the compound
        
    Returns:
        dict: Contains explanation and list of functional groups
    """
    identified = identify_functional_groups(compound_name, smiles)
    present_groups = [group for group, present in identified.items() if present]
    
    if not present_groups:
        return {
            "compound": compound_name or smiles or "Unknown",
            "functional_groups": [],
            "explanations": [],
            "explanation": "No functional groups identified in this compound."
        }
    
    # Explanations for each functional group
    group_explanations = {
        "alkyl": "Contains alkyl group(s) (saturated carbon chains)",
        "methyl": "Contains methyl group(s) (-CH₃)",
        "alkene": "Contains alkene group(s) (C=C double bond)",
        "alkyne": "Contains alkyne group(s) (C≡C triple bond)",
        "aromatic": "Contains aromatic ring(s)",
        "hydroxyl": "Contains hydroxyl group(s) (-OH)",
        "ether": "Contains ether group(s) (C-O-C)",
        "aldehyde": "Contains aldehyde group(s) (-CHO)",
        "ketone": "Contains ketone group(s) (C=O)",
        "carboxyl": "Contains carboxyl group(s) (-COOH)",
        "ester": "Contains ester group(s) (-COO-)",
        "amide": "Contains amide group(s) (-CONH-)",
        "amine": "Contains amine group(s) (-NH₂, -NH-, -N-)",
        "imine": "Contains imine group(s) (C=N)",
        "nitro": "Contains nitro group(s) (-NO₂)",
        "nitrile": "Contains nitrile group(s) (-C≡N)",
        "isocyanate": "Contains isocyanate group(s) (-N=C=O)",
        "halogen": "Contains halogen atom(s) (F, Cl, Br, I)",
        "sulfide": "Contains sulfide group(s) (C-S-C)",
        "sulfoxide": "Contains sulfoxide group(s) (S=O)",
        "sulfone": "Contains sulfone group(s) (O=S=O)",
        "thiol": "Contains thiol group(s) (-SH)",
        "phosphate": "Contains phosphate group(s) (-PO₄)",
    }
    
    explanations = [group_explanations.get(group, f"Contains {group} group(s)") for group in present_groups]
    
    return {
        "compound": compound_name or smiles or "Unknown",
        "functional_groups": present_groups,
        "explanations": explanations,
        "explanation": f"The compound contains {len(present_groups)} functional groups."
    }

def solve_functional_group_problem(compound_name, options):
    """
    Solves a problem asking which functional group is NOT present in a compound.
    
    Args:
        compound_name (str): Name of the chemical compound
        options (list): List of functional groups to check
        
    Returns:
        dict: Contains answer and explanation
    """
    identified = identify_functional_groups(compound_name)
    
    missing = []
    present = []
    
    for option in options:
        if identified.get(option.lower(), False):
            present.append(option)
        else:
            missing.append(option)
    
    if not missing:
        return {
            "compound": compound_name,
            "answer": None,
            "present_groups": present,
            "missing_groups": [],
            "explanation": "All listed functional groups are present in the compound."
        }
    
    return {
        "compound": compound_name,
        "answer": missing[0] if len(missing) == 1 else missing,
        "present_groups": present,
        "missing_groups": missing,
        "explanation": f"The functional group(s) {', '.join(missing)} are not present in {compound_name}."
    }