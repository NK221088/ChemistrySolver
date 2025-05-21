"""
Acid-Base Analysis Module for Chemistry Problem Solver
"""

def identify_acid_base(compound):
    """
    Identifies if a compound is likely an acid, base, or neutral using Brønsted definition.
    
    Args:
        compound (str): Chemical formula of the compound
    
    Returns:
        dict: Contains classification and explanation
    """
    # Clean the input
    compound = compound.strip()
    
    # Common acids dictionary
    common_acids = {
        "HCl": "Hydrochloric acid - donates H+ in solution",
        "HBr": "Hydrobromic acid - donates H+ in solution",
        "HI": "Hydroiodic acid - donates H+ in solution", 
        "HF": "Hydrofluoric acid - donates H+ in solution",
        "H2SO4": "Sulfuric acid - can donate H+ ions",
        "H2SO3": "Sulfurous acid - can donate H+ ions",
        "HNO3": "Nitric acid - donates H+ in solution",
        "HNO2": "Nitrous acid - donates H+ in solution",
        "H3PO4": "Phosphoric acid - can donate H+ ions",
        "CH3COOH": "Acetic acid - carboxylic acid that donates H+ from -COOH group",
        "HCOOH": "Formic acid - carboxylic acid that donates H+ from -COOH group",
        "H2CO3": "Carbonic acid - can donate H+ ions"
    }
    
    # Common bases dictionary
    common_bases = {
        "NaOH": "Sodium hydroxide - releases OH- which accepts H+",
        "KOH": "Potassium hydroxide - releases OH- which accepts H+",
        "NH3": "Ammonia - accepts H+ with its lone pair",
        "Ca(OH)2": "Calcium hydroxide - releases OH- which accepts H+",
        "Mg(OH)2": "Magnesium hydroxide - releases OH- which accepts H+",
        "NaH": "Sodium hydride - contains H- which acts as a proton acceptor",
        "KH": "Potassium hydride - contains H- which acts as a proton acceptor",
        "LiH": "Lithium hydride - contains H- which acts as a proton acceptor",
        "NaNH2": "Sodium amide - strong base",
        "NaHCO3": "Sodium bicarbonate - weak base"
    }
    
    # Common neutral compounds
    common_neutral = {
        "CH4": "Methane - C-H bonds are not acidic enough to donate protons under normal conditions",
        "C2H6": "Ethane - hydrocarbon with non-acidic C-H bonds",
        "H2O": "Water - amphoteric (can act as both acid and base)",
        "CO2": "Carbon dioxide - forms carbonic acid in water but molecule itself is neutral",
        "N2": "Nitrogen gas - inert molecule"
    }
    
    # First check if it's a common compound we know
    if compound in common_acids:
        return {"classification": "Acid", "explanation": common_acids[compound]}
    
    if compound in common_bases:
        return {"classification": "Base", "explanation": common_bases[compound]}
    
    if compound in common_neutral:
        return {"classification": "Neutral", "explanation": common_neutral[compound]}
    
    # Rules for identifying acids
    # 1. Inorganic acids often start with H
    if compound.startswith("H") and any(char.isupper() for char in compound[1:]):
        return {
            "classification": "Likely Acid", 
            "explanation": f"Compound {compound} starts with H followed by non-metals, suggesting it may donate H+ ions."
        }
    
    # 2. Carboxylic acids end with COOH
    if "COOH" in compound:
        return {
            "classification": "Acid", 
            "explanation": f"Compound {compound} contains a carboxylic acid group (-COOH) which can donate H+."
        }
    
    # Rules for identifying bases
    # 1. Metal hydroxides
    metal_symbols = ["Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba"]
    if any(compound.startswith(metal) for metal in metal_symbols) and "OH" in compound:
        return {
            "classification": "Base", 
            "explanation": f"Compound {compound} appears to be a metal hydroxide, which releases OH- ions that accept H+."
        }
    
    # 2. Metal hydrides
    if any(compound.startswith(metal) and compound[len(metal):] == "H" for metal in metal_symbols):
        return {
            "classification": "Base", 
            "explanation": f"Compound {compound} appears to be a metal hydride, which contains H- ions that accept H+."
        }
    
    # Alcohols (weak acids, often treated as neutral in intro chemistry)
    if "OH" in compound and not any(compound.startswith(metal) for metal in metal_symbols):
        return {
            "classification": "Very Weak Acid / Practically Neutral", 
            "explanation": f"Compound {compound} appears to be an alcohol with -OH group. These can technically donate H+ but are extremely weak acids."
        }
    
    # Default case - need more information
    return {
        "classification": "Unknown", 
        "explanation": f"Cannot confidently classify {compound} without more information about its structure and properties."
    }

def analyze_compound_list(compounds):
    """
    Analyzes a list of compounds and determines if they are all acids.
    
    Args:
        compounds (list): List of chemical formulas
    
    Returns:
        dict: Results for each compound and overall analysis
    """
    results = []
    all_acids = True
    
    for compound in compounds:
        result = identify_acid_base(compound)
        results.append({
            "compound": compound,
            "classification": result["classification"],
            "explanation": result["explanation"]
        })
        
        # Check if this compound is NOT an acid
        if "Acid" not in result["classification"] or "Neutral" in result["classification"]:
            all_acids = False
    
    return {
        "compounds": results,
        "all_acids": all_acids
    }