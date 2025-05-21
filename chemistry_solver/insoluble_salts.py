"""
Module for solving qualitative analysis problems involving precipitation reactions.
"""

# Solubility data for common cations with different anions
# True = soluble, False = insoluble
SOLUBILITY_DATA = {
    # Format: 'cation': {'anion': solubility_boolean}
    'Ag+': {
        'Cl-': False,
        'Br-': False,
        'I-': False,
        'OH-': False,
        'SO4^2-': True,  # Silver sulfate is slightly soluble
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Ba^2+': {
        'Cl-': True,
        'Br-': True,
        'I-': True,
        'OH-': True,
        'SO4^2-': False,  # Barium sulfate is insoluble
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Ca^2+': {
        'Cl-': True,
        'Br-': True,
        'I-': True,
        'OH-': False,
        'SO4^2-': True,  # Slightly soluble
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Cu^2+': {
        'Cl-': True,
        'Br-': True,
        'I-': False,
        'OH-': False,
        'SO4^2-': True,
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Fe^2+': {
        'Cl-': True,
        'Br-': True,
        'I-': True,
        'OH-': False,
        'SO4^2-': True,
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Fe^3+': {
        'Cl-': True,
        'Br-': True,
        'I-': True,
        'OH-': False,
        'SO4^2-': True,
        'CO3^3-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Pb^2+': {
        'Cl-': False,  # Slightly soluble in cold water
        'Br-': False,
        'I-': False,
        'OH-': False,
        'SO4^2-': False,  # Lead sulfate is insoluble
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Mg^2+': {
        'Cl-': True,
        'Br-': True,
        'I-': True,
        'OH-': False,
        'SO4^2-': True,
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Hg2^2+': {
        'Cl-': False,
        'Br-': False,
        'I-': False,
        'OH-': False,
        'SO4^2-': True,
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    },
    'Zn^2+': {
        'Cl-': True,
        'Br-': True,
        'I-': True,
        'OH-': False,
        'SO4^2-': True,
        'CO3^2-': False,
        'PO4^3-': False,
        'S^2-': False,
        'NO3-': True,
        'CH3COO-': True
    }
}

# Additional information for reagents
REAGENTS = {
    'H2SO4': {'provides': 'SO4^2-'},
    'NaOH': {'provides': 'OH-'},
    'HCl': {'provides': 'Cl-'},
    'NaCl': {'provides': 'Cl-'},
    'KCl': {'provides': 'Cl-'},
    'AgNO3': {'provides': 'Ag+'},
    'BaCl2': {'provides': 'Ba^2+'}
}

def identify_cation_from_precipitations(candidates, precipitates_with, no_precipitate_with):
    """
    Identify a cation based on which reagents it precipitates with and which it doesn't.
    
    Args:
        candidates (list): List of candidate cations
        precipitates_with (list): List of reagents that cause precipitation
        no_precipitate_with (list): List of reagents that don't cause precipitation
        
    Returns:
        dict: Result containing identified cation and analysis
    """
    results = []
    steps = []
    steps.append("Analyzing precipitation patterns:")
    
    # Check each candidate cation
    for cation in candidates:
        steps.append(f"\nAnalyzing {cation}:")
        matches_all_criteria = True
        
        # Check reagents that should cause precipitation
        for reagent in precipitates_with:
            if reagent not in REAGENTS:
                steps.append(f"  Warning: Reagent {reagent} not in database")
                continue
                
            anion = REAGENTS[reagent]['provides']
            if cation not in SOLUBILITY_DATA:
                steps.append(f"  Warning: Cation {cation} not in database")
                continue
                
            if anion not in SOLUBILITY_DATA[cation]:
                steps.append(f"  Warning: Solubility data for {cation} with {anion} not available")
                continue
                
            is_soluble = SOLUBILITY_DATA[cation][anion]
            precipitates = not is_soluble
            steps.append(f"  {cation} {'precipitates' if precipitates else 'does not precipitate'} with {anion} from {reagent}")
            
            if not precipitates:
                matches_all_criteria = False
                steps.append(f"  ❌ {cation} should precipitate with {reagent}, but doesn't according to data")
                break
        
        # If already failed, skip to next candidate
        if not matches_all_criteria:
            continue
            
        # Check reagents that should NOT cause precipitation
        for reagent in no_precipitate_with:
            if reagent not in REAGENTS:
                steps.append(f"  Warning: Reagent {reagent} not in database")
                continue
                
            anion = REAGENTS[reagent]['provides']
            if cation not in SOLUBILITY_DATA:
                steps.append(f"  Warning: Cation {cation} not in database")
                continue
                
            if anion not in SOLUBILITY_DATA[cation]:
                steps.append(f"  Warning: Solubility data for {cation} with {anion} not available")
                continue
                
            is_soluble = SOLUBILITY_DATA[cation][anion]
            precipitates = not is_soluble
            steps.append(f"  {cation} {'precipitates' if precipitates else 'does not precipitate'} with {anion} from {reagent}")
            
            if precipitates:
                matches_all_criteria = False
                steps.append(f"  ❌ {cation} should NOT precipitate with {reagent}, but does according to data")
                break
        
        # If cation matches all criteria, add to results
        if matches_all_criteria:
            steps.append(f"  ✓ {cation} matches all precipitation patterns")
            results.append(cation)
    
    # Format the conclusion
    if len(results) == 1:
        conclusion = f"The solution must contain {results[0]}. This is the only cation among the candidates that precipitates with {', '.join(precipitates_with)} but not with {', '.join(no_precipitate_with)}."
    elif len(results) > 1:
        conclusion = f"The solution could contain any of these cations: {', '.join(results)}. Further tests would be needed to distinguish between them."
    else:
        conclusion = "No cation matches the given precipitation pattern. There may be an error in the data or the tests."
    
    return {
        "identified_cations": results,
        "steps": steps,
        "conclusion": conclusion
    }

def solve_qualitative_analysis_problem(cation_candidates, precipitates_with=None, no_precipitate_with=None):
    """
    Solves a qualitative analysis problem with given constraints.
    
    Args:
        cation_candidates (list): List of possible cations
        precipitates_with (list): List of reagents that cause precipitation
        no_precipitate_with (list): List of reagents that don't cause precipitation
        
    Returns:
        dict: Results and analysis
    """
    if precipitates_with is None:
        precipitates_with = []
    if no_precipitate_with is None:
        no_precipitate_with = []
        
    # Validate inputs
    for cation in cation_candidates:
        if cation not in SOLUBILITY_DATA:
            return {
                "error": f"Cation {cation} not found in database",
                "available_cations": list(SOLUBILITY_DATA.keys())
            }
    
    for reagent in precipitates_with + no_precipitate_with:
        if reagent not in REAGENTS:
            return {
                "error": f"Reagent {reagent} not found in database",
                "available_reagents": list(REAGENTS.keys())
            }
    
    # Identify cation based on precipitation patterns
    result = identify_cation_from_precipitations(
        cation_candidates, 
        precipitates_with, 
        no_precipitate_with
    )
    
    return result

def analyze_specific_scenario(scenario_id):
    """
    Analyze a specific predefined scenario.
    
    Args:
        scenario_id (str): Identifier for the scenario
        
    Returns:
        dict: Analysis result
    """
    if scenario_id == "W20_8":
        # Scenario from problem W20_8
        return solve_qualitative_analysis_problem(
            cation_candidates=["Ag+", "Ba^2+", "Pb^2+"],
            precipitates_with=["H2SO4"],
            no_precipitate_with=["NaOH"]
        )
    else:
        return {"error": f"Scenario {scenario_id} not found"}

def get_available_cations():
    """
    Get a list of all available cations in the database.
    
    Returns:
        list: List of cation names
    """
    return list(SOLUBILITY_DATA.keys())

def get_available_reagents():
    """
    Get a list of all available reagents in the database.
    
    Returns:
        list: List of reagent names
    """
    return list(REAGENTS.keys())