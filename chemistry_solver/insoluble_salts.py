"""
Extended module for solving qualitative analysis problems and Ksp-based solubility calculations.
"""
import math

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

# Ksp values for common slightly soluble compounds
KSP_DATA = {
    'AgCl': {'ksp': 1.6e-10, 'molar_mass': 143.32, 'formula': 'AgCl', 'type': 'AB'},
    'AgBr': {'ksp': 5.0e-13, 'molar_mass': 187.77, 'formula': 'AgBr', 'type': 'AB'},
    'AgI': {'ksp': 8.3e-17, 'molar_mass': 234.77, 'formula': 'AgI', 'type': 'AB'},
    'BaSO4': {'ksp': 1.1e-10, 'molar_mass': 233.39, 'formula': 'BaSO4', 'type': 'AB'},
    'CaCO3': {'ksp': 3.3e-9, 'molar_mass': 100.09, 'formula': 'CaCO3', 'type': 'AB'},
    'PbCl2': {'ksp': 1.9e-4, 'molar_mass': 278.11, 'formula': 'PbCl2', 'type': 'AB2'},
    'PbSO4': {'ksp': 2.5e-8, 'molar_mass': 303.26, 'formula': 'PbSO4', 'type': 'AB'},
    'Mg(OH)2': {'ksp': 5.6e-12, 'molar_mass': 58.32, 'formula': 'Mg(OH)2', 'type': 'AB2'},
    'Ca(OH)2': {'ksp': 5.5e-6, 'molar_mass': 74.09, 'formula': 'Ca(OH)2', 'type': 'AB2'},
    'CaF2': {'ksp': 3.9e-11, 'molar_mass': 78.07, 'formula': 'CaF2', 'type': 'AB2'},
    'Ag2CrO4': {'ksp': 1.1e-12, 'molar_mass': 331.73, 'formula': 'Ag2CrO4', 'type': 'A2B'},
    'Fe(OH)3': {'ksp': 2.8e-39, 'molar_mass': 106.87, 'formula': 'Fe(OH)3', 'type': 'AB3'},
    'Al(OH)3': {'ksp': 3.0e-34, 'molar_mass': 78.00, 'formula': 'Al(OH)3', 'type': 'AB3'}
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

def calculate_solubility_from_ksp(compound, ksp_value=None, temperature=25):
    """
    Calculate the molar and mass solubility of a compound from its Ksp value.
    
    Args:
        compound (str): Chemical formula of the compound
        ksp_value (float): Ksp value (if None, uses value from database)
        temperature (float): Temperature in Celsius (default 25°C)
        
    Returns:
        dict: Results containing molar solubility, mass solubility, and calculation steps
    """
    steps = []
    steps.append(f"Calculating solubility of {compound}")
    
    # Get compound data
    if compound in KSP_DATA and ksp_value is None:
        compound_data = KSP_DATA[compound]
        ksp = compound_data['ksp']
        molar_mass = compound_data['molar_mass']
        compound_type = compound_data['type']
        steps.append(f"Using database Ksp value: {ksp:.2e}")
    elif ksp_value is not None:
        ksp = ksp_value
        if compound in KSP_DATA:
            molar_mass = KSP_DATA[compound]['molar_mass']
            compound_type = KSP_DATA[compound]['type']
        else:
            # Default values if compound not in database
            molar_mass = 100.0  # placeholder
            compound_type = 'AB'  # assume simple 1:1 compound
            steps.append(f"Warning: Compound not in database, using default molar mass of {molar_mass} g/mol")
        steps.append(f"Using provided Ksp value: {ksp:.2e}")
    else:
        return {
            "error": f"Compound {compound} not found in database and no Ksp value provided",
            "available_compounds": list(KSP_DATA.keys())
        }
    
    steps.append(f"Compound type: {compound_type}")
    steps.append(f"Molar mass: {molar_mass} g/mol")
    
    # Calculate molar solubility based on compound type
    if compound_type == 'AB':  # e.g., AgCl → Ag+ + Cl-
        steps.append(f"\nFor AB type compound: {compound} ⇌ A+ + B-")
        steps.append("If solubility = s mol/L, then [A+] = s and [B-] = s")
        steps.append(f"Ksp = [A+][B-] = s × s = s²")
        steps.append(f"s² = {ksp:.2e}")
        molar_solubility = math.sqrt(ksp)
        steps.append(f"s = √({ksp:.2e}) = {molar_solubility:.2e} mol/L")
        
    elif compound_type == 'AB2':  # e.g., PbCl2 → Pb2+ + 2Cl-
        steps.append(f"\nFor AB2 type compound: {compound} ⇌ A2+ + 2B-")
        steps.append("If solubility = s mol/L, then [A2+] = s and [B-] = 2s")
        steps.append(f"Ksp = [A2+][B-]² = s × (2s)² = 4s³")
        steps.append(f"4s³ = {ksp:.2e}")
        steps.append(f"s³ = {ksp/4:.2e}")
        molar_solubility = (ksp / 4) ** (1/3)
        steps.append(f"s = ∛({ksp/4:.2e}) = {molar_solubility:.2e} mol/L")
        
    elif compound_type == 'A2B':  # e.g., Ag2CrO4 → 2Ag+ + CrO4²-
        steps.append(f"\nFor A2B type compound: {compound} ⇌ 2A+ + B2-")
        steps.append("If solubility = s mol/L, then [A+] = 2s and [B2-] = s")
        steps.append(f"Ksp = [A+]²[B2-] = (2s)² × s = 4s³")
        steps.append(f"4s³ = {ksp:.2e}")
        steps.append(f"s³ = {ksp/4:.2e}")
        molar_solubility = (ksp / 4) ** (1/3)
        steps.append(f"s = ∛({ksp/4:.2e}) = {molar_solubility:.2e} mol/L")
        
    elif compound_type == 'AB3':  # e.g., Fe(OH)3 → Fe3+ + 3OH-
        steps.append(f"\nFor AB3 type compound: {compound} ⇌ A3+ + 3B-")
        steps.append("If solubility = s mol/L, then [A3+] = s and [B-] = 3s")
        steps.append(f"Ksp = [A3+][B-]³ = s × (3s)³ = 27s⁴")
        steps.append(f"27s⁴ = {ksp:.2e}")
        steps.append(f"s⁴ = {ksp/27:.2e}")
        molar_solubility = (ksp / 27) ** (1/4)
        steps.append(f"s = ⁴√({ksp/27:.2e}) = {molar_solubility:.2e} mol/L")
        
    else:
        return {"error": f"Unsupported compound type: {compound_type}"}
    
    # Calculate mass solubility
    mass_solubility = molar_solubility * molar_mass
    steps.append(f"\nMass solubility = molar solubility × molar mass")
    steps.append(f"Mass solubility = {molar_solubility:.2e} mol/L × {molar_mass} g/mol")
    steps.append(f"Mass solubility = {mass_solubility:.2e} g/L")
    
    # Convert to different units
    mass_solubility_mg_L = mass_solubility * 1000
    mass_solubility_g_100mL = mass_solubility / 10
    
    return {
        "compound": compound,
        "ksp": ksp,
        "molar_mass": molar_mass,
        "molar_solubility": molar_solubility,
        "mass_solubility_g_L": mass_solubility,
        "mass_solubility_mg_L": mass_solubility_mg_L,
        "mass_solubility_g_100mL": mass_solubility_g_100mL,
        "steps": steps,
        "temperature": temperature
    }

def solve_ksp_problem(compound, ksp_value, problem_type="solubility"):
    """
    Solve various types of Ksp problems.
    
    Args:
        compound (str): Chemical formula
        ksp_value (float): Ksp value
        problem_type (str): Type of problem ("solubility", "common_ion", etc.)
        
    Returns:
        dict: Problem solution
    """
    if problem_type == "solubility":
        return calculate_solubility_from_ksp(compound, ksp_value)
    else:
        return {"error": f"Problem type '{problem_type}' not yet implemented"}

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
    elif scenario_id == "S21_11":
        # Scenario from problem S21_11 - AgCl solubility
        return calculate_solubility_from_ksp("AgCl", 1.6e-10)
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

def get_available_compounds():
    """
    Get a list of all compounds with Ksp data.
    
    Returns:
        list: List of compound names
    """
    return list(KSP_DATA.keys())

# Example usage for the specific problem from the image
if __name__ == "__main__":
    # Solve the AgCl solubility problem from S21_11
    print("Solving AgCl solubility problem (S21_11):")
    print("="*50)
    
    result = calculate_solubility_from_ksp("AgCl", 1.6e-10)
    
    if "error" in result:
        print(f"Error: {result['error']}")
    else:
        # Display calculation steps
        for step in result["steps"]:
            print(step)
        
        print(f"\nFinal Results:")
        print(f"Molar solubility: {result['molar_solubility']:.2e} mol/L")
        print(f"Mass solubility: {result['mass_solubility_g_L']:.2e} g/L")
        print(f"Mass solubility: {result['mass_solubility_mg_L']:.2e} mg/L")
        
        # The answer choices are in scientific notation with g/L
        print(f"\nAnswer: {result['mass_solubility_g_L']:.2e} g/L")
        
        # Check which answer choice this matches
        choices = [1.82e-3, 3.82e-3, 1.82e-5, 3.82e-5, 1.82e-7]
        closest_choice = min(choices, key=lambda x: abs(x - result['mass_solubility_g_L']))
        print(f"Closest answer choice: {closest_choice:.2e} g/L")