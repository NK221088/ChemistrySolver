"""
Extended module for solving qualitative analysis problems and Ksp-based solubility calculations.
"""
import math
from .solubility_data import (
    SOLUBILITY_DATA, 
    KSP_DATA, 
    REAGENTS, 
    SCENARIOS,
    get_available_cations,
    get_available_reagents,
    get_available_compounds,
    get_available_scenarios,
    get_scenario_data
)

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
            "available_compounds": get_available_compounds()
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

def identify_absent_cation(candidates, precipitates_with, no_precipitate_with):
    """
    Identify which cation is NOT present based on precipitation tests.
    This is the inverse of the normal identification - we find which cation
    doesn't match the observed pattern.
    
    Args:
        candidates (list): List of candidate cations
        precipitates_with (list): List of reagents that cause precipitation
        no_precipitate_with (list): List of reagents that don't cause precipitation
        
    Returns:
        dict: Result containing absent cations and analysis
    """
    present_cations = []
    absent_cations = []
    steps = []
    steps.append("Determining which cation is NOT present based on precipitation patterns:")
    
    # Check each candidate cation
    for cation in candidates:
        steps.append(f"\nAnalyzing {cation}:")
        matches_pattern = True
        
        # Check reagents that caused precipitation
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
                matches_pattern = False
                steps.append(f"  ❌ {cation} does not precipitate with {reagent}, so it could explain the observed pattern")
        
        # Check reagents that did NOT cause precipitation
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
                matches_pattern = False
                steps.append(f"  ❌ {cation} precipitates with {reagent}, but no precipitation was observed")
        
        # Classify the cation
        if matches_pattern:
            steps.append(f"  ✓ {cation} matches the observed precipitation pattern and could be present")
            present_cations.append(cation)
        else:
            steps.append(f"  ✗ {cation} does NOT match the observed pattern and is likely absent")
            absent_cations.append(cation)
    
    return {
        "present_cations": present_cations,
        "absent_cations": absent_cations,
        "steps": steps
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
                "available_cations": get_available_cations()
            }
    
    for reagent in precipitates_with + no_precipitate_with:
        if reagent not in REAGENTS:
            return {
                "error": f"Reagent {reagent} not found in database",
                "available_reagents": get_available_reagents()
            }
    
    # Identify cation based on precipitation patterns
    result = identify_cation_from_precipitations(
        cation_candidates, 
        precipitates_with, 
        no_precipitate_with
    )
    
    return result

def solve_waste_water_problem():
    """
    Solve the specific waste water problem from the image.
    
    Returns:
        dict: Analysis result showing which cation is NOT present
    """
    scenario_data = get_scenario_data("waste_water")
    if not scenario_data:
        return {"error": "Waste water scenario not found"}
    
    candidates = scenario_data["cation_candidates"]
    precipitates_with = scenario_data["precipitates_with"]
    no_precipitate_with = scenario_data["no_precipitate_with"]
    
    # Use the identify_absent_cation function
    result = identify_absent_cation(candidates, precipitates_with, no_precipitate_with)
    
    # Add specific analysis for this problem
    steps = result["steps"]
    steps.append("\n" + "="*50)
    steps.append("PROBLEM ANALYSIS:")
    steps.append("- HCl causes no precipitation → rules out Ag+ (forms AgCl precipitate)")
    steps.append("- H2SO4 causes precipitation → Sr^2+ forms SrSO4 precipitate")
    steps.append("- Na+ and K+ are always soluble with both HCl and H2SO4")
    steps.append("- Only Sr^2+ explains the observed pattern (no ppt with HCl, ppt with H2SO4)")
    
    conclusion = "Based on the precipitation tests:"
    conclusion += f"\n- Present: {', '.join(result['present_cations'])}"
    conclusion += f"\n- Absent: {', '.join(result['absent_cations'])}"
    conclusion += f"\n\nThe cation that is definitely NOT in the waste water is: {result['absent_cations'][0] if result['absent_cations'] else 'None identified'}"
    
    return {
        "present_cations": result["present_cations"],
        "absent_cations": result["absent_cations"],
        "steps": steps,
        "conclusion": conclusion,
        "answer": result['absent_cations'][0] if result['absent_cations'] else None
    }

def analyze_specific_scenario(scenario_id):
    """
    Analyze a specific predefined scenario.
    
    Args:
        scenario_id (str): Identifier for the scenario
        
    Returns:
        dict: Analysis result
    """
    scenario_data = get_scenario_data(scenario_id)
    if not scenario_data:
        return {
            "error": f"Scenario {scenario_id} not found",
            "available_scenarios": get_available_scenarios()
        }
    
    if scenario_id == "W20_8":
        return solve_qualitative_analysis_problem(
            cation_candidates=scenario_data["cation_candidates"],
            precipitates_with=scenario_data["precipitates_with"],
            no_precipitate_with=scenario_data["no_precipitate_with"]
        )
    elif scenario_id == "S21_11":
        return calculate_solubility_from_ksp(
            scenario_data["compound"], 
            scenario_data["ksp"]
        )
    elif scenario_id == "waste_water":
        return solve_waste_water_problem()
    else:
        return {"error": f"Analysis not implemented for scenario {scenario_id}"}

# Example usage and main execution
if __name__ == "__main__":
    print("Solving Waste Water Analysis Problem:")
    print("="*50)
    
    result = solve_waste_water_problem()
    
    if "error" in result:
        print(f"Error: {result['error']}")
    else:
        # Display analysis steps
        for step in result["steps"]:
            print(step)
        
        print(f"\n{result['conclusion']}")
        
        if result["answer"]:
            print(f"\nFINAL ANSWER: {result['answer']}")