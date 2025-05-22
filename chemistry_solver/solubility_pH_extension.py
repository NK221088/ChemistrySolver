"""
Extension to Acid-Base Analysis Module for Solubility Equilibrium Problems

This module adds functionality to solve problems involving solubility equilibria,
particularly for calculating pH of solutions with sparingly soluble salts like metal hydroxides.
"""
import math

def solve_hydroxide_salt_pH(formula, ksp):
    """
    Calculate the pH of a water solution saturated with a sparingly soluble metal hydroxide.
    
    Args:
        formula (str): Chemical formula of the metal hydroxide (e.g., "Mg(OH)2", "Ca(OH)2")
        ksp (float): Solubility product constant
        
    Returns:
        dict: Contains calculated values including pH and concentrations
    """
    # Parse the formula to identify the metal and determine the number of OH groups
    metal, n_oh = _parse_hydroxide_formula(formula)
    
    if metal is None or n_oh == 0:
        return {
            "error": f"{formula} does not appear to be a valid metal hydroxide formula."
        }
    
    # Calculate the solubility (s) from Ksp
    # For M(OH)n, Ksp = [M^n+][OH^-]^n
    # If s is the solubility in mol/L, then [M^n+] = s and [OH^-] = n*s
    # This gives Ksp = s * (n*s)^n = s * n^n * s^n = n^n * s^(n+1)
    # Solving for s: s = (Ksp/n^n)^(1/(n+1))
    
    s = (ksp / (n_oh ** n_oh)) ** (1 / (n_oh + 1))
    
    # Calculate [OH-] concentration from the solubility
    oh_concentration = n_oh * s
    
    # Calculate pOH and pH
    poh = -math.log10(oh_concentration)
    ph = 14 - poh  # at 25°C, pH + pOH = 14
    
    return {
        "formula": formula,
        "ksp": ksp,
        "metal": metal,
        "hydroxide_groups": n_oh,
        "solubility": s,
        "metal_ion_concentration": s,
        "hydroxide_concentration": oh_concentration,
        "ph": ph,
        "poh": poh
    }

def _parse_hydroxide_formula(formula):
    """
    Parse a metal hydroxide formula to determine the metal and number of OH groups.
    
    Args:
        formula (str): Chemical formula of the metal hydroxide
        
    Returns:
        tuple: (metal symbol, number of OH groups)
    """
    import re
    
    # Common metal hydroxide patterns
    patterns = [
        # Pattern for M(OH)n where n might be omitted if 1
        r"([A-Z][a-z]*)(?:\(OH\)(\d*))",
        # Pattern for MOH (without parentheses)
        r"([A-Z][a-z]*)OH"
    ]
    
    for pattern in patterns:
        match = re.match(pattern, formula)
        if match:
            metal = match.group(1)
            
            # Get the number of OH groups
            if len(match.groups()) > 1 and match.group(2):
                n_oh = int(match.group(2))
            else:
                n_oh = 1
                
            return metal, n_oh
    
    # If no match found
    return None, 0

def analyze_hydroxide_pH_options(formula, ksp, options):
    """
    Analyzes a multiple choice question about the pH of a metal hydroxide solution
    
    Args:
        formula (str): Chemical formula of the metal hydroxide
        ksp (float): Solubility product constant
        options (list): List of pH value options
        
    Returns:
        dict: Analysis of each option and the correct answer
    """
    # Calculate the actual pH
    result = solve_hydroxide_salt_pH(formula, ksp)
    
    if "error" in result:
        return {"error": result["error"]}
    
    actual_ph = result["ph"]
    
    # Analyze each option
    analysis = {}
    closest_option = None
    smallest_diff = float('inf')
    
    for i, option_str in enumerate(options):
        # Extract the numeric value from the option string
        # Options might be in formats like "~10.4" or "approximately 10.4"
        import re
        ph_match = re.search(r"~?\s*(\d+\.\d+)", option_str)
        
        if ph_match:
            option_ph = float(ph_match.group(1))
            
            # Calculate difference from actual pH
            diff = abs(option_ph - actual_ph)
            
            # Track the closest option
            if diff < smallest_diff:
                smallest_diff = diff
                closest_option = i + 1
            
            is_correct = diff < 0.05  # Consider it correct if within 0.05 pH units
            
            analysis[f"option_{i+1}"] = {
                "ph_value": option_ph,
                "difference": diff,
                "is_correct": is_correct,
                "explanation": f"The calculated pH is {actual_ph:.4f}, which is {'close to' if is_correct else 'significantly different from'} the option value of {option_ph}."
            }
        else:
            analysis[f"option_{i+1}"] = {
                "error": f"Could not extract a pH value from option: {option_str}",
                "is_correct": False
            }
    
    # Find the correct answer
    correct_options = [i+1 for i, data in enumerate(analysis.values()) if "is_correct" in data and data["is_correct"]]
    
    if not correct_options and closest_option:
        # If no exact match found, suggest the closest
        analysis[f"option_{closest_option}"]["suggestion"] = "This is the closest option to the calculated pH."
    
    return {
        "calculated_result": result,
        "options_analysis": analysis,
        "correct_options": correct_options,
        "closest_option": closest_option,
        "summary": f"For a saturated solution of {formula} with Ksp = {ksp}, the calculated pH is {actual_ph:.4f}."
    }

def solve_ksp_from_pH(formula, pH):
    """
    Calculate the solubility product constant (Ksp) from the pH of a saturated metal hydroxide solution.
    
    Args:
        formula (str): Chemical formula of the metal hydroxide
        pH (float): Measured pH of the saturated solution
        
    Returns:
        dict: Contains calculated values including Ksp
    """
    # Parse the formula
    metal, n_oh = _parse_hydroxide_formula(formula)
    
    if metal is None or n_oh == 0:
        return {
            "error": f"{formula} does not appear to be a valid metal hydroxide formula."
        }
    
    # Calculate [OH-] from pH
    poh = 14 - pH  # at 25°C
    oh_concentration = 10 ** (-poh)
    
    # For a metal hydroxide M(OH)n:
    # [M^n+] = s
    # [OH-] = n*s
    # Therefore, s = [OH-]/n
    s = oh_concentration / n_oh
    
    # Calculate Ksp from solubility
    # Ksp = [M^n+][OH-]^n = s * (n*s)^n = s * n^n * s^n = n^n * s^(n+1)
    ksp = (n_oh ** n_oh) * (s ** (n_oh + 1))
    
    return {
        "formula": formula,
        "ph": pH,
        "poh": poh,
        "metal": metal,
        "hydroxide_groups": n_oh,
        "solubility": s,
        "metal_ion_concentration": s,
        "hydroxide_concentration": oh_concentration,
        "ksp": ksp
    }

def analyze_salt_solubility_questions(question_type, formula, **kwargs):
    """
    General function to analyze solubility equilibrium questions
    
    Args:
        question_type (str): Type of question ('pH_from_ksp', 'ksp_from_pH', etc.)
        formula (str): Chemical formula of the salt
        **kwargs: Additional parameters depending on question type
        
    Returns:
        dict: Analysis with calculated values
    """
    if question_type == 'pH_from_ksp':
        # Calculate pH from Ksp
        if 'ksp' not in kwargs:
            return {"error": "Ksp value is required for this calculation."}
        
        ksp = kwargs['ksp']
        
        if 'OH' in formula:
            # For metal hydroxides
            result = solve_hydroxide_salt_pH(formula, ksp)
        else:
            # For other salts (future implementation)
            return {"error": "Currently only metal hydroxide calculations are supported."}
        
        # Handle multiple choice options if provided
        if 'options' in kwargs:
            result['options_analysis'] = analyze_hydroxide_pH_options(formula, ksp, kwargs['options'])
            
        return result
    
    elif question_type == 'ksp_from_pH':
        # Calculate Ksp from pH
        if 'pH' not in kwargs:
            return {"error": "pH value is required for this calculation."}
        
        pH = kwargs['pH']
        
        if 'OH' in formula:
            # For metal hydroxides
            result = solve_ksp_from_pH(formula, pH)
        else:
            # For other salts (future implementation)
            return {"error": "Currently only metal hydroxide calculations are supported."}
            
        return result
        
    else:
        return {"error": f"Unknown question type: {question_type}"}

# Add this function to solve the specific question from the problem
def solve_specific_problem():
    """
    Solve the specific problem: Mg(OH)2 with Ksp = 5.61 × 10^-12
    """
    ksp = 5.61e-12
    formula = "Mg(OH)2"
    
    # Options given in the problem
    options = ["~8.4", "~9.4", "~10.4", "~11.4", "~12.4"]
    
    # Call the general analysis function
    result = analyze_salt_solubility_questions(
        'pH_from_ksp', 
        formula, 
        ksp=ksp,
        options=options
    )
    
    return result