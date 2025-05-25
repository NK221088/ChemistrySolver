"""
Oxidation State Calculator Module for Chemistry Problem Solver
"""
import re

def parse_formula(formula):
    """
    Parse a chemical formula into a dictionary of elements and their counts.
    Now handles charged species like SO4^2-, NH4^+, etc.
    
    Args:
        formula (str): Chemical formula (e.g., "H2O", "Fe2O3", "SO4^2-", "NH4^+")
        
    Returns:
        tuple: (elements_dict, charge)
            - elements_dict: Dictionary with element symbols as keys and counts as values
            - charge: Integer representing the charge of the compound/ion (0 for neutral)
    """
    # Extract charge information first
    charge = 0
    formula_clean = formula
    
    # Look for charge patterns like ^2-, ^+, ^3+, etc.
    charge_pattern = r'\^(\d*)([+-])'
    charge_match = re.search(charge_pattern, formula)
    
    if charge_match:
        charge_num = charge_match.group(1)
        charge_sign = charge_match.group(2)
        
        # If no number is specified, assume 1
        charge_magnitude = int(charge_num) if charge_num else 1
        charge = charge_magnitude if charge_sign == '+' else -charge_magnitude
        
        # Remove charge notation from formula for element parsing
        formula_clean = re.sub(charge_pattern, '', formula)
    
    # Regular expression to match elements and their counts
    pattern = r'([A-Z][a-z]*)(\d*)'
    matches = re.findall(pattern, formula_clean)
    
    # Convert matches to dictionary
    elements = {}
    for element, count in matches:
        count = int(count) if count else 1
        if element in elements:
            elements[element] += count
        else:
            elements[element] = count
            
    return elements, charge

def calculate_oxidation_number(compound, element):
    """
    Calculate the oxidation number of an element in a compound.
    Now properly handles charged species.
    
    Args:
        compound (str): Chemical formula (e.g., "CrO2Cl2", "SO4^2-", "NH4^+")
        element (str): Element to find oxidation number for (e.g., "Cr", "S", "N")
        
    Returns:
        dict: Dictionary with calculation details
    """
    # Parse the compound to get elements and charge
    elements, total_charge = parse_formula(compound)
    
    # Verify that the element exists in the compound
    if element not in elements:
        raise ValueError(f"Element '{element}' not found in compound '{compound}'")
    
    # Known oxidation numbers for common elements
    fixed_oxidation = {
        'O': -2,  # Oxygen usually has -2 (except in peroxides, etc.)
        'F': -1,  # Fluorine always has -1
        'Cl': -1, # Chlorine usually has -1 in compounds
        'Br': -1, # Bromine usually has -1 in compounds
        'I': -1,  # Iodine usually has -1 in compounds
        'Na': 1,  # Sodium usually has +1
        'K': 1,   # Potassium usually has +1
        'Li': 1,  # Lithium usually has +1
        'H': 1    # Hydrogen usually has +1 (except in metal hydrides)
    }
    
    # Handle special cases - more accurately detect peroxides and other special cases
    # Remove charge notation for comparison with special case lists
    compound_neutral = re.sub(r'\^[+-]?\d*[+-]?', '', compound)
    
    peroxides = ["H2O2", "Na2O2", "BaO2", "K2O2", "Li2O2"]
    superoxides = ["KO2", "NaO2", "RbO2", "CsO2"]
    
    if compound_neutral in peroxides:
        # In peroxides, oxygen has -1 oxidation state
        fixed_oxidation['O'] = -1
    elif compound_neutral in superoxides:
        # In superoxides, oxygen has -1/2 oxidation state
        fixed_oxidation['O'] = -0.5
    
    # Special case for metal hydrides where H has -1 oxidation state
    metal_hydrides = ["NaH", "KH", "LiH", "CaH2", "MgH2"]
    if compound_neutral in metal_hydrides:
        fixed_oxidation['H'] = -1
    
    # Calculate the sum of known oxidation numbers
    known_sum = 0
    for elem, count in elements.items():
        if elem != element and elem in fixed_oxidation:
            known_sum += fixed_oxidation[elem] * count
    
    # Calculate the unknown oxidation number
    # Sum of all oxidation numbers must equal the total charge
    unknown_oxidation = (total_charge - known_sum) / elements[element]
    
    # Prepare calculation steps for display
    steps = [
        f"1. Identify elements in {compound}: {', '.join([f'{e} (count: {c})' for e, c in elements.items()])}",
        f"2. Assign known oxidation numbers:"
    ]
    
    for elem, count in elements.items():
        if elem != element and elem in fixed_oxidation:
            steps.append(f"   - {elem}: {fixed_oxidation[elem]}")
    
    charge_description = f"compound/ion: {total_charge}" if total_charge != 0 else f"compound: {total_charge}"
    steps.append(f"3. Total charge of {charge_description}")
    steps.append(f"4. Sum of known oxidation contributions:")
    
    for elem, count in elements.items():
        if elem != element and elem in fixed_oxidation:
            contribution = fixed_oxidation[elem] * count
            steps.append(f"   - {elem}: {fixed_oxidation[elem]} × {count} = {contribution}")
    
    steps.append(f"5. Sum of known oxidation numbers: {known_sum}")
    steps.append(f"6. Calculate oxidation number for {element}:")
    steps.append(f"   - {total_charge} - ({known_sum}) = {element} oxidation × {elements[element]}")
    steps.append(f"   - {element} oxidation = ({total_charge} - ({known_sum})) ÷ {elements[element]} = {unknown_oxidation}")
    
    # Check if the result is an integer or very close to one
    if abs(unknown_oxidation - round(unknown_oxidation)) < 0.01:
        unknown_oxidation = int(round(unknown_oxidation))
        
    return {
        "compound": compound,
        "element": element,
        "oxidation_number": unknown_oxidation,
        "steps": steps
    }