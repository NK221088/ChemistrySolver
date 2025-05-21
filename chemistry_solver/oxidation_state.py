"""
Oxidation State Calculator Module for Chemistry Problem Solver
"""
import re

def parse_formula(formula):
    """
    Parse a chemical formula into a dictionary of elements and their counts.
    
    Args:
        formula (str): Chemical formula (e.g., "H2O", "Fe2O3")
        
    Returns:
        dict: Dictionary with element symbols as keys and counts as values
    """
    # Regular expression to match elements and their counts
    pattern = r'([A-Z][a-z]*)(\d*)'
    matches = re.findall(pattern, formula)
    
    # Convert matches to dictionary
    elements = {}
    for element, count in matches:
        count = int(count) if count else 1
        if element in elements:
            elements[element] += count
        else:
            elements[element] = count
            
    return elements

def calculate_oxidation_number(compound, element):
    """
    Calculate the oxidation number of an element in a compound.
    
    Args:
        compound (str): Chemical formula (e.g., "CrO2Cl2", "Fe2O3")
        element (str): Element to find oxidation number for (e.g., "Cr", "Fe")
        
    Returns:
        dict: Dictionary with calculation details
    """
    # Parse the compound
    elements = parse_formula(compound)
    
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
    peroxides = ["H2O2", "Na2O2", "BaO2", "K2O2", "Li2O2"]
    superoxides = ["KO2", "NaO2", "RbO2", "CsO2"]
    
    if compound in peroxides:
        # In peroxides, oxygen has -1 oxidation state
        fixed_oxidation['O'] = -1
    elif compound in superoxides:
        # In superoxides, oxygen has -1/2 oxidation state
        fixed_oxidation['O'] = -0.5
    
    # Special case for metal hydrides where H has -1 oxidation state
    metal_hydrides = ["NaH", "KH", "LiH", "CaH2", "MgH2"]
    if compound in metal_hydrides:
        fixed_oxidation['H'] = -1
    
    # Sum of all oxidation numbers must equal the total charge (usually 0 for neutral compounds)
    total_charge = 0  # Assuming a neutral compound
    
    # Calculate the sum of known oxidation numbers
    known_sum = 0
    for elem, count in elements.items():
        if elem != element and elem in fixed_oxidation:
            known_sum += fixed_oxidation[elem] * count
    
    # Calculate the unknown oxidation number
    unknown_oxidation = (total_charge - known_sum) / elements[element]
    
    # Prepare calculation steps for display
    steps = [
        f"1. Identify elements in {compound}: {', '.join([f'{e} (count: {c})' for e, c in elements.items()])}",
        f"2. Assign known oxidation numbers:"
    ]
    
    for elem, count in elements.items():
        if elem != element and elem in fixed_oxidation:
            steps.append(f"   - {elem}: {fixed_oxidation[elem]}")
    
    steps.append(f"3. Total charge of compound: {total_charge}")
    steps.append(f"4. Sum of known oxidation contributions:")
    
    for elem, count in elements.items():
        if elem != element and elem in fixed_oxidation:
            contribution = fixed_oxidation[elem] * count
            steps.append(f"   - {elem}: {fixed_oxidation[elem]} × {count} = {contribution}")
    
    steps.append(f"5. Sum of known oxidation numbers: {known_sum}")
    steps.append(f"6. Calculate oxidation number for {element}:")
    steps.append(f"   - {total_charge} - {known_sum} = {element} oxidation × {elements[element]}")
    steps.append(f"   - {element} oxidation = ({total_charge} - {known_sum}) ÷ {elements[element]} = {unknown_oxidation}")
    
    # Check if the result is an integer or very close to one
    if abs(unknown_oxidation - round(unknown_oxidation)) < 0.01:
        unknown_oxidation = int(round(unknown_oxidation))
        
    return {
        "compound": compound,
        "element": element,
        "oxidation_number": unknown_oxidation,
        "steps": steps
    }