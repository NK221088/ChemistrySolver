"""
Redox Reactions Module for Chemistry Problem Solver
"""

import math
from chemistry_solver.redox_data import (
    STANDARD_REDUCTION_POTENTIALS,
    COMMON_METALS,
    ACTIVITY_SERIES,
    COMMON_REDOX_REACTIONS,
    COMMON_OXIDIZING_AGENTS,
    COMMON_REDUCING_AGENTS,
    COMMON_REDOX_PAIRS,
    get_half_reaction_potential,
    get_redox_pair,
    get_common_redox_reaction,
    get_oxidizing_agent,
    get_reducing_agent
)


def parse_redox_reaction(reaction):
    """
    Parse a redox reaction and identify the oxidation and reduction half-reactions.
    
    Args:
        reaction (str): Chemical equation for a redox reaction
        
    Returns:
        dict: Dictionary containing oxidizing agent, reducing agent, and half-reactions
    """
    # Clean and standardize the reaction
    reaction = reaction.replace("→", "->").replace("⟶", "->")
    
    # First check if this is a common redox reaction in our database
    for name, data in COMMON_REDOX_REACTIONS.items():
        if reaction in data["reaction"] or reaction in data["balanced_equation"]:
            return {
                "oxidizing_agent": get_oxidizing_agent_from_reaction(data),
                "reducing_agent": get_reducing_agent_from_reaction(data),
                "oxidation_half": data["oxidation_half"],
                "reduction_half": data["reduction_half"],
                "electron_transfer": get_electron_transfer_from_half_reaction(data["oxidation_half"]),
                "is_common_reaction": True,
                "common_name": name
            }
    
    # Split into reactants and products
    parts = reaction.split("->")
    if len(parts) != 2:
        raise ValueError("Invalid reaction format. Use '->' or '→' to separate reactants and products.")
        
    reactants = [r.strip() for r in parts[0].split("+")]
    products = [p.strip() for p in parts[1].split("+")]
    
    # Try to identify which species is oxidized (loses electrons) and which is reduced (gains electrons)
    metal_reactants = []
    metal_products = []
    ion_reactants = []
    ion_products = []
    
    # Look for metals and metal ions
    for reactant in reactants:
        if any(metal in reactant and "(s)" in reactant for metal in COMMON_METALS):
            metal_reactants.append(reactant)
        if any(metal in reactant and any(ion_state in reactant for ion_state in ["(aq)", "+"]) for metal in COMMON_METALS):
            ion_reactants.append(reactant)
    
    for product in products:
        if any(metal in product and "(s)" in product for metal in COMMON_METALS):
            metal_products.append(product)
        if any(metal in product and any(ion_state in product for ion_state in ["(aq)", "+"]) for metal in COMMON_METALS):
            ion_products.append(product)
    
    # Identify oxidizing and reducing agents
    oxidizing_agent = None
    reducing_agent = None
    oxidation_half = None
    reduction_half = None
    electron_transfer = None
    
    # Check if this matches any common redox pairs
    for name, data in COMMON_REDOX_PAIRS.items():
        if reaction in data["net_reaction"]:
            return {
                "oxidizing_agent": get_oxidizing_agent_from_pair(data),
                "reducing_agent": get_reducing_agent_from_pair(data),
                "oxidation_half": data["oxidation"],
                "reduction_half": data["reduction"],
                "electron_transfer": get_electron_transfer_from_half_reaction(data["oxidation"]),
                "is_common_pair": True,
                "common_name": name,
                "e_cell": data["e_cell"],
                "favorable": data["favorable"]
            }
    
    # Common case: metal displacement reaction (metal + metal ion -> metal ion + metal)
    if len(metal_reactants) == 1 and len(ion_reactants) == 1 and len(metal_products) == 1 and len(ion_products) == 1:
        # The metal reactant is being oxidized (losing electrons)
        reducing_agent = metal_reactants[0]
        # The ion reactant is being reduced (gaining electrons)
        oxidizing_agent = ion_reactants[0]
        
        # Extract the bare metal symbols (without states)
        reducing_metal = ''.join([c for c in reducing_agent if c.isalpha()])
        oxidizing_ion = ''.join([c for c in oxidizing_agent if c.isalpha()])
        
        # Estimate electron transfer (common case: 2)
        electron_transfer = 2  # Default assumption for most simple metal displacement
        
        # Create half-reactions
        oxidation_half = f"{reducing_metal}(s) -> {reducing_metal}²⁺(aq) + 2e⁻"
        reduction_half = f"{oxidizing_ion}²⁺(aq) + 2e⁻ -> {oxidizing_ion}(s)"
    else:
        # For more complex reactions, we would need more sophisticated analysis
        oxidation_half = "[Oxidation half-reaction could not be automatically determined]"
        reduction_half = "[Reduction half-reaction could not be automatically determined]"
    
    return {
        "oxidizing_agent": oxidizing_agent,
        "reducing_agent": reducing_agent,
        "oxidation_half": oxidation_half,
        "reduction_half": reduction_half,
        "electron_transfer": electron_transfer
    }


def get_oxidizing_agent_from_reaction(reaction_data):
    """Helper function to extract oxidizing agent from reaction data"""
    # In a common redox reaction, the species being reduced is the oxidizing agent
    # It typically appears on the left side of the reduction half-reaction
    reduction_half = reaction_data["reduction_half"]
    parts = reduction_half.split("->")[0].strip().split("+")
    return parts[0].strip()


def get_reducing_agent_from_reaction(reaction_data):
    """Helper function to extract reducing agent from reaction data"""
    # In a common redox reaction, the species being oxidized is the reducing agent
    # It typically appears on the left side of the oxidation half-reaction
    oxidation_half = reaction_data["oxidation_half"]
    parts = oxidation_half.split("->")[0].strip().split("+")
    return parts[0].strip()


def get_oxidizing_agent_from_pair(pair_data):
    """Helper function to extract oxidizing agent from redox pair data"""
    # In a redox pair, the species being reduced is the oxidizing agent
    # It typically appears on the left side of the reduction half-reaction
    reduction_half = pair_data["reduction"]
    parts = reduction_half.split("->")[0].strip().split("+")
    return parts[0].strip()


def get_reducing_agent_from_pair(pair_data):
    """Helper function to extract reducing agent from redox pair data"""
    # In a redox pair, the species being oxidized is the reducing agent
    # It typically appears on the left side of the oxidation half-reaction
    oxidation_half = pair_data["oxidation"]
    parts = oxidation_half.split("->")[0].strip().split("+")
    return parts[0].strip()


def get_electron_transfer_from_half_reaction(half_reaction):
    """Extract the number of electrons transferred from a half-reaction"""
    if "e-" in half_reaction or "e⁻" in half_reaction:
        # Standardize the electron notation
        half_reaction = half_reaction.replace("e⁻", "e-")
        
        # Try to extract the coefficient
        parts = half_reaction.split("e-")
        left_part = parts[0].strip()
        
        # Check if there's a number before e-
        electron_part = left_part.split()[-1].strip()
        if electron_part.isdigit():
            return int(electron_part)
        elif electron_part == "":
            return 1  # Single electron with no coefficient
    
    # Default to 2 if we can't determine
    return 2


def find_standard_reduction_potential(half_reaction):
    """
    Find the standard reduction potential for a given half-reaction.
    
    Args:
        half_reaction (str): Half-reaction string
        
    Returns:
        float or None: Standard reduction potential in volts, or None if not found
    """
    # Use the function from redox_data if possible
    potential = get_half_reaction_potential(half_reaction)
    if potential is not None:
        return potential
    
    # Try direct lookup
    if half_reaction in STANDARD_REDUCTION_POTENTIALS:
        return STANDARD_REDUCTION_POTENTIALS[half_reaction]
    
    # Try normalized lookup (remove spaces)
    normalized = half_reaction.replace(" ", "")
    for key, value in STANDARD_REDUCTION_POTENTIALS.items():
        if key.replace(" ", "") == normalized:
            return value
    
    # Try to extract the metal/ion for simple metal reduction
    # Example: "Cu²⁺(aq) + 2e⁻ → Cu(s)" or equivalent
    if "→" in half_reaction or "->" in half_reaction:
        half_reaction = half_reaction.replace("→", "->")
        parts = half_reaction.split("->")
        if len(parts) == 2:
            left_side = parts[0].strip()
            right_side = parts[1].strip()
            
            # Look for similar patterns in the database
            for key, value in STANDARD_REDUCTION_POTENTIALS.items():
                key_normalized = key.replace(" ", "").replace("→", "->")
                key_parts = key_normalized.split("->")
                if len(key_parts) == 2:
                    key_left = key_parts[0].strip()
                    key_right = key_parts[1].strip()
                    
                    # Check if the core elements match
                    if (any(metal in left_side and metal in key_left for metal in COMMON_METALS) and
                        any(metal in right_side and metal in key_right for metal in COMMON_METALS)):
                        return value
    
    return None


def determine_redox_favorability(reaction):
    """
    Determine if a redox reaction is favorable under standard conditions.
    
    Args:
        reaction (str): Chemical equation for a redox reaction
        
    Returns:
        dict: Results containing favorability information and calculation details
    """
    try:
        # First check if this is a common redox pair with known favorability
        for name, data in COMMON_REDOX_PAIRS.items():
            if reaction in data["net_reaction"]:
                steps = [
                    f"1. Identified reaction: {name} (common redox pair)",
                    f"2. Oxidation half-reaction: {data['oxidation']}",
                    f"3. Reduction half-reaction: {data['reduction']}",
                    f"4. Standard cell potential: E°cell = {data['e_cell']} V",
                    f"5. {'Favorable' if data['favorable'] else 'Not favorable'} reaction (E°cell {'>' if data['favorable'] else '<'} 0)"
                ]
                
                return {
                    "favorable": data["favorable"],
                    "message": "Favorable" if data["favorable"] else "Not favorable",
                    "oxidation_half": data["oxidation"],
                    "reduction_half": data["reduction"],
                    "e_cell": data["e_cell"],
                    "steps": steps
                }
        
        # Try to parse the reaction
        parsed = parse_redox_reaction(reaction)
        
        # Check if we identified a common reaction with parsed data
        if parsed.get("is_common_reaction"):
            common_name = parsed.get("common_name")
            reaction_data = get_common_redox_reaction(common_name)
            
            # Try to calculate cell potential
            reduction_half = reaction_data["reduction_half"]
            oxidation_half = reaction_data["oxidation_half"]
            
            reduction_potential = find_standard_reduction_potential(reduction_half)
            oxidation_potential = find_standard_reduction_potential(
                oxidation_half.replace("->", "←").replace("2e-", "")
            )
            
            if reduction_potential is not None and oxidation_potential is not None:
                e_cell = reduction_potential - oxidation_potential
                is_favorable = e_cell > 0
                
                steps = [
                    f"1. Identified reaction: {common_name} (common redox reaction)",
                    f"2. Oxidation half-reaction: {oxidation_half}",
                    f"3. Reduction half-reaction: {reduction_half}",
                    f"4. Standard reduction potential for reduction: {reduction_potential} V",
                    f"5. Standard reduction potential for oxidation: {oxidation_potential} V",
                    f"6. Standard cell potential: E°cell = {reduction_potential} V - ({oxidation_potential} V) = {e_cell} V",
                    f"7. {'Favorable' if is_favorable else 'Not favorable'} reaction (E°cell {'>' if is_favorable else '<'} 0)"
                ]
                
                return {
                    "favorable": is_favorable,
                    "message": "Favorable" if is_favorable else "Not favorable",
                    "oxidation_half": oxidation_half,
                    "reduction_half": reduction_half,
                    "e_cell": e_cell,
                    "steps": steps
                }
        
        # For simple metal displacement, we can use activity series logic
        is_displacement = parsed["oxidizing_agent"] and parsed["reducing_agent"]
        
        if is_displacement:
            # Extract metal symbols from the agents
            reducing_metal = ''.join([c for c in parsed["reducing_agent"] if c.isalpha()])
            oxidizing_ion = ''.join([c for c in parsed["oxidizing_agent"] if c.isalpha()])
            
            # Find each metal's position in the activity series
            try:
                red_position = ACTIVITY_SERIES.index(reducing_metal)
                oxd_position = ACTIVITY_SERIES.index(oxidizing_ion)
                
                # If the reducing metal is higher in the activity series,
                # the reaction is favorable
                is_favorable = red_position < oxd_position
                
                steps = [
                    f"1. Identified reaction type: Metal displacement reaction",
                    f"2. Reducing agent (loses electrons): {parsed['reducing_agent']}",
                    f"3. Oxidizing agent (gains electrons): {parsed['oxidizing_agent']}",
                    f"4. Position in activity series: {reducing_metal} is at position {red_position}",
                    f"5. Position in activity series: {oxidizing_ion} is at position {oxd_position}",
                    f"6. Using the activity series: {'Favorable' if is_favorable else 'Not favorable'} " +
                    f"(a metal can displace metals lower in the activity series)"
                ]
                
                # Try to calculate cell potential as well
                reduction_potential_oxd = find_standard_reduction_potential(parsed["reduction_half"])
                reduction_potential_red = find_standard_reduction_potential(
                    parsed["oxidation_half"].replace("->", "←").replace("+", "-").replace("2e⁻", "")
                )
                
                if reduction_potential_oxd is not None and reduction_potential_red is not None:
                    cell_potential = reduction_potential_oxd - reduction_potential_red
                    is_favorable_by_potential = cell_potential > 0
                    
                    steps.extend([
                        f"7. Standard reduction potential for {oxidizing_ion}: {reduction_potential_oxd} V",
                        f"8. Standard reduction potential for {reducing_metal}: {reduction_potential_red} V",
                        f"9. Standard cell potential: E°cell = {reduction_potential_oxd} V - ({reduction_potential_red} V) = {cell_potential} V",
                        f"10. Using cell potential: {'Favorable' if is_favorable_by_potential else 'Not favorable'} " +
                        f"(E°cell {'>' if is_favorable_by_potential else '<'} 0)"
                    ])
                    
                    # Check if the two methods agree
                    if is_favorable != is_favorable_by_potential:
                        steps.append("Note: The activity series and cell potential methods give different results. " +
                                    "The cell potential method is generally more accurate.")
                        is_favorable = is_favorable_by_potential
                
                return {
                    "favorable": is_favorable,
                    "message": "Favorable" if is_favorable else "Not favorable",
                    "oxidizing_agent": parsed["oxidizing_agent"],
                    "reducing_agent": parsed["reducing_agent"],
                    "oxidation_half": parsed["oxidation_half"],
                    "reduction_half": parsed["reduction_half"],
                    "steps": steps
                }
                
            except ValueError:
                # One of the metals wasn't found in the activity series
                pass
        
        # Attempt calculation using standard reduction potentials
        # First, we need the proper half-reactions
        oxidation_half = parsed["oxidation_half"]
        reduction_half = parsed["reduction_half"]
        
        # Get the standard reduction potentials
        e_oxidation = find_standard_reduction_potential(
            oxidation_half.replace("->", "←").replace("+", "-").replace("2e⁻", "")
        )
        e_reduction = find_standard_reduction_potential(reduction_half)
        
        if e_oxidation is None or e_reduction is None:
            # If we can't find both potentials, we need to ask the user
            steps = [
                "Could not automatically determine standard reduction potentials.",
                "Please provide the standard reduction potentials for the half-reactions:",
                f"Oxidation half-reaction: {oxidation_half}",
                f"Reduction half-reaction: {reduction_half}"
            ]
            
            return {
                "favorable": None,
                "message": "Insufficient information",
                "oxidation_half": oxidation_half,
                "reduction_half": reduction_half,
                "steps": steps
            }
        
        # Calculate the standard cell potential
        e_cell = e_reduction - e_oxidation
        is_favorable = e_cell > 0
        
        steps = [
            f"1. Oxidation half-reaction: {oxidation_half}",
            f"2. Reduction half-reaction: {reduction_half}",
            f"3. Standard reduction potential for oxidation: {e_oxidation} V",
            f"4. Standard reduction potential for reduction: {e_reduction} V",
            f"5. Standard cell potential: E°cell = {e_reduction} V - ({e_oxidation} V) = {e_cell} V",
            f"6. {'Favorable' if is_favorable else 'Not favorable'} reaction (E°cell {'>' if is_favorable else '<'} 0)"
        ]
        
        return {
            "favorable": is_favorable,
            "message": "Favorable" if is_favorable else "Not favorable",
            "oxidation_half": oxidation_half,
            "reduction_half": reduction_half,
            "e_cell": e_cell,
            "steps": steps
        }
        
    except Exception as e:
        return {
            "favorable": None,
            "message": f"Error: {str(e)}",
            "steps": [f"An error occurred: {str(e)}"]
        }


def determine_favorability_from_potentials(oxidation_potential, reduction_potential):
    """
    Determine if a redox reaction is favorable based on provided potentials.
    
    Args:
        oxidation_potential (float): Standard reduction potential for the species being oxidized
        reduction_potential (float): Standard reduction potential for the species being reduced
        
    Returns:
        dict: Results containing favorability information and calculation details
    """
    # Calculate the standard cell potential
    e_cell = reduction_potential - oxidation_potential
    is_favorable = e_cell > 0
    
    steps = [
        f"1. Standard reduction potential for oxidation: {oxidation_potential} V",
        f"2. Standard reduction potential for reduction: {reduction_potential} V",
        f"3. Standard cell potential: E°cell = {reduction_potential} V - ({oxidation_potential} V) = {e_cell} V",
        f"4. {'Favorable' if is_favorable else 'Not favorable'} reaction (E°cell {'>' if is_favorable else '<'} 0)"
    ]
    
    return {
        "favorable": is_favorable,
        "message": "Favorable" if is_favorable else "Not favorable",
        "e_cell": e_cell,
        "steps": steps
    }


def calculate_nernst_equation(e_cell_standard, n, concentration_ratio, temperature=298.15):
    """
    Calculate the cell potential using the Nernst equation.
    
    Args:
        e_cell_standard (float): Standard cell potential in volts
        n (int): Number of electrons transferred
        concentration_ratio (float): [products]/[reactants] concentration ratio
        temperature (float): Temperature in Kelvin (default: 298.15 K, which is 25°C)
        
    Returns:
        float: Cell potential in volts
    """
    # Constants
    R = 8.314  # J/(mol·K)
    F = 96485  # C/mol
    
    # Calculate the Nernst equation
    e_cell = e_cell_standard - ((R * temperature) / (n * F)) * 2.303 * math.log10(concentration_ratio)
    
    return e_cell


def get_half_reaction_by_element(element, state=None):
    """
    Find a half-reaction containing a specific element in the database.
    
    Args:
        element (str): Chemical element symbol
        state (str, optional): State of matter (e.g., 's', 'aq', 'g')
        
    Returns:
        list: List of matching half-reactions
    """
    results = []
    
    search_term = element
    if state:
        search_term += f"({state})"
    
    for half_reaction, potential in STANDARD_REDUCTION_POTENTIALS.items():
        if search_term in half_reaction:
            results.append({
                "half_reaction": half_reaction,
                "potential": potential
            })
    
    return results