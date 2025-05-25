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

from chemistry_solver.oxidation_state import calculate_oxidation_number, parse_formula

def parse_complex_redox_reaction(reaction):
    """
    Parse a complex redox reaction involving polyatomic ions and multiple
    elements changing oxidation states.
    
    Args:
        reaction (str): Chemical equation for a complex redox reaction
        
    Returns:
        dict: Dictionary containing oxidation and reduction half-reactions,
              along with other details if available
    """
    # Clean and standardize the reaction
    reaction = reaction.replace("→", "->").replace("⟶", "->")
    reaction = reaction.replace("^2-", "²⁻").replace("^2+", "²⁺")
    reaction = reaction.replace("^3-", "³⁻").replace("^3+", "³⁺")
    reaction = reaction.replace(" - ", " + ")  # Standardize minus signs as plus signs
    
    # Common complex redox reactions
    # Permanganate reduction (acidic)
    if "MnO4-" in reaction and "Mn2+" in reaction and "H+" in reaction:
        if "SO2" in reaction and "SO4" in reaction:
            # Example: SO2 + MnO4- + H2O -> SO4²⁻ + Mn2+ + H+
            return {
                "oxidation_half": "SO2 + 2H2O → SO4²⁻ + 4H⁺ + 2e⁻",
                "reduction_half": "MnO4⁻ + 8H⁺ + 5e⁻ → Mn²⁺ + 4H2O",
                "balanced_equation": "5SO2 + 2MnO4⁻ + 2H2O → 5SO4²⁻ + 2Mn²⁺ + 4H⁺",
                "oxidizing_agent": "MnO4⁻ (Permanganate ion)",
                "reducing_agent": "SO2 (Sulfur dioxide)",
                "electron_transfer": "10 electrons (5 × 2 = 10)",
                "notes": "In this reaction, sulfur is oxidized from +4 to +6, while manganese is reduced from +7 to +2."
            }
    
    # Dichromate reduction (acidic)
    if "Cr2O7" in reaction and "Cr3+" in reaction:
        if "Fe2+" in reaction and "Fe3+" in reaction:
            # Example: Cr2O7²⁻ + Fe²⁺ + H⁺ → Cr³⁺ + Fe³⁺ + H2O
            return {
                "oxidation_half": "Fe²⁺ → Fe³⁺ + e⁻",
                "reduction_half": "Cr2O7²⁻ + 14H⁺ + 6e⁻ → 2Cr³⁺ + 7H2O",
                "balanced_equation": "Cr2O7²⁻ + 6Fe²⁺ + 14H⁺ → 2Cr³⁺ + 6Fe³⁺ + 7H2O",
                "oxidizing_agent": "Cr2O7²⁻ (Dichromate ion)",
                "reducing_agent": "Fe²⁺ (Iron(II) ion)",
                "electron_transfer": "6 electrons",
                "notes": "This is a common titration reaction. Iron is oxidized from +2 to +3, while chromium is reduced from +6 to +3."
            }
        if "SO2" in reaction or "HSO3" in reaction or "SO3" in reaction:
            # Example: Cr2O7²⁻ + SO2 + H⁺ → Cr³⁺ + SO4²⁻ + H2O
            return {
                "oxidation_half": "SO2 + 2H2O → SO4²⁻ + 4H⁺ + 2e⁻",
                "reduction_half": "Cr2O7²⁻ + 14H⁺ + 6e⁻ → 2Cr³⁺ + 7H2O",
                "balanced_equation": "Cr2O7²⁻ + 3SO2 + 2H⁺ → 2Cr³⁺ + 3SO4²⁻ + H2O",
                "oxidizing_agent": "Cr2O7²⁻ (Dichromate ion)",
                "reducing_agent": "SO2 (Sulfur dioxide)",
                "electron_transfer": "6 electrons",
                "notes": "Sulfur is oxidized from +4 to +6, while chromium is reduced from +6 to +3."
            }
    
    # Disproportionation of hydrogen peroxide
    if "H2O2" in reaction and "H2O" in reaction and "O2" in reaction:
        return {
            "oxidation_half": "H2O2 → O2 + 2H⁺ + 2e⁻",
            "reduction_half": "H2O2 + 2H⁺ + 2e⁻ → 2H2O",
            "balanced_equation": "2H2O2 → 2H2O + O2",
            "oxidizing_agent": "H2O2 (Hydrogen peroxide)",
            "reducing_agent": "H2O2 (Hydrogen peroxide)",
            "electron_transfer": "2 electrons",
            "notes": "This is a disproportionation reaction where hydrogen peroxide acts as both the oxidizing and reducing agent. Oxygen changes from -1 to 0 (oxidation) and from -1 to -2 (reduction)."
        }
    
    # For general undefined complex reactions, provide template for manual completion
    return {
        "oxidation_half": "[Please complete the oxidation half-reaction]",
        "reduction_half": "[Please complete the reduction half-reaction]",
        "balanced_equation": "[Please provide the balanced equation]",
        "oxidizing_agent": "[Please identify the oxidizing agent]",
        "reducing_agent": "[Please identify the reducing agent]",
        "notes": "This complex redox reaction needs manual balancing using the half-reaction method."
    }
    
# Update the existing parse_redox_reaction function to include complex reaction handling
def enhanced_parse_redox_reaction(reaction):
    """
    Enhanced version of parse_redox_reaction that also handles complex reactions.
    
    Args:
        reaction (str): Chemical equation for a redox reaction
        
    Returns:
        dict: Dictionary containing oxidation and reduction half-reactions,
              along with other details if available
    """
    # First try the regular parser
    result = parse_redox_reaction(reaction)
    
    # If we couldn't determine the half-reactions properly, try the complex parser
    if (result["oxidation_half"] == "[Oxidation half-reaction could not be automatically determined]" or
        result["reduction_half"] == "[Reduction half-reaction could not be automatically determined]"):
        complex_result = parse_complex_redox_reaction(reaction)
        
        # Merge the results, with complex parser taking precedence
        result.update(complex_result)
    
    return result


# This function would be used to enhance the determine_redox_favorability function
def enhanced_determine_redox_favorability(reaction):
    """
    Enhanced version of determine_redox_favorability that also handles complex reactions.
    
    Args:
        reaction (str): Chemical equation for a redox reaction
        
    Returns:
        dict: Results containing favorability information and calculation details
    """
    try:
        # First check if this is a complex reaction
        complex_result = parse_complex_redox_reaction(reaction)
        
        # If we have specific half-reactions for this complex reaction
        if not complex_result["oxidation_half"].startswith("[Please"):
            # Create steps for the analysis
            steps = [
                f"1. Identified reaction: Complex redox reaction",
                f"2. Oxidation half-reaction: {complex_result['oxidation_half']}",
                f"3. Reduction half-reaction: {complex_result['reduction_half']}",
                f"4. Balanced equation: {complex_result['balanced_equation']}",
                f"5. Oxidizing agent: {complex_result['oxidizing_agent']}",
                f"6. Reducing agent: {complex_result['reducing_agent']}",
            ]
            
            if "notes" in complex_result:
                steps.append(f"7. Notes: {complex_result['notes']}")
            
            # For known reactions, we can determine favorability (generally all examples are favorable)
            return {
                "favorable": True,
                "message": "Favorable",
                "oxidation_half": complex_result["oxidation_half"],
                "reduction_half": complex_result["reduction_half"],
                "steps": steps,
                "oxidizing_agent": complex_result.get("oxidizing_agent"),
                "reducing_agent": complex_result.get("reducing_agent"),
                "balanced_equation": complex_result.get("balanced_equation"),
                "notes": complex_result.get("notes")
            }
        
        # If complex reaction parsing didn't give specific results, fall back to the original function
        return determine_redox_favorability(reaction)
        
    except Exception as e:
        return {
            "favorable": None,
            "message": f"Error: {str(e)}",
            "steps": [f"An error occurred: {str(e)}"]
        }

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

def identify_redox_from_oxidation_states(reaction):
    """
    Identify oxidation and reduction by analyzing oxidation state changes.
    
    Args:
        reaction (str): Chemical equation for a redox reaction
        
    Returns:
        dict: Dictionary containing oxidation state analysis and redox identification
    """
    # Clean and standardize the reaction
    reaction = reaction.replace("→", "->").replace("⟶", "->")
    
    # Split into reactants and products
    parts = reaction.split("->")
    if len(parts) != 2:
        raise ValueError("Invalid reaction format. Use '->' or '→' to separate reactants and products.")
        
    reactants = [r.strip() for r in parts[0].split("+")]
    products = [p.strip() for p in parts[1].split("+")]
    
    # Track oxidation state changes
    oxidation_changes = {}
    element_compounds = {}
    
    # Analyze reactants
    for compound in reactants:
        # Remove state notations like (s), (aq), (g), (l)
        clean_compound = compound.replace("(s)", "").replace("(aq)", "").replace("(g)", "").replace("(l)", "").strip()
        
        try:
            elements, charge = parse_formula(clean_compound)
            for element in elements:
                try:
                    ox_result = calculate_oxidation_number(clean_compound, element)
                    ox_state = ox_result["oxidation_number"]
                    
                    if element not in oxidation_changes:
                        oxidation_changes[element] = {"reactant": [], "product": []}
                        element_compounds[element] = {"reactant": [], "product": []}
                    
                    oxidation_changes[element]["reactant"].append(ox_state)
                    element_compounds[element]["reactant"].append(compound)
                except (ValueError, KeyError):
                    continue
        except Exception:
            continue
    
    # Analyze products
    for compound in products:
        # Remove state notations
        clean_compound = compound.replace("(s)", "").replace("(aq)", "").replace("(g)", "").replace("(l)", "").strip()
        
        try:
            elements, charge = parse_formula(clean_compound)
            for element in elements:
                try:
                    ox_result = calculate_oxidation_number(clean_compound, element)
                    ox_state = ox_result["oxidation_number"]
                    
                    if element not in oxidation_changes:
                        oxidation_changes[element] = {"reactant": [], "product": []}
                        element_compounds[element] = {"reactant": [], "product": []}
                    
                    oxidation_changes[element]["product"].append(ox_state)
                    element_compounds[element]["product"].append(compound)
                except (ValueError, KeyError):
                    continue
        except Exception:
            continue
    
    # Identify which elements are oxidized and reduced
    oxidized_elements = []
    reduced_elements = []
    oxidation_analysis = []
    
    for element, states in oxidation_changes.items():
        if states["reactant"] and states["product"]:
            reactant_state = states["reactant"][0]  # Take first occurrence
            product_state = states["product"][0]    # Take first occurrence
            
            if reactant_state < product_state:
                # Oxidation (increase in oxidation state)
                oxidized_elements.append({
                    "element": element,
                    "from_state": reactant_state,
                    "to_state": product_state,
                    "change": product_state - reactant_state,
                    "reactant_compound": element_compounds[element]["reactant"][0],
                    "product_compound": element_compounds[element]["product"][0]
                })
                oxidation_analysis.append(f"{element}: {reactant_state} → {product_state} (oxidized, loses {product_state - reactant_state} electrons)")
                
            elif reactant_state > product_state:
                # Reduction (decrease in oxidation state)
                reduced_elements.append({
                    "element": element,
                    "from_state": reactant_state,
                    "to_state": product_state,
                    "change": reactant_state - product_state,
                    "reactant_compound": element_compounds[element]["reactant"][0],
                    "product_compound": element_compounds[element]["product"][0]
                })
                oxidation_analysis.append(f"{element}: {reactant_state} → {product_state} (reduced, gains {reactant_state - product_state} electrons)")
            else:
                oxidation_analysis.append(f"{element}: {reactant_state} → {product_state} (no change)")
    
    return {
        "oxidized_elements": oxidized_elements,
        "reduced_elements": reduced_elements,
        "oxidation_analysis": oxidation_analysis,
        "is_redox": len(oxidized_elements) > 0 and len(reduced_elements) > 0
    }

# Add this enhanced version of parse_redox_reaction that uses oxidation states
def enhanced_parse_redox_reaction_with_oxidation_states(reaction):
    """
    Enhanced version that uses oxidation state analysis to identify redox reactions.
    
    Args:
        reaction (str): Chemical equation for a redox reaction
        
    Returns:
        dict: Enhanced dictionary with oxidation state analysis
    """
    # First get the original parsing result
    original_result = parse_redox_reaction(reaction)
    
    # Add oxidation state analysis
    try:
        ox_analysis = identify_redox_from_oxidation_states(reaction)
        
        # Enhance the result with oxidation state information
        enhanced_result = original_result.copy()
        enhanced_result.update({
            "oxidation_state_analysis": ox_analysis["oxidation_analysis"],
            "oxidized_elements": ox_analysis["oxidized_elements"],
            "reduced_elements": ox_analysis["reduced_elements"],
            "confirmed_redox": ox_analysis["is_redox"]
        })
        
        # If original parsing couldn't identify agents, try using oxidation state analysis
        if (not original_result.get("oxidizing_agent") or not original_result.get("reducing_agent")) and ox_analysis["is_redox"]:
            if ox_analysis["oxidized_elements"]:
                enhanced_result["reducing_agent"] = ox_analysis["oxidized_elements"][0]["reactant_compound"]
            if ox_analysis["reduced_elements"]:
                enhanced_result["oxidizing_agent"] = ox_analysis["reduced_elements"][0]["reactant_compound"]
        
        return enhanced_result
        
    except Exception as e:
        # If oxidation state analysis fails, return original result with error note
        enhanced_result = original_result.copy()
        enhanced_result["oxidation_state_error"] = str(e)
        return enhanced_result

# Add this function to validate if a reaction is actually a redox reaction
def validate_redox_reaction(reaction):
    """
    Validate if a given reaction is actually a redox reaction by checking oxidation state changes.
    
    Args:
        reaction (str): Chemical equation to validate
        
    Returns:
        dict: Validation results with detailed analysis
    """
    try:
        ox_analysis = identify_redox_from_oxidation_states(reaction)
        
        validation_steps = [
            f"1. Analyzing oxidation states in reaction: {reaction}",
            "2. Oxidation state changes found:"
        ]
        
        if ox_analysis["oxidation_analysis"]:
            validation_steps.extend([f"   - {analysis}" for analysis in ox_analysis["oxidation_analysis"]])
        else:
            validation_steps.append("   - No oxidation state changes detected")
        
        if ox_analysis["is_redox"]:
            validation_steps.extend([
                "3. This IS a redox reaction:",
                f"   - Elements oxidized: {', '.join([elem['element'] for elem in ox_analysis['oxidized_elements']])}",
                f"   - Elements reduced: {', '.join([elem['element'] for elem in ox_analysis['reduced_elements']])}"
            ])
        else:
            validation_steps.append("3. This is NOT a redox reaction (no electron transfer detected)")
        
        return {
            "is_redox": ox_analysis["is_redox"],
            "oxidized_elements": ox_analysis["oxidized_elements"],
            "reduced_elements": ox_analysis["reduced_elements"],
            "validation_steps": validation_steps,
            "message": "Valid redox reaction" if ox_analysis["is_redox"] else "Not a redox reaction"
        }
        
    except Exception as e:
        return {
            "is_redox": None,
            "message": f"Error during validation: {str(e)}",
            "validation_steps": [f"Error occurred: {str(e)}"]
        }

# Replace the existing enhanced_determine_redox_favorability function with this version
def enhanced_determine_redox_favorability_with_oxidation_states(reaction):
    """
    Enhanced version that includes oxidation state analysis in favorability determination.
    
    Args:
        reaction (str): Chemical equation for a redox reaction
        
    Returns:
        dict: Results with oxidation state analysis and favorability information
    """
    try:
        # First validate that this is actually a redox reaction
        validation = validate_redox_reaction(reaction)
        
        if not validation["is_redox"]:
            return {
                "favorable": None,
                "message": "Not a redox reaction",
                "validation": validation,
                "steps": validation["validation_steps"] + ["Cannot determine favorability for non-redox reactions."]
            }
        
        # Get enhanced parsing with oxidation states
        enhanced_parsing = enhanced_parse_redox_reaction_with_oxidation_states(reaction)
        
        # Try to get favorability using existing methods
        favorability_result = enhanced_determine_redox_favorability(reaction)
        
        # Combine results
        combined_result = favorability_result.copy()
        combined_result.update({
            "validation": validation,
            "oxidation_state_analysis": enhanced_parsing.get("oxidation_state_analysis", []),
            "oxidized_elements": enhanced_parsing.get("oxidized_elements", []),
            "reduced_elements": enhanced_parsing.get("reduced_elements", [])
        })
        
        # Enhance steps with oxidation state information
        if "steps" in combined_result:
            oxidation_steps = [
                "Oxidation State Analysis:",
                *[f"   - {analysis}" for analysis in enhanced_parsing.get("oxidation_state_analysis", [])]
            ]
            combined_result["steps"] = oxidation_steps + [""] + combined_result["steps"]
        
        return combined_result
        
    except Exception as e:
        return {
            "favorable": None,
            "message": f"Error: {str(e)}",
            "steps": [f"An error occurred: {str(e)}"]
        }

# Add this function to balance redox equations using oxidation state information
def balance_redox_equation_with_oxidation_states(reaction):
    """
    Attempt to balance a redox equation using oxidation state information.
    
    Args:
        reaction (str): Unbalanced chemical equation
        
    Returns:
        dict: Balancing information and suggested balanced equation
    """
    try:
        ox_analysis = identify_redox_from_oxidation_states(reaction)
        
        balancing_steps = [
            f"1. Analyzing unbalanced equation: {reaction}",
            "2. Oxidation state changes:"
        ]
        
        if ox_analysis["oxidation_analysis"]:
            balancing_steps.extend([f"   - {analysis}" for analysis in ox_analysis["oxidation_analysis"]])
        
        if not ox_analysis["is_redox"]:
            return {
                "balanced_equation": reaction,
                "message": "Not a redox reaction - no balancing needed",
                "steps": balancing_steps + ["3. No electron transfer detected"]
            }
        
        # Calculate electron transfer
        total_electrons_lost = sum(elem["change"] for elem in ox_analysis["oxidized_elements"])
        total_electrons_gained = sum(elem["change"] for elem in ox_analysis["reduced_elements"])
        
        balancing_steps.extend([
            f"3. Total electrons lost: {total_electrons_lost}",
            f"4. Total electrons gained: {total_electrons_gained}",
            "5. For a balanced equation, electrons lost must equal electrons gained"
        ])
        
        if total_electrons_lost != total_electrons_gained:
            balancing_steps.append("6. Equation needs balancing - use half-reaction method or algebraic method")
        else:
            balancing_steps.append("6. Electron transfer is already balanced")
        
        return {
            "balanced_equation": "[Use half-reaction method to balance]",
            "electrons_lost": total_electrons_lost,
            "electrons_gained": total_electrons_gained,
            "needs_balancing": total_electrons_lost != total_electrons_gained,
            "steps": balancing_steps,
            "oxidized_elements": ox_analysis["oxidized_elements"],
            "reduced_elements": ox_analysis["reduced_elements"]
        }
        
    except Exception as e:
        return {
            "balanced_equation": reaction,
            "message": f"Error during balancing analysis: {str(e)}",
            "steps": [f"Error occurred: {str(e)}"]
        }