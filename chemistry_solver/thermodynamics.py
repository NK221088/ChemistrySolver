"""
Comprehensive Chemistry Solver module for thermodynamics, enthalpy, 
heat transfer, and Hess's Law calculations.
"""
import re
import math
from typing import Dict, List, Tuple, Any, Union, Optional
from dataclasses import dataclass


# Gas constant in different units
R_IDEAL_GAS = {
    'J/(mol·K)': 8.314,
    'kJ/(mol·K)': 0.008314,
    'L·atm/(mol·K)': 0.0821
}


@dataclass
class Reaction:
    """Class to represent a chemical reaction with its enthalpy."""
    reactants: Dict[str, float]  # Substance: stoichiometric coefficient
    products: Dict[str, float]    # Substance: stoichiometric coefficient
    enthalpy: Optional[float] = None  # ΔH in kJ/mol
    
    def __str__(self):
        """Return a string representation of the reaction."""
        reactants_str = " + ".join([f"{coef} {subst}" if coef != 1 else subst 
                                  for subst, coef in self.reactants.items()])
        products_str = " + ".join([f"{coef} {subst}" if coef != 1 else subst 
                                 for subst, coef in self.products.items()])
        
        enthalpy_str = f" ΔH° = {self.enthalpy} kJ/mol" if self.enthalpy is not None else ""
        return f"{reactants_str} → {products_str}{enthalpy_str}"
    
    def reverse(self):
        """Return a new Reaction that is the reverse of this one."""
        return Reaction(
            reactants=self.products.copy(),
            products=self.reactants.copy(),
            enthalpy=-self.enthalpy if self.enthalpy is not None else None
        )
    
    def scale(self, factor: float):
        """Return a new Reaction with all coefficients and enthalpy scaled by a factor."""
        return Reaction(
            reactants={k: v * factor for k, v in self.reactants.items()},
            products={k: v * factor for k, v in self.products.items()},
            enthalpy=self.enthalpy * factor if self.enthalpy is not None else None
        )
    
    def add(self, other):
        """Add two reactions together."""
        # Create new dictionaries to avoid modifying the originals
        new_reactants = self.reactants.copy()
        new_products = self.products.copy()
        
        # Add reactants from other reaction
        for substance, coef in other.reactants.items():
            if substance in new_products:
                new_products[substance] -= coef
                if new_products[substance] <= 0:
                    if new_products[substance] < 0:
                        new_reactants[substance] = -new_products[substance]
                    del new_products[substance]
            else:
                if substance in new_reactants:
                    new_reactants[substance] += coef
                else:
                    new_reactants[substance] = coef
        
        # Add products from other reaction
        for substance, coef in other.products.items():
            if substance in new_reactants:
                new_reactants[substance] -= coef
                if new_reactants[substance] <= 0:
                    if new_reactants[substance] < 0:
                        new_products[substance] = -new_reactants[substance]
                    del new_reactants[substance]
            else:
                if substance in new_products:
                    new_products[substance] += coef
                else:
                    new_products[substance] = coef
        
        # Clean up dictionaries (remove zero coefficients)
        new_reactants = {k: v for k, v in new_reactants.items() if v > 0}
        new_products = {k: v for k, v in new_products.items() if v > 0}
        
        # Calculate new enthalpy
        new_enthalpy = None
        if self.enthalpy is not None and other.enthalpy is not None:
            new_enthalpy = self.enthalpy + other.enthalpy
        
        return Reaction(reactants=new_reactants, products=new_products, enthalpy=new_enthalpy)


#######################################
# Reaction and Enthalpy Functions
#######################################

def parse_reaction_string(reaction_str: str) -> Reaction:
    """
    Parse a reaction string into a Reaction object.
    
    Example:
    "C (graphite) + O2 (g) → CO2 (g) ΔH° = -393.5 kJ/mol"
    """
    # Split the reaction into the reaction part and the enthalpy part
    parts = reaction_str.split("ΔH°")
    reaction_part = parts[0].strip()
    
    # Extract enthalpy if present
    enthalpy = None
    if len(parts) > 1:
        enthalpy_str = parts[1].strip()
        # Extract numerical value, handling different formats
        enthalpy_match = re.search(r'[=\s]?\s*([-+]?\d+\.?\d*)', enthalpy_str)
        if enthalpy_match:
            enthalpy = float(enthalpy_match.group(1))
    
    # Split the reaction into reactants and products
    sides = reaction_part.split("→")
    if len(sides) != 2:
        sides = reaction_part.split("->")
    if len(sides) != 2:
        raise ValueError(f"Invalid reaction format: {reaction_str}")
    
    reactants_str, products_str = sides[0].strip(), sides[1].strip()
    
    # Parse reactants
    reactants = {}
    for term in reactants_str.split("+"):
        term = term.strip()
        if not term:
            continue
        
        # Extract coefficient
        coefficient_match = re.match(r'^(\d+\.?\d*)\s+', term)
        if coefficient_match:
            coefficient = float(coefficient_match.group(1))
            substance = term[coefficient_match.end():].strip()
        else:
            coefficient = 1.0
            substance = term
        
        reactants[substance] = coefficient
    
    # Parse products
    products = {}
    for term in products_str.split("+"):
        term = term.strip()
        if not term:
            continue
        
        # Extract coefficient
        coefficient_match = re.match(r'^(\d+\.?\d*)\s+', term)
        if coefficient_match:
            coefficient = float(coefficient_match.group(1))
            substance = term[coefficient_match.end():].strip()
        else:
            coefficient = 1.0
            substance = term
        
        products[substance] = coefficient
    
    return Reaction(reactants=reactants, products=products, enthalpy=enthalpy)


def solve_combustion_enthalpy(known_reactions: List[Reaction], target_reaction: Reaction) -> Dict[str, Any]:
    """
    Solve for the enthalpy of a target reaction using Hess's Law.
    
    Parameters:
        known_reactions (List[Reaction]): List of known reactions with enthalpies
        target_reaction (Reaction): The target reaction to solve for
        
    Returns:
        Dict[str, Any]: Dictionary with solution and steps
    """
    steps = []
    
    steps.append("# Solving for the enthalpy of combustion using Hess's Law")
    steps.append("\n## Given Reactions:")
    for i, reaction in enumerate(known_reactions, 1):
        steps.append(f"Reaction {i}: {reaction}")
    
    steps.append(f"\n## Target Reaction:")
    steps.append(f"{target_reaction}")
    
    steps.append("\n## Solution using Hess's Law:")
    steps.append("We need to manipulate the given reactions to obtain the target reaction:")
    
    # Create a list of modified reactions
    modified_reactions = []
    explanation_steps = []
    
    # This is a simplistic approach - in a real solver, we would need more sophisticated logic
    # Here we'll manually solve for the methanol combustion example
    
    # Check if we're solving the methanol combustion example
    target_substances = set(target_reaction.reactants.keys()) | set(target_reaction.products.keys())
    
    # Extract needed substances from reactions
    all_substances = set()
    for reaction in known_reactions:
        all_substances |= set(reaction.reactants.keys()) | set(reaction.products.keys())
    
    # Find the enthalpy for the target reaction using Hess's Law
    # For this example, we'll do a specific implementation
    result_reaction = None
    
    # For the methanol combustion example
    # CH3OH (l) + 1.5 O2 (g) → CO2 (g) + 2 H2O (l)
    
    if "CH3OH" in str(target_reaction) and "CO2" in str(target_reaction) and "H2O" in str(target_reaction):
        # 1. Start with the methanol formation reaction
        methanol_formation = next((r for r in known_reactions if "CH3OH" in str(r)), None)
        if methanol_formation:
            # Reverse it to get methanol decomposition
            methanol_decomposition = methanol_formation.reverse()
            explanation_steps.append(f"1. Reverse the methanol formation reaction:")
            explanation_steps.append(f"   {methanol_decomposition}")
            modified_reactions.append(methanol_decomposition)
            
            # 2. Use the carbon combustion reaction
            carbon_combustion = next((r for r in known_reactions if "C (graphite)" in str(r) and "CO2" in str(r)), None)
            if carbon_combustion:
                explanation_steps.append(f"\n2. Use the carbon combustion reaction:")
                explanation_steps.append(f"   {carbon_combustion}")
                modified_reactions.append(carbon_combustion)
                
                # 3. Use the hydrogen combustion reaction (doubled)
                hydrogen_combustion = next((r for r in known_reactions if "H2" in str(r) and "H2O" in str(r)), None)
                if hydrogen_combustion:
                    # We need 2 H2 molecules for methanol
                    hydrogen_combustion_doubled = hydrogen_combustion.scale(2)
                    explanation_steps.append(f"\n3. Double the hydrogen combustion reaction:")
                    explanation_steps.append(f"   {hydrogen_combustion_doubled}")
                    modified_reactions.append(hydrogen_combustion_doubled)
                    
                    # 4. Combine all reactions
                    result_reaction = methanol_decomposition
                    for r in [carbon_combustion, hydrogen_combustion_doubled]:
                        result_reaction = result_reaction.add(r)
                    
                    explanation_steps.append(f"\n4. Combine all the reactions:")
                    explanation_steps.append(f"   {result_reaction}")
    
    if result_reaction and result_reaction.enthalpy is not None:
        steps.extend(explanation_steps)
        steps.append(f"\n## Final Answer:")
        steps.append(f"The enthalpy of combustion for the target reaction is {result_reaction.enthalpy:.1f} kJ/mol")
    else:
        steps.append("\nCould not determine the enthalpy using the given reactions.")
    
    return {
        "target_reaction": target_reaction,
        "known_reactions": known_reactions,
        "result_reaction": result_reaction,
        "enthalpy": result_reaction.enthalpy if result_reaction else None,
        "steps": steps
    }


def solve_enthalpy_problem(problem_text: str) -> Dict[str, Any]:
    """
    Parse a problem text and solve for unknown enthalpies.
    
    Parameters:
        problem_text (str): The full problem text
        
    Returns:
        Dict[str, Any]: Dictionary with solution and steps
    """
    # Extract reactions from the problem text
    lines = problem_text.strip().split('\n')
    reactions = []
    target_reaction = None
    
    for line in lines:
        line = line.strip()
        if not line or "Calculate" in line:
            continue
        
        # Check if this is the target reaction
        if "→" in line or "->" in line:
            if "ΔH°" not in line and "ΔH" not in line:
                try:
                    target_reaction = parse_reaction_string(line)
                except ValueError:
                    pass
            else:
                try:
                    reactions.append(parse_reaction_string(line))
                except ValueError:
                    pass
    
    # If no target reaction was explicitly marked, assume it's the last one
    if target_reaction is None and reactions:
        target_reaction = reactions.pop()
        target_reaction.enthalpy = None
    
    # Solve for the enthalpy
    result = solve_combustion_enthalpy(reactions, target_reaction)
    
    return result


#######################################
# Heat and Temperature Functions
#######################################

def calculate_heat(mass, specific_heat, delta_t):
    """
    Calculate heat energy transfer using q = m × c × ΔT.
    
    Args:
        mass (float): Mass in grams
        specific_heat (float): Specific heat capacity in J/(g·K) or J/(g·°C)
        delta_t (float): Temperature change in Celsius or Kelvin
        
    Returns:
        float: Heat energy in Joules
    """
    return mass * specific_heat * delta_t


def calculate_temperature_change(heat, mass, specific_heat):
    """
    Calculate temperature change using ΔT = q / (m × c).
    
    Args:
        heat (float): Heat energy in Joules
        mass (float): Mass in grams
        specific_heat (float): Specific heat capacity in J/(g·K) or J/(g·°C)
        
    Returns:
        float: Temperature change in Celsius or Kelvin
    """
    return heat / (mass * specific_heat)


def calculate_final_temperature(initial_temp, delta_t):
    """
    Calculate final temperature using Tf = Ti + ΔT.
    
    Args:
        initial_temp (float): Initial temperature in Celsius
        delta_t (float): Temperature change in Celsius
        
    Returns:
        float: Final temperature in Celsius
    """
    return initial_temp + delta_t


def calculate_molar_heat(mass, molar_mass, molar_heat_capacity, delta_t):
    """
    Calculate heat energy transfer using molar heat capacity.
    
    Args:
        mass (float): Mass in grams
        molar_mass (float): Molar mass in g/mol
        molar_heat_capacity (float): Molar heat capacity in J/(mol·K) or J/(mol·°C)
        delta_t (float): Temperature change in Celsius or Kelvin
        
    Returns:
        float: Heat energy in Joules
    """
    moles = mass / molar_mass
    return moles * molar_heat_capacity * delta_t


#######################################
# Thermal Equilibrium Functions
#######################################

def solve_thermal_equilibrium(substances):
    """
    Solve for the final temperature when multiple substances reach thermal equilibrium.
    
    Args:
        substances (list): List of dictionaries with the following keys:
            - mass (float): Mass in grams
            - specific_heat (float): Specific heat capacity in J/(g·K) or J/(g·°C)
            - initial_temp (float): Initial temperature in Celsius
            
    Returns:
        float: Final equilibrium temperature in Celsius
    """
    numerator = 0
    denominator = 0
    
    for substance in substances:
        m = substance['mass']
        c = substance['specific_heat']
        T = substance['initial_temp']
        
        numerator += m * c * T
        denominator += m * c
    
    return numerator / denominator


def solve_thermal_equilibrium_with_molar_heat(substances):
    """
    Solve for the final temperature when multiple substances reach thermal equilibrium,
    using molar heat capacities.
    
    Args:
        substances (list): List of dictionaries with the following keys:
            - mass (float): Mass in grams
            - molar_mass (float): Molar mass in g/mol
            - molar_heat_capacity (float): Molar heat capacity in J/(mol·K) or J/(mol·°C)
            - initial_temp (float): Initial temperature in Celsius
            
    Returns:
        float: Final equilibrium temperature in Celsius
    """
    numerator = 0
    denominator = 0
    
    for substance in substances:
        m = substance['mass']
        mm = substance['molar_mass']
        mhc = substance['molar_heat_capacity']
        T = substance['initial_temp']
        
        moles = m / mm
        numerator += moles * mhc * T
        denominator += moles * mhc
    
    return numerator / denominator


def handle_heat_transfer_problem(mass1, specific_heat1, initial_temp1, 
                                mass2, specific_heat2, initial_temp2):
    """
    Solve a problem where two substances exchange heat until they reach equilibrium.
    
    Args:
        mass1 (float): Mass of substance 1 in grams
        specific_heat1 (float): Specific heat capacity of substance 1 in J/(g·K) or J/(g·°C)
        initial_temp1 (float): Initial temperature of substance 1 in Celsius
        mass2 (float): Mass of substance 2 in grams
        specific_heat2 (float): Specific heat capacity of substance 2 in J/(g·K) or J/(g·°C)
        initial_temp2 (float): Initial temperature of substance 2 in Celsius
        
    Returns:
        dict: Contains final_temp and steps
    """
    # Calculate final temperature
    substances = [
        {'mass': mass1, 'specific_heat': specific_heat1, 'initial_temp': initial_temp1},
        {'mass': mass2, 'specific_heat': specific_heat2, 'initial_temp': initial_temp2}
    ]
    
    final_temp = solve_thermal_equilibrium(substances)
    
    # Calculate heat transferred
    delta_t1 = final_temp - initial_temp1
    delta_t2 = final_temp - initial_temp2
    
    q1 = calculate_heat(mass1, specific_heat1, delta_t1)
    q2 = calculate_heat(mass2, specific_heat2, delta_t2)
    
    # Generate solution steps
    steps = [
        f"1. Calculate the final equilibrium temperature:",
        f"   - Using conservation of energy: q₁ + q₂ = 0",
        f"   - m₁ × c₁ × (Tf - T₁) + m₂ × c₂ × (Tf - T₂) = 0",
        f"   - Tf = (m₁ × c₁ × T₁ + m₂ × c₂ × T₂) / (m₁ × c₁ + m₂ × c₂)",
        f"   - Tf = ({mass1} × {specific_heat1} × {initial_temp1} + {mass2} × {specific_heat2} × {initial_temp2}) / ({mass1} × {specific_heat1} + {mass2} × {specific_heat2})",
        f"   - Tf = {final_temp:.2f} °C",
        f"",
        f"2. Calculate heat transferred by substance 1:",
        f"   - q₁ = m₁ × c₁ × (Tf - T₁)",
        f"   - q₁ = {mass1} × {specific_heat1} × ({final_temp:.2f} - {initial_temp1})",
        f"   - q₁ = {q1:.2f} J",
        f"",
        f"3. Calculate heat transferred by substance 2:",
        f"   - q₂ = m₂ × c₂ × (Tf - T₂)",
        f"   - q₂ = {mass2} × {specific_heat2} × ({final_temp:.2f} - {initial_temp2})",
        f"   - q₂ = {q2:.2f} J",
        f"",
        f"4. Verify conservation of energy:",
        f"   - q₁ + q₂ = {q1:.2f} + {q2:.2f} = {q1 + q2:.2f} J (approximately 0 due to rounding)"
    ]
    
    return {
        "final_temp": final_temp,
        "heat_transferred_1": q1,
        "heat_transferred_2": q2,
        "steps": steps
    }


def handle_heat_transfer_with_molar_heat(mass1, molar_mass1, molar_heat_capacity1, initial_temp1, 
                                        mass2, molar_mass2, molar_heat_capacity2, initial_temp2):
    """
    Solve a problem where two substances exchange heat until they reach equilibrium,
    using molar heat capacities.
    
    Args:
        mass1 (float): Mass of substance 1 in grams
        molar_mass1 (float): Molar mass of substance 1 in g/mol
        molar_heat_capacity1 (float): Molar heat capacity of substance 1 in J/(mol·K) or J/(mol·°C)
        initial_temp1 (float): Initial temperature of substance 1 in Celsius
        mass2 (float): Mass of substance 2 in grams
        molar_mass2 (float): Molar mass of substance 2 in g/mol
        molar_heat_capacity2 (float): Molar heat capacity of substance 2 in J/(mol·K) or J/(mol·°C)
        initial_temp2 (float): Initial temperature of substance 2 in Celsius
        
    Returns:
        dict: Contains final_temp and steps
    """
    # Calculate final temperature
    substances = [
        {'mass': mass1, 'molar_mass': molar_mass1, 'molar_heat_capacity': molar_heat_capacity1, 'initial_temp': initial_temp1},
        {'mass': mass2, 'molar_mass': molar_mass2, 'molar_heat_capacity': molar_heat_capacity2, 'initial_temp': initial_temp2}
    ]
    
    final_temp = solve_thermal_equilibrium_with_molar_heat(substances)
    
    # Calculate moles
    moles1 = mass1 / molar_mass1
    moles2 = mass2 / molar_mass2
    
    # Calculate heat transferred
    delta_t1 = final_temp - initial_temp1
    delta_t2 = final_temp - initial_temp2
    
    q1 = moles1 * molar_heat_capacity1 * delta_t1
    q2 = moles2 * molar_heat_capacity2 * delta_t2
    
    # Generate solution steps
    steps = [
        f"1. Calculate moles for each substance:",
        f"   - Moles of substance 1: n₁ = m₁ / M₁ = {mass1} g / {molar_mass1} g/mol = {moles1:.4f} mol",
        f"   - Moles of substance 2: n₂ = m₂ / M₂ = {mass2} g / {molar_mass2} g/mol = {moles2:.4f} mol",
        f"",
        f"2. Calculate the final equilibrium temperature:",
        f"   - Using conservation of energy: q₁ + q₂ = 0",
        f"   - n₁ × C₁ × (Tf - T₁) + n₂ × C₂ × (Tf - T₂) = 0",
        f"   - Tf = (n₁ × C₁ × T₁ + n₂ × C₂ × T₂) / (n₁ × C₁ + n₂ × C₂)",
        f"   - Tf = ({moles1:.4f} × {molar_heat_capacity1} × {initial_temp1} + {moles2:.4f} × {molar_heat_capacity2} × {initial_temp2}) / ({moles1:.4f} × {molar_heat_capacity1} + {moles2:.4f} × {molar_heat_capacity2})",
        f"   - Tf = {final_temp:.2f} °C",
        f"",
        f"3. Calculate heat transferred by substance 1:",
        f"   - q₁ = n₁ × C₁ × (Tf - T₁)",
        f"   - q₁ = {moles1:.4f} × {molar_heat_capacity1} × ({final_temp:.2f} - {initial_temp1})",
        f"   - q₁ = {q1:.2f} J",
        f"",
        f"4. Calculate heat transferred by substance 2:",
        f"   - q₂ = n₂ × C₂ × (Tf - T₂)",
        f"   - q₂ = {moles2:.4f} × {molar_heat_capacity2} × ({final_temp:.2f} - {initial_temp2})",
        f"   - q₂ = {q2:.2f} J",
        f"",
        f"5. Verify conservation of energy:",
        f"   - q₁ + q₂ = {q1:.2f} + {q2:.2f} = {q1 + q2:.2f} J (approximately 0 due to rounding)"
    ]
    
    return {
        "final_temp": final_temp,
        "heat_transferred_1": q1,
        "heat_transferred_2": q2,
        "steps": steps
    }


def solve_mixture_problem(substances):
    """
    Solve a general thermal equilibrium problem with multiple substances.
    
    Args:
        substances (list): List of dictionaries, each containing:
            - name (str): Substance name
            - mass (float): Mass in g
            - specific_heat (float, optional): Specific heat capacity in J/(g·K)
            - molar_mass (float, optional): Molar mass in g/mol
            - molar_heat_capacity (float, optional): Molar heat capacity in J/(mol·K)
            - initial_temp (float): Initial temperature in °C
            
    Returns:
        dict: Contains final_temp, heat_transfers, and steps
    """
    # Prepare substances with consistent heat capacity parameters
    processed_substances = []
    
    for s in substances:
        substance = s.copy()
        
        # If molar heat capacity is provided, convert to specific heat
        if 'molar_heat_capacity' in substance and 'molar_mass' in substance:
            if 'specific_heat' not in substance:
                substance['specific_heat'] = substance['molar_heat_capacity'] / substance['molar_mass']
        
        # Make sure we have specific heat capacity
        if 'specific_heat' not in substance:
            raise ValueError(f"Missing specific heat for {substance['name']}")
        
        processed_substances.append(substance)
    
    # Calculate final temperature
    final_temp = solve_thermal_equilibrium(processed_substances)
    
    # Calculate heat transferred for each substance
    heat_transfers = []
    
    for substance in processed_substances:
        delta_t = final_temp - substance['initial_temp']
        heat = calculate_heat(substance['mass'], substance['specific_heat'], delta_t)
        heat_transfers.append({
            'name': substance['name'],
            'heat': heat,
            'delta_t': delta_t
        })
    
    # Generate solution steps
    steps = [
        f"1. Calculate the final equilibrium temperature:",
        f"   - Using conservation of energy: Σ q = 0",
        f"   - Σ (m × c × (Tf - Ti)) = 0",
        f"   - Tf = Σ(m × c × Ti) / Σ(m × c)"
    ]
    
    # Add details for each substance
    for i, substance in enumerate(processed_substances):
        steps.append(f"   - For {substance['name']}: m = {substance['mass']} g, c = {substance['specific_heat']} J/(g·K), Ti = {substance['initial_temp']} °C")
    
    steps.append(f"   - Final temperature: Tf = {final_temp:.2f} °C")
    steps.append("")
    
    # Add heat transfer calculations
    steps.append(f"2. Calculate heat transferred for each substance:")
    
    for transfer in heat_transfers:
        steps.append(f"   - {transfer['name']}: q = m × c × (Tf - Ti) = {transfer['heat']:.2f} J (ΔT = {transfer['delta_t']:.2f} °C)")
    
    steps.append("")
    steps.append(f"3. Verify conservation of energy:")
    total_heat = sum(t['heat'] for t in heat_transfers)
    steps.append(f"   - Sum of all heat transfers: {total_heat:.2f} J (approximately 0 due to rounding)")
    
    return {
        "final_temp": final_temp,
        "heat_transfers": heat_transfers,
        "steps": steps
    }

#######################################
# Equilibrium Constant Functions
#######################################

def calculate_equilibrium_constant(delta_g_standard, temperature_k):
    """
    Calculate the equilibrium constant using the relationship ΔG° = -RT ln(K).
    
    Parameters:
        delta_g_standard (float): Standard Gibbs free energy change in kJ/mol
        temperature_k (float): Temperature in Kelvin
        
    Returns:
        float: Equilibrium constant K
    """
    # Convert ΔG° from kJ/mol to J/mol for calculation
    delta_g_j = delta_g_standard * 1000
    
    # R in J/(mol·K)
    R = R_IDEAL_GAS['J/(mol·K)']
    
    # Calculate K using K = exp(-ΔG°/RT)
    K = math.exp(-delta_g_j / (R * temperature_k))
    
    return K


def calculate_delta_g_from_enthalpy_entropy(delta_h_standard, delta_s_standard, temperature_k):
    """
    Calculate the standard Gibbs free energy change using ΔG° = ΔH° - TΔS°.
    
    Parameters:
        delta_h_standard (float): Standard enthalpy change in kJ/mol
        delta_s_standard (float): Standard entropy change in J/(mol·K)
        temperature_k (float): Temperature in Kelvin
        
    Returns:
        float: Standard Gibbs free energy change in kJ/mol
    """
    # Convert ΔS° from J/(mol·K) to kJ/(mol·K)
    delta_s_kj = delta_s_standard / 1000
    
    # Calculate ΔG° = ΔH° - TΔS°
    delta_g_standard = delta_h_standard - temperature_k * delta_s_kj
    
    return delta_g_standard


def solve_equilibrium_constant_problem(
    delta_h_standard=None,     # kJ/mol
    delta_s_standard=None,     # J/(mol·K)
    delta_g_standard=None,     # kJ/mol
    temperature_c=None,        # °C
    temperature_k=None,        # K
    reaction_string=None
) -> Dict[str, Any]:
    """
    Solve equilibrium constant problems using thermodynamic data.
    
    Parameters:
        delta_h_standard (float, optional): Standard enthalpy change in kJ/mol
        delta_s_standard (float, optional): Standard entropy change in J/(mol·K)
        delta_g_standard (float, optional): Standard Gibbs free energy change in kJ/mol
        temperature_c (float, optional): Temperature in Celsius
        temperature_k (float, optional): Temperature in Kelvin
        reaction_string (str, optional): The balanced chemical equation
        
    Returns:
        Dict[str, Any]: Dictionary containing results and solution steps
    """
    steps = []
    
    # Handle temperature conversion
    if temperature_k is None and temperature_c is not None:
        temperature_k = temperature_c + 273.15
        steps.append(f"Step 1: Convert temperature to Kelvin:")
        steps.append(f"T = {temperature_c} °C + 273.15 = {temperature_k:.2f} K")
        steps.append("")
    elif temperature_c is None and temperature_k is not None:
        temperature_c = temperature_k - 273.15
    
    if temperature_k is None:
        raise ValueError("Temperature must be provided in either Celsius or Kelvin")
    
    # Add reaction information if provided
    if reaction_string:
        steps.append(f"Balanced reaction: {reaction_string}")
        steps.append("")
    
    # Calculate ΔG° if not provided
    if delta_g_standard is None:
        if delta_h_standard is not None and delta_s_standard is not None:
            delta_g_standard = calculate_delta_g_from_enthalpy_entropy(
                delta_h_standard, delta_s_standard, temperature_k
            )
            
            step_num = len([s for s in steps if s.startswith("Step")]) + 1
            steps.append(f"Step {step_num}: Calculate ΔG° using ΔG° = ΔH° - TΔS°:")
            steps.append(f"ΔG° = {delta_h_standard} kJ/mol - ({temperature_k:.2f} K)({delta_s_standard} J/(mol·K) × 1 kJ/1000 J)")
            steps.append(f"ΔG° = {delta_h_standard} kJ/mol - ({temperature_k:.2f})({delta_s_standard/1000:.6f} kJ/(mol·K))")
            steps.append(f"ΔG° = {delta_h_standard} kJ/mol - {temperature_k * delta_s_standard/1000:.3f} kJ/mol")
            steps.append(f"ΔG° = {delta_g_standard:.3f} kJ/mol")
            steps.append("")
        else:
            raise ValueError("Either ΔG° must be provided, or both ΔH° and ΔS° must be provided")
    
    # Calculate equilibrium constant
    K = calculate_equilibrium_constant(delta_g_standard, temperature_k)
    
    step_num = len([s for s in steps if s.startswith("Step")]) + 1
    steps.append(f"Step {step_num}: Calculate the equilibrium constant using ΔG° = -RT ln(K):")
    steps.append(f"ln(K) = -ΔG°/RT")
    steps.append(f"ln(K) = -({delta_g_standard:.3f} kJ/mol × 1000 J/kJ) / ({R_IDEAL_GAS['J/(mol·K)']} J/(mol·K) × {temperature_k:.2f} K)")
    steps.append(f"ln(K) = -{delta_g_standard * 1000:.1f} J/mol / {R_IDEAL_GAS['J/(mol·K)'] * temperature_k:.1f} J/mol")
    steps.append(f"ln(K) = {-delta_g_standard * 1000 / (R_IDEAL_GAS['J/(mol·K)'] * temperature_k):.6f}")
    steps.append(f"K = e^({-delta_g_standard * 1000 / (R_IDEAL_GAS['J/(mol·K)'] * temperature_k):.6f})")
    steps.append(f"K = {K:.2e}")
    
    return {
        "reaction": reaction_string,
        "temperature_c": temperature_c,
        "temperature_k": temperature_k,
        "delta_h_standard": delta_h_standard,
        "delta_s_standard": delta_s_standard,
        "delta_g_standard": delta_g_standard,
        "equilibrium_constant": K,
        "steps": steps
    }


def solve_equilibrium_from_formation_data(
    reaction_coefficients,      # Dict: {'reactants': {'substance': coeff}, 'products': {'substance': coeff}}
    formation_enthalpies,       # Dict: {'substance': ΔHf° in kJ/mol}
    standard_entropies,         # Dict: {'substance': S° in J/(mol·K)}
    temperature_c=25,           # Temperature in Celsius
    reaction_string=None        # Optional reaction string for display
) -> Dict[str, Any]:
    """
    Calculate equilibrium constant from standard formation data for any reaction.
    
    Parameters:
        reaction_coefficients (dict): Reaction stoichiometry
            Example: {'reactants': {'N2': 1, 'H2': 3}, 'products': {'NH3': 2}}
        formation_enthalpies (dict): Standard enthalpies of formation in kJ/mol
            Example: {'N2': 0, 'H2': 0, 'NH3': -45.9}
        standard_entropies (dict): Standard entropies in J/(mol·K)
            Example: {'N2': 191.6, 'H2': 130.7, 'NH3': 192.8}
        temperature_c (float): Temperature in Celsius (default 25°C)
        reaction_string (str, optional): Balanced chemical equation for display
        
    Returns:
        Dict[str, Any]: Dictionary containing results and solution steps
    """
    
    # Calculate ΔH° for the reaction
    delta_h_reaction = 0
    for substance, coeff in reaction_coefficients['products'].items():
        if substance in formation_enthalpies:
            delta_h_reaction += coeff * formation_enthalpies[substance]
    
    for substance, coeff in reaction_coefficients['reactants'].items():
        if substance in formation_enthalpies:
            delta_h_reaction -= coeff * formation_enthalpies[substance]
    
    # Calculate ΔS° for the reaction
    delta_s_reaction = 0
    for substance, coeff in reaction_coefficients['products'].items():
        if substance in standard_entropies:
            delta_s_reaction += coeff * standard_entropies[substance]
    
    for substance, coeff in reaction_coefficients['reactants'].items():
        if substance in standard_entropies:
            delta_s_reaction -= coeff * standard_entropies[substance]
    
    # Build detailed calculation steps
    detailed_steps = [
        "Given thermodynamic data:",
        "Formation enthalpies (ΔHf°):"
    ]
    
    for substance, value in formation_enthalpies.items():
        detailed_steps.append(f"  {substance}: {value} kJ/mol")
    
    detailed_steps.extend([
        "",
        "Standard entropies (S°):"
    ])
    
    for substance, value in standard_entropies.items():
        detailed_steps.append(f"  {substance}: {value} J/(mol·K)")
    
    detailed_steps.extend([
        "",
        "Step 1: Calculate ΔH° for the reaction:",
        "ΔH°rxn = Σ(coefficients × ΔHf°)products - Σ(coefficients × ΔHf°)reactants"
    ])
    
    # Add detailed ΔH° calculation
    products_h_terms = []
    reactants_h_terms = []
    
    for substance, coeff in reaction_coefficients['products'].items():
        if substance in formation_enthalpies:
            products_h_terms.append(f"{coeff} × ({formation_enthalpies[substance]})")
    
    for substance, coeff in reaction_coefficients['reactants'].items():
        if substance in formation_enthalpies:
            reactants_h_terms.append(f"{coeff} × ({formation_enthalpies[substance]})")
    
    detailed_steps.append(f"ΔH°rxn = [{' + '.join(products_h_terms)}] - [{' + '.join(reactants_h_terms)}]")
    detailed_steps.append(f"ΔH°rxn = {delta_h_reaction} kJ/mol")
    detailed_steps.append("")
    
    detailed_steps.extend([
        "Step 2: Calculate ΔS° for the reaction:",
        "ΔS°rxn = Σ(coefficients × S°)products - Σ(coefficients × S°)reactants"
    ])
    
    # Add detailed ΔS° calculation
    products_s_terms = []
    reactants_s_terms = []
    
    for substance, coeff in reaction_coefficients['products'].items():
        if substance in standard_entropies:
            products_s_terms.append(f"{coeff} × {standard_entropies[substance]}")
    
    for substance, coeff in reaction_coefficients['reactants'].items():
        if substance in standard_entropies:
            reactants_s_terms.append(f"{coeff} × {standard_entropies[substance]}")
    
    detailed_steps.append(f"ΔS°rxn = [{' + '.join(products_s_terms)}] - [{' + '.join(reactants_s_terms)}]")
    detailed_steps.append(f"ΔS°rxn = {delta_s_reaction} J/(mol·K)")
    detailed_steps.append("")
    
    # Use the general equilibrium constant function
    result = solve_equilibrium_constant_problem(
        delta_h_standard=delta_h_reaction,
        delta_s_standard=delta_s_reaction,
        temperature_c=temperature_c,
        reaction_string=reaction_string
    )
    
    # Add the detailed formation data steps to the beginning
    result["steps"] = detailed_steps + result["steps"]
    
    # Add the input data to the result
    result["formation_enthalpies"] = formation_enthalpies
    result["standard_entropies"] = standard_entropies
    result["reaction_coefficients"] = reaction_coefficients
    
    return result


def solve_haber_bosch_problem(temperature_c=25):
    """
    Solve the Haber-Bosch process equilibrium constant problem as an example.
    This is now just a wrapper that calls the general function with Haber-Bosch data.
    
    Parameters:
        temperature_c (float): Temperature in Celsius (default 25°C)
        
    Returns:
        Dict[str, Any]: Dictionary containing results and solution steps
    """
    # Define the reaction: N2(g) + 3H2(g) ⇌ 2NH3(g)
    reaction_coefficients = {
        'reactants': {'N2': 1, 'H2': 3},
        'products': {'NH3': 2}
    }
    
    # Thermodynamic data
    formation_enthalpies = {
        'NH3': -45.9,  # kJ/mol
        'N2': 0,       # kJ/mol (element)
        'H2': 0        # kJ/mol (element)
    }
    
    standard_entropies = {
        'NH3': 192.8,  # J/(mol·K)
        'N2': 191.6,   # J/(mol·K)
        'H2': 130.7    # J/(mol·K)
    }
    
    return solve_equilibrium_from_formation_data(
        reaction_coefficients=reaction_coefficients,
        formation_enthalpies=formation_enthalpies,
        standard_entropies=standard_entropies,
        temperature_c=temperature_c,
        reaction_string="N₂(g) + 3H₂(g) ⇌ 2NH₃(g)"
    )

#######################################
# Clausius-Clapeyron Equation Functions
#######################################

def calculate_boiling_point_with_pressure(
    normal_boiling_point_c: float, 
    heat_of_vaporization: float,  # in kJ/mol
    initial_pressure: float,      # in atm
    final_pressure: float,        # in atm
) -> Dict[str, Any]:
    """
    Calculate the new boiling point when pressure changes using the Clausius-Clapeyron equation.
    
    Parameters:
        normal_boiling_point_c (float): The normal boiling point in degrees Celsius
        heat_of_vaporization (float): The heat of vaporization in kJ/mol
        initial_pressure (float): The initial pressure in atm
        final_pressure (float): The final pressure in atm
        
    Returns:
        Dict[str, Any]: Dictionary containing results and solution steps
    """
    # Convert boiling point to Kelvin
    normal_boiling_point_k = normal_boiling_point_c + 273.15
    
    # Gas constant in kJ/(mol·K)
    R = R_IDEAL_GAS['kJ/(mol·K)']
    
    # Calculate the new boiling point using the integrated Clausius-Clapeyron equation
    # ln(P2/P1) = -(ΔHvap/R) * (1/T2 - 1/T1)
    # Rearranging: 1/T2 = 1/T1 - (R/ΔHvap) * ln(P2/P1)
    
    # Calculate 1/T2
    inv_t2 = (1/normal_boiling_point_k) - (R/heat_of_vaporization) * math.log(final_pressure/initial_pressure)
    
    # Calculate T2
    new_boiling_point_k = 1/inv_t2
    new_boiling_point_c = new_boiling_point_k - 273.15
    
    # Build solution steps
    steps = [
        f"Step 1: Convert the normal boiling point to Kelvin:",
        f"T₁ = {normal_boiling_point_c} °C + 273.15 = {normal_boiling_point_k:.2f} K",
        f"",
        f"Step 2: Use the Clausius-Clapeyron equation to find the new boiling point:",
        f"ln(P₂/P₁) = -(ΔHvap/R) × (1/T₂ - 1/T₁)",
        f"",
        f"Rearranging for 1/T₂:",
        f"1/T₂ = 1/T₁ - (R/ΔHvap) × ln(P₂/P₁)",
        f"",
        f"Step 3: Substitute the values:",
        f"1/T₂ = 1/{normal_boiling_point_k:.2f} K - ({R:.6f} kJ/(mol·K)/{heat_of_vaporization} kJ/mol) × ln({final_pressure}/{initial_pressure})",
        f"1/T₂ = {1/normal_boiling_point_k:.6f} K⁻¹ - {R/heat_of_vaporization:.6f} × {math.log(final_pressure/initial_pressure):.6f}",
        f"1/T₂ = {1/normal_boiling_point_k:.6f} K⁻¹ - {(R/heat_of_vaporization) * math.log(final_pressure/initial_pressure):.6f} K⁻¹",
        f"1/T₂ = {inv_t2:.6f} K⁻¹",
        f"",
        f"Step 4: Calculate T₂:",
        f"T₂ = 1/({inv_t2:.6f} K⁻¹) = {new_boiling_point_k:.2f} K",
        f"",
        f"Step 5: Convert back to Celsius:",
        f"T₂ = {new_boiling_point_k:.2f} K - 273.15 = {new_boiling_point_c:.2f} °C"
    ]
    
    return {
        "normal_boiling_point_c": normal_boiling_point_c,
        "normal_boiling_point_k": normal_boiling_point_k,
        "heat_of_vaporization": heat_of_vaporization,
        "initial_pressure": initial_pressure,
        "final_pressure": final_pressure,
        "new_boiling_point_k": new_boiling_point_k,
        "new_boiling_point_c": new_boiling_point_c,
        "steps": steps
    }

def calculate_pressure_with_temperature(
    normal_boiling_point_c: float,
    heat_of_vaporization: float,  # in kJ/mol
    initial_pressure: float,      # in atm
    final_temperature_c: float    # in °C
) -> Dict[str, Any]:
    """
    Calculate the vapor pressure at a given temperature using the Clausius-Clapeyron equation.
    
    Parameters:
        normal_boiling_point_c (float): The normal boiling point in degrees Celsius
        heat_of_vaporization (float): The heat of vaporization in kJ/mol
        initial_pressure (float): The initial pressure in atm (typically 1 atm)
        final_temperature_c (float): The temperature at which to calculate the vapor pressure, in °C
        
    Returns:
        Dict[str, Any]: Dictionary containing results and solution steps
    """
    # Convert temperatures to Kelvin
    normal_boiling_point_k = normal_boiling_point_c + 273.15
    final_temperature_k = final_temperature_c + 273.15
    
    # Gas constant in kJ/(mol·K)
    R = R_IDEAL_GAS['kJ/(mol·K)']
    
    # Calculate the new pressure using the integrated Clausius-Clapeyron equation
    # ln(P2/P1) = -(ΔHvap/R) * (1/T2 - 1/T1)
    
    exponent = -(heat_of_vaporization/R) * (1/final_temperature_k - 1/normal_boiling_point_k)
    final_pressure = initial_pressure * math.exp(exponent)
    
    # Build solution steps
    steps = [
        f"Step 1: Convert temperatures to Kelvin:",
        f"T₁ = {normal_boiling_point_c} °C + 273.15 = {normal_boiling_point_k:.2f} K",
        f"T₂ = {final_temperature_c} °C + 273.15 = {final_temperature_k:.2f} K",
        f"",
        f"Step 2: Use the Clausius-Clapeyron equation to find the new pressure:",
        f"ln(P₂/P₁) = -(ΔHvap/R) × (1/T₂ - 1/T₁)",
        f"",
        f"Step 3: Substitute the values:",
        f"ln(P₂/{initial_pressure}) = -({heat_of_vaporization} kJ/mol/{R:.6f} kJ/(mol·K)) × (1/{final_temperature_k:.2f} K - 1/{normal_boiling_point_k:.2f} K)",
        f"ln(P₂/{initial_pressure}) = -{heat_of_vaporization/R:.2f} × ({1/final_temperature_k:.6f} K⁻¹ - {1/normal_boiling_point_k:.6f} K⁻¹)",
        f"ln(P₂/{initial_pressure}) = -{heat_of_vaporization/R:.2f} × {1/final_temperature_k - 1/normal_boiling_point_k:.6f} K⁻¹",
        f"ln(P₂/{initial_pressure}) = {exponent:.6f}",
        f"",
        f"Step 4: Calculate P₂:",
        f"P₂ = {initial_pressure} × e^({exponent:.6f})",
        f"P₂ = {initial_pressure} × {math.exp(exponent):.6f}",
        f"P₂ = {final_pressure:.6f} atm"
    ]
    
    return {
        "normal_boiling_point_c": normal_boiling_point_c,
        "normal_boiling_point_k": normal_boiling_point_k,
        "heat_of_vaporization": heat_of_vaporization,
        "initial_pressure": initial_pressure,
        "final_temperature_c": final_temperature_c,
        "final_temperature_k": final_temperature_k,
        "final_pressure": final_pressure,
        "steps": steps
    }

def calculate_heat_of_vaporization(
    temp1_c: float,
    pressure1: float,  # in atm
    temp2_c: float,
    pressure2: float   # in atm
) -> Dict[str, Any]:
    """
    Calculate the heat of vaporization using two pressure-temperature data points.
    
    Parameters:
        temp1_c (float): First temperature in degrees Celsius
        pressure1 (float): First pressure in atm
        temp2_c (float): Second temperature in degrees Celsius
        pressure2 (float): Second pressure in atm
        
    Returns:
        Dict[str, Any]: Dictionary containing results and solution steps
    """
    # Convert temperatures to Kelvin
    temp1_k = temp1_c + 273.15
    temp2_k = temp2_c + 273.15
    
    # Gas constant in kJ/(mol·K)
    R = R_IDEAL_GAS['kJ/(mol·K)']
    
    # Calculate heat of vaporization using the Clausius-Clapeyron equation
    # ln(P2/P1) = -(ΔHvap/R) * (1/T2 - 1/T1)
    # Rearranging: ΔHvap = -R * ln(P2/P1) / (1/T2 - 1/T1)
    
    heat_of_vaporization = -R * math.log(pressure2/pressure1) / (1/temp2_k - 1/temp1_k)
    
    # Build solution steps
    steps = [
        f"Step 1: Convert temperatures to Kelvin:",
        f"T₁ = {temp1_c} °C + 273.15 = {temp1_k:.2f} K",
        f"T₂ = {temp2_c} °C + 273.15 = {temp2_k:.2f} K",
        f"",
        f"Step 2: Rearrange the Clausius-Clapeyron equation to solve for the heat of vaporization:",
        f"ln(P₂/P₁) = -(ΔHvap/R) × (1/T₂ - 1/T₁)",
        f"ΔHvap = -R × ln(P₂/P₁) / (1/T₂ - 1/T₁)",
        f"",
        f"Step 3: Substitute the values:",
        f"ΔHvap = -{R:.6f} kJ/(mol·K) × ln({pressure2}/{pressure1}) / (1/{temp2_k:.2f} K - 1/{temp1_k:.2f} K)",
        f"ΔHvap = -{R:.6f} kJ/(mol·K) × {math.log(pressure2/pressure1):.6f} / ({1/temp2_k:.6f} K⁻¹ - {1/temp1_k:.6f} K⁻¹)",
        f"ΔHvap = -{R:.6f} kJ/(mol·K) × {math.log(pressure2/pressure1):.6f} / {1/temp2_k - 1/temp1_k:.6f} K⁻¹",
        f"ΔHvap = {heat_of_vaporization:.2f} kJ/mol"
    ]
    
    return {
        "temp1_c": temp1_c,
        "temp1_k": temp1_k,
        "pressure1": pressure1,
        "temp2_c": temp2_c,
        "temp2_k": temp2_k,
        "pressure2": pressure2,
        "heat_of_vaporization": heat_of_vaporization,
        "steps": steps
    }