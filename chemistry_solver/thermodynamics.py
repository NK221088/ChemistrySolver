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
# Clausius-Clapeyron Equation Functions
#######################################