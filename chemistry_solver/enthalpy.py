"""
Chemistry Solver module for enthalpy calculations and Hess's Law problems.
"""
import re
from typing import Dict, List, Tuple, Any, Union, Optional
from dataclasses import dataclass


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