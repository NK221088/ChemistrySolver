"""
Refactored Redox Reactions Module for Chemistry Problem Solver

This module provides comprehensive functionality for analyzing and working with redox reactions,
including parsing, balancing, favorability determination, and oxidation state analysis.
"""

import math
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass

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


@dataclass
class RedoxElement:
    """Represents an element undergoing oxidation or reduction."""
    element: str
    from_state: int
    to_state: int
    change: int
    reactant_compound: str
    product_compound: str
    
    @property
    def is_oxidized(self) -> bool:
        """Returns True if element is oxidized (loses electrons)."""
        return self.from_state < self.to_state
    
    @property
    def is_reduced(self) -> bool:
        """Returns True if element is reduced (gains electrons)."""
        return self.from_state > self.to_state


@dataclass
class RedoxAnalysis:
    """Complete analysis of a redox reaction."""
    oxidized_elements: List[RedoxElement]
    reduced_elements: List[RedoxElement]
    oxidation_analysis: List[str]
    is_redox: bool
    
    @property
    def total_electrons_lost(self) -> int:
        """Total electrons lost by all oxidized elements."""
        return sum(elem.change for elem in self.oxidized_elements)
    
    @property
    def total_electrons_gained(self) -> int:
        """Total electrons gained by all reduced elements."""
        return sum(elem.change for elem in self.reduced_elements)
    
    @property
    def is_balanced(self) -> bool:
        """Returns True if electrons lost equals electrons gained."""
        return self.total_electrons_lost == self.total_electrons_gained


@dataclass
class FavorabilityResult:
    """Result of redox favorability analysis."""
    favorable: Optional[bool]
    message: str
    steps: List[str]
    oxidation_half: Optional[str] = None
    reduction_half: Optional[str] = None
    e_cell: Optional[float] = None
    oxidizing_agent: Optional[str] = None
    reducing_agent: Optional[str] = None
    balanced_equation: Optional[str] = None
    notes: Optional[str] = None


class ReactionParser:
    """Handles parsing and standardization of chemical reactions."""
    
    @staticmethod
    def standardize_reaction(reaction: str) -> str:
        """Standardize reaction notation."""
        reaction = reaction.replace("→", "->").replace("⟶", "->")
        reaction = reaction.replace("^2-", "²⁻").replace("^2+", "²⁺")
        reaction = reaction.replace("^3-", "³⁻").replace("^3+", "³⁺")
        reaction = reaction.replace(" - ", " + ")
        return reaction.strip()
    
    @staticmethod
    def split_reaction(reaction: str) -> Tuple[List[str], List[str]]:
        """Split reaction into reactants and products."""
        standardized = ReactionParser.standardize_reaction(reaction)
        parts = standardized.split("->")
        
        if len(parts) != 2:
            raise ValueError("Invalid reaction format. Use '->' or '→' to separate reactants and products.")
        
        reactants = [r.strip() for r in parts[0].split("+")]
        products = [p.strip() for p in parts[1].split("+")]
        
        return reactants, products


class ComplexReactionHandler:
    """Handles complex redox reactions with predefined patterns."""
    
    COMPLEX_PATTERNS = {
        "permanganate_so2": {
            "conditions": lambda r: "MnO4-" in r and "Mn2+" in r and "H+" in r and "SO2" in r and "SO4" in r,
            "data": {
                "oxidation_half": "SO2 + 2H2O → SO4²⁻ + 4H⁺ + 2e⁻",
                "reduction_half": "MnO4⁻ + 8H⁺ + 5e⁻ → Mn²⁺ + 4H2O",
                "balanced_equation": "5SO2 + 2MnO4⁻ + 2H2O → 5SO4²⁻ + 2Mn²⁺ + 4H⁺",
                "oxidizing_agent": "MnO4⁻ (Permanganate ion)",
                "reducing_agent": "SO2 (Sulfur dioxide)",
                "electron_transfer": "10 electrons (5 × 2 = 10)",
                "notes": "Sulfur is oxidized from +4 to +6, while manganese is reduced from +7 to +2."
            }
        },
        "dichromate_fe": {
            "conditions": lambda r: "Cr2O7" in r and "Cr3+" in r and "Fe2+" in r and "Fe3+" in r,
            "data": {
                "oxidation_half": "Fe²⁺ → Fe³⁺ + e⁻",
                "reduction_half": "Cr2O7²⁻ + 14H⁺ + 6e⁻ → 2Cr³⁺ + 7H2O",
                "balanced_equation": "Cr2O7²⁻ + 6Fe²⁺ + 14H⁺ → 2Cr³⁺ + 6Fe³⁺ + 7H2O",
                "oxidizing_agent": "Cr2O7²⁻ (Dichromate ion)",
                "reducing_agent": "Fe²⁺ (Iron(II) ion)",
                "electron_transfer": "6 electrons",
                "notes": "Common titration reaction. Iron is oxidized from +2 to +3, while chromium is reduced from +6 to +3."
            }
        },
        "dichromate_so2": {
            "conditions": lambda r: "Cr2O7" in r and ("SO2" in r or "HSO3" in r or "SO3" in r),
            "data": {
                "oxidation_half": "SO2 + 2H2O → SO4²⁻ + 4H⁺ + 2e⁻",
                "reduction_half": "Cr2O7²⁻ + 14H⁺ + 6e⁻ → 2Cr³⁺ + 7H2O",
                "balanced_equation": "Cr2O7²⁻ + 3SO2 + 2H⁺ → 2Cr³⁺ + 3SO4²⁻ + H2O",
                "oxidizing_agent": "Cr2O7²⁻ (Dichromate ion)",
                "reducing_agent": "SO2 (Sulfur dioxide)",
                "electron_transfer": "6 electrons",
                "notes": "Sulfur is oxidized from +4 to +6, while chromium is reduced from +6 to +3."
            }
        },
        "h2o2_disproportionation": {
            "conditions": lambda r: "H2O2" in r and "H2O" in r and "O2" in r,
            "data": {
                "oxidation_half": "H2O2 → O2 + 2H⁺ + 2e⁻",
                "reduction_half": "H2O2 + 2H⁺ + 2e⁻ → 2H2O",
                "balanced_equation": "2H2O2 → 2H2O + O2",
                "oxidizing_agent": "H2O2 (Hydrogen peroxide)",
                "reducing_agent": "H2O2 (Hydrogen peroxide)",
                "electron_transfer": "2 electrons",
                "notes": "Disproportionation reaction where hydrogen peroxide acts as both oxidizing and reducing agent. Oxygen changes from -1 to 0 (oxidation) and from -1 to -2 (reduction)."
            }
        }
    }
    
    @classmethod
    def parse_complex_reaction(cls, reaction: str) -> Dict:
        """Parse complex redox reactions using predefined patterns."""
        for pattern_name, pattern_info in cls.COMPLEX_PATTERNS.items():
            if pattern_info["conditions"](reaction):
                return pattern_info["data"]
        
        # Return template for unknown complex reactions
        return {
            "oxidation_half": "[Please complete the oxidation half-reaction]",
            "reduction_half": "[Please complete the reduction half-reaction]",
            "balanced_equation": "[Please provide the balanced equation]",
            "oxidizing_agent": "[Please identify the oxidizing agent]",
            "reducing_agent": "[Please identify the reducing agent]",
            "notes": "This complex redox reaction needs manual balancing using the half-reaction method."
        }


class OxidationStateAnalyzer:
    """Analyzes oxidation state changes in chemical reactions."""
    
    @staticmethod
    def clean_compound(compound: str) -> str:
        """Remove state notations from compound formula."""
        return compound.replace("(s)", "").replace("(aq)", "").replace("(g)", "").replace("(l)", "").strip()
    
    @classmethod
    def analyze_reaction(cls, reaction: str) -> RedoxAnalysis:
        """Analyze oxidation state changes in a reaction."""
        reactants, products = ReactionParser.split_reaction(reaction)
        
        oxidation_changes = {}
        element_compounds = {}
        
        # Analyze reactants
        cls._analyze_compounds(reactants, oxidation_changes, element_compounds, "reactant")
        
        # Analyze products
        cls._analyze_compounds(products, oxidation_changes, element_compounds, "product")
        
        # Identify oxidized and reduced elements
        oxidized_elements = []
        reduced_elements = []
        oxidation_analysis = []
        
        for element, states in oxidation_changes.items():
            if states["reactant"] and states["product"]:
                reactant_state = states["reactant"][0]
                product_state = states["product"][0]
                
                if reactant_state < product_state:
                    # Oxidation
                    redox_elem = RedoxElement(
                        element=element,
                        from_state=reactant_state,
                        to_state=product_state,
                        change=product_state - reactant_state,
                        reactant_compound=element_compounds[element]["reactant"][0],
                        product_compound=element_compounds[element]["product"][0]
                    )
                    oxidized_elements.append(redox_elem)
                    oxidation_analysis.append(
                        f"{element}: {reactant_state} → {product_state} (oxidized, loses {product_state - reactant_state} electrons)"
                    )
                    
                elif reactant_state > product_state:
                    # Reduction
                    redox_elem = RedoxElement(
                        element=element,
                        from_state=reactant_state,
                        to_state=product_state,
                        change=reactant_state - product_state,
                        reactant_compound=element_compounds[element]["reactant"][0],
                        product_compound=element_compounds[element]["product"][0]
                    )
                    reduced_elements.append(redox_elem)
                    oxidation_analysis.append(
                        f"{element}: {reactant_state} → {product_state} (reduced, gains {reactant_state - product_state} electrons)"
                    )
                else:
                    oxidation_analysis.append(f"{element}: {reactant_state} → {product_state} (no change)")
        
        return RedoxAnalysis(
            oxidized_elements=oxidized_elements,
            reduced_elements=reduced_elements,
            oxidation_analysis=oxidation_analysis,
            is_redox=len(oxidized_elements) > 0 and len(reduced_elements) > 0
        )
    
    @classmethod
    def _analyze_compounds(cls, compounds: List[str], oxidation_changes: Dict, element_compounds: Dict, side: str):
        """Analyze compounds for oxidation states."""
        for compound in compounds:
            clean_compound = cls.clean_compound(compound)
            
            try:
                elements, charge = parse_formula(clean_compound)
                for element in elements:
                    try:
                        ox_result = calculate_oxidation_number(clean_compound, element)
                        ox_state = ox_result["oxidation_number"]
                        
                        if element not in oxidation_changes:
                            oxidation_changes[element] = {"reactant": [], "product": []}
                            element_compounds[element] = {"reactant": [], "product": []}
                        
                        oxidation_changes[element][side].append(ox_state)
                        element_compounds[element][side].append(compound)
                    except (ValueError, KeyError):
                        continue
            except Exception:
                continue


class PotentialCalculator:
    """Handles standard reduction potential calculations."""
    
    @staticmethod
    def find_standard_reduction_potential(half_reaction: str) -> Optional[float]:
        """Find standard reduction potential for a half-reaction."""
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
    
    @staticmethod
    def calculate_nernst_equation(e_cell_standard: float, n: int, concentration_ratio: float, temperature: float = 298.15) -> float:
        """Calculate cell potential using the Nernst equation."""
        R = 8.314  # J/(mol·K)
        F = 96485  # C/mol
        
        e_cell = e_cell_standard - ((R * temperature) / (n * F)) * 2.303 * math.log10(concentration_ratio)
        return e_cell


class ActivitySeriesAnalyzer:
    """Analyzes reactions using the activity series."""
    
    @staticmethod
    def can_displace(reducing_metal: str, oxidizing_ion: str) -> bool:
        """Check if reducing metal can displace oxidizing ion based on activity series."""
        try:
            red_position = ACTIVITY_SERIES.index(reducing_metal)
            oxd_position = ACTIVITY_SERIES.index(oxidizing_ion)
            return red_position < oxd_position
        except ValueError:
            return None
    
    @staticmethod
    def get_activity_positions(reducing_metal: str, oxidizing_ion: str) -> Tuple[Optional[int], Optional[int]]:
        """Get positions of metals in activity series."""
        try:
            red_position = ACTIVITY_SERIES.index(reducing_metal)
        except ValueError:
            red_position = None
        
        try:
            oxd_position = ACTIVITY_SERIES.index(oxidizing_ion)
        except ValueError:
            oxd_position = None
        
        return red_position, oxd_position


class RedoxParser:
    """Main class for parsing redox reactions."""
    
    def __init__(self):
        self.complex_handler = ComplexReactionHandler()
        self.ox_analyzer = OxidationStateAnalyzer()
        self.potential_calc = PotentialCalculator()
        self.activity_analyzer = ActivitySeriesAnalyzer()
    
    def parse_reaction(self, reaction: str) -> Dict:
        """Parse a redox reaction and identify components."""
        standardized_reaction = ReactionParser.standardize_reaction(reaction)
        
        # First check if this is a common redox reaction
        common_result = self._check_common_reactions(standardized_reaction)
        if common_result:
            return common_result
        
        # Check for complex reactions
        complex_result = self.complex_handler.parse_complex_reaction(standardized_reaction)
        if not complex_result["oxidation_half"].startswith("[Please"):
            return {**complex_result, "is_complex_reaction": True}
        
        # Try standard parsing
        return self._parse_standard_reaction(standardized_reaction)
    
    def _check_common_reactions(self, reaction: str) -> Optional[Dict]:
        """Check if reaction matches common redox reactions."""
        for name, data in COMMON_REDOX_REACTIONS.items():
            if reaction in data["reaction"] or reaction in data["balanced_equation"]:
                return {
                    "oxidizing_agent": self._get_oxidizing_agent_from_reaction(data),
                    "reducing_agent": self._get_reducing_agent_from_reaction(data),
                    "oxidation_half": data["oxidation_half"],
                    "reduction_half": data["reduction_half"],
                    "electron_transfer": self._get_electron_transfer(data["oxidation_half"]),
                    "is_common_reaction": True,
                    "common_name": name
                }
        
        # Check common redox pairs
        for name, data in COMMON_REDOX_PAIRS.items():
            if reaction in data["net_reaction"]:
                return {
                    "oxidizing_agent": self._get_oxidizing_agent_from_pair(data),
                    "reducing_agent": self._get_reducing_agent_from_pair(data),
                    "oxidation_half": data["oxidation"],
                    "reduction_half": data["reduction"],
                    "electron_transfer": self._get_electron_transfer(data["oxidation"]),
                    "is_common_pair": True,
                    "common_name": name,
                    "e_cell": data["e_cell"],
                    "favorable": data["favorable"]
                }
        
        return None
    
    def _parse_standard_reaction(self, reaction: str) -> Dict:
        """Parse standard metal displacement reactions."""
        reactants, products = ReactionParser.split_reaction(reaction)
        
        # Look for metals and metal ions
        metal_reactants = [r for r in reactants if any(metal in r and "(s)" in r for metal in COMMON_METALS)]
        metal_products = [p for p in products if any(metal in p and "(s)" in p for metal in COMMON_METALS)]
        ion_reactants = [r for r in reactants if any(metal in r and any(ion_state in r for ion_state in ["(aq)", "+"]) for metal in COMMON_METALS)]
        ion_products = [p for p in products if any(metal in p and any(ion_state in p for ion_state in ["(aq)", "+"]) for metal in COMMON_METALS)]
        
        # Common case: metal displacement reaction
        if len(metal_reactants) == 1 and len(ion_reactants) == 1 and len(metal_products) == 1 and len(ion_products) == 1:
            reducing_agent = metal_reactants[0]
            oxidizing_agent = ion_reactants[0]
            
            # Extract metal symbols
            reducing_metal = ''.join([c for c in reducing_agent if c.isalpha()])
            oxidizing_ion = ''.join([c for c in oxidizing_agent if c.isalpha()])
            
            # Create half-reactions (assuming 2+ ions for simplicity)
            oxidation_half = f"{reducing_metal}(s) -> {reducing_metal}²⁺(aq) + 2e⁻"
            reduction_half = f"{oxidizing_ion}²⁺(aq) + 2e⁻ -> {oxidizing_ion}(s)"
            
            return {
                "oxidizing_agent": oxidizing_agent,
                "reducing_agent": reducing_agent,
                "oxidation_half": oxidation_half,
                "reduction_half": reduction_half,
                "electron_transfer": 2
            }
        
        # For complex cases, return placeholder
        return {
            "oxidizing_agent": None,
            "reducing_agent": None,
            "oxidation_half": "[Oxidation half-reaction could not be automatically determined]",
            "reduction_half": "[Reduction half-reaction could not be automatically determined]",
            "electron_transfer": None
        }
    
    def _get_oxidizing_agent_from_reaction(self, reaction_data: Dict) -> str:
        """Extract oxidizing agent from reaction data."""
        reduction_half = reaction_data["reduction_half"]
        parts = reduction_half.split("->")[0].strip().split("+")
        return parts[0].strip()
    
    def _get_reducing_agent_from_reaction(self, reaction_data: Dict) -> str:
        """Extract reducing agent from reaction data."""
        oxidation_half = reaction_data["oxidation_half"]
        parts = oxidation_half.split("->")[0].strip().split("+")
        return parts[0].strip()
    
    def _get_oxidizing_agent_from_pair(self, pair_data: Dict) -> str:
        """Extract oxidizing agent from redox pair data."""
        reduction_half = pair_data["reduction"]
        parts = reduction_half.split("->")[0].strip().split("+")
        return parts[0].strip()
    
    def _get_reducing_agent_from_pair(self, pair_data: Dict) -> str:
        """Extract reducing agent from redox pair data."""
        oxidation_half = pair_data["oxidation"]
        parts = oxidation_half.split("->")[0].strip().split("+")
        return parts[0].strip()
    
    def _get_electron_transfer(self, half_reaction: str) -> int:
        """Extract number of electrons from half-reaction."""
        if "e-" in half_reaction or "e⁻" in half_reaction:
            half_reaction = half_reaction.replace("e⁻", "e-")
            parts = half_reaction.split("e-")
            left_part = parts[0].strip()
            electron_part = left_part.split()[-1].strip()
            
            if electron_part.isdigit():
                return int(electron_part)
            elif electron_part == "":
                return 1
        
        return 2  # Default assumption


class FavorabilityAnalyzer:
    """Analyzes redox reaction favorability."""
    
    def __init__(self):
        self.parser = RedoxParser()
        self.potential_calc = PotentialCalculator()
        self.activity_analyzer = ActivitySeriesAnalyzer()
        self.ox_analyzer = OxidationStateAnalyzer()
    
    def determine_favorability(self, reaction: str) -> FavorabilityResult:
        """Determine if a redox reaction is favorable."""
        try:
            # First validate it's a redox reaction
            validation = self._validate_redox_reaction(reaction)
            if not validation["is_redox"]:
                return FavorabilityResult(
                    favorable=None,
                    message="Not a redox reaction",
                    steps=validation["validation_steps"] + ["Cannot determine favorability for non-redox reactions."]
                )
            
            # Check common redox pairs first
            common_result = self._check_common_pairs(reaction)
            if common_result:
                return common_result
            
            # Parse the reaction
            parsed = self.parser.parse_reaction(reaction)
            
            # Try different analysis methods
            if parsed.get("is_common_reaction"):
                return self._analyze_common_reaction(parsed)
            
            if self._is_metal_displacement(parsed):
                return self._analyze_metal_displacement(parsed, reaction)
            
            # Try potential calculation
            return self._analyze_by_potentials(parsed)
            
        except Exception as e:
            return FavorabilityResult(
                favorable=None,
                message=f"Error: {str(e)}",
                steps=[f"An error occurred: {str(e)}"]
            )
    
    def _validate_redox_reaction(self, reaction: str) -> Dict:
        """Validate if reaction is actually redox."""
        try:
            analysis = self.ox_analyzer.analyze_reaction(reaction)
            
            validation_steps = [
                f"1. Analyzing oxidation states in reaction: {reaction}",
                "2. Oxidation state changes found:"
            ]
            
            if analysis.oxidation_analysis:
                validation_steps.extend([f"   - {item}" for item in analysis.oxidation_analysis])
            else:
                validation_steps.append("   - No oxidation state changes detected")
            
            if analysis.is_redox:
                validation_steps.extend([
                    "3. This IS a redox reaction:",
                    f"   - Elements oxidized: {', '.join([elem.element for elem in analysis.oxidized_elements])}",
                    f"   - Elements reduced: {', '.join([elem.element for elem in analysis.reduced_elements])}"
                ])
            else:
                validation_steps.append("3. This is NOT a redox reaction (no electron transfer detected)")
            
            return {
                "is_redox": analysis.is_redox,
                "validation_steps": validation_steps,
                "analysis": analysis
            }
            
        except Exception as e:
            return {
                "is_redox": None,
                "validation_steps": [f"Error during validation: {str(e)}"]
            }
    
    def _check_common_pairs(self, reaction: str) -> Optional[FavorabilityResult]:
        """Check if reaction matches common redox pairs."""
        for name, data in COMMON_REDOX_PAIRS.items():
            if reaction in data["net_reaction"]:
                steps = [
                    f"1. Identified reaction: {name} (common redox pair)",
                    f"2. Oxidation half-reaction: {data['oxidation']}",
                    f"3. Reduction half-reaction: {data['reduction']}",
                    f"4. Standard cell potential: E°cell = {data['e_cell']} V",
                    f"5. {'Favorable' if data['favorable'] else 'Not favorable'} reaction (E°cell {'>' if data['favorable'] else '<'} 0)"
                ]
                
                return FavorabilityResult(
                    favorable=data["favorable"],
                    message="Favorable" if data["favorable"] else "Not favorable",
                    steps=steps,
                    oxidation_half=data["oxidation"],
                    reduction_half=data["reduction"],
                    e_cell=data["e_cell"]
                )
        
        return None
    
    def _analyze_common_reaction(self, parsed: Dict) -> FavorabilityResult:
        """Analyze common reactions."""
        common_name = parsed.get("common_name")
        reaction_data = get_common_redox_reaction(common_name)
        
        reduction_half = reaction_data["reduction_half"]
        oxidation_half = reaction_data["oxidation_half"]
        
        reduction_potential = self.potential_calc.find_standard_reduction_potential(reduction_half)
        oxidation_potential = self.potential_calc.find_standard_reduction_potential(
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
            
            return FavorabilityResult(
                favorable=is_favorable,
                message="Favorable" if is_favorable else "Not favorable",
                steps=steps,
                oxidation_half=oxidation_half,
                reduction_half=reduction_half,
                e_cell=e_cell
            )
        
        return FavorabilityResult(
            favorable=None,
            message="Insufficient data for common reaction",
            steps=["Could not find standard reduction potentials for this common reaction"]
        )
    
    def _is_metal_displacement(self, parsed: Dict) -> bool:
        """Check if this is a metal displacement reaction."""
        return parsed.get("oxidizing_agent") and parsed.get("reducing_agent")
    
    def _analyze_metal_displacement(self, parsed: Dict, reaction: str) -> FavorabilityResult:
        """Analyze metal displacement reactions."""
        reducing_metal = ''.join([c for c in parsed["reducing_agent"] if c.isalpha()])
        oxidizing_ion = ''.join([c for c in parsed["oxidizing_agent"] if c.isalpha()])
        
        # Get activity series positions
        red_position, oxd_position = self.activity_analyzer.get_activity_positions(reducing_metal, oxidizing_ion)
        
        if red_position is not None and oxd_position is not None:
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
            
            return FavorabilityResult(
                favorable=is_favorable,
                message="Favorable" if is_favorable else "Not favorable",
                steps=steps,
                oxidizing_agent=parsed["oxidizing_agent"],
                reducing_agent=parsed["reducing_agent"],
                oxidation_half=parsed["oxidation_half"],
                reduction_half=parsed["reduction_half"]
            )
        
        return FavorabilityResult(
            favorable=None,
            message="Metals not found in activity series",
            steps=["Could not find both metals in the activity series"]
        )
    
    def _analyze_by_potentials(self, parsed: Dict) -> FavorabilityResult:
        """Analyze using standard reduction potentials."""
        oxidation_half = parsed["oxidation_half"]
        reduction_half = parsed["reduction_half"]
        
        if oxidation_half.startswith("[") or reduction_half.startswith("["):
            return FavorabilityResult(
                favorable=None,
                message="Insufficient information",
                steps=[
                    "Could not automatically determine half-reactions.",
                    "Please provide the half-reactions or standard reduction potentials."
                ]
            )
        
        e_oxidation = self.potential_calc.find_standard_reduction_potential(
            oxidation_half.replace("->", "←").replace("+", "-").replace("2e⁻", "")
        )
        e_reduction = self.potential_calc.find_standard_reduction_potential(reduction_half)
        
        if e_oxidation is None or e_reduction is None:
            steps = [
                "Could not automatically determine standard reduction potentials.",
                "Please provide the standard reduction potentials for the half-reactions:",
                f"Oxidation half-reaction: {oxidation_half}",
                f"Reduction half-reaction: {reduction_half}"
            ]
            
            return FavorabilityResult(
                favorable=None,
                message="Insufficient information",
                steps=steps,
                oxidation_half=oxidation_half,
                reduction_half=reduction_half
            )
        
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
        
        return FavorabilityResult(
            favorable=is_favorable,
            message="Favorable" if is_favorable else "Not favorable",
            steps=steps,
            oxidation_half=oxidation_half,
            reduction_half=reduction_half,
            e_cell=e_cell
        )


# Convenience functions to maintain API compatibility
def parse_redox_reaction(reaction: str) -> Dict:
    """Parse a redox reaction - backward compatibility function."""
    parser = RedoxParser()
    return parser.parse_reaction(reaction)


def determine_redox_favorability(reaction: str) -> Dict:
    """Determine redox favorability - backward compatibility function."""
    analyzer = FavorabilityAnalyzer()
    result = analyzer.determine_favorability(reaction)
    
    # Convert to dict format for backward compatibility
    return {
        "favorable": result.favorable,
        "message": result.message,
        "steps": result.steps,
        "oxidation_half": result.oxidation_half,
        "reduction_half": result.reduction_half,
        "e_cell": result.e_cell,
        "oxidizing_agent": result.oxidizing_agent,
        "reducing_agent": result.reducing_agent,
        "balanced_equation": result.balanced_equation,
        "notes": result.notes
    }


def identify_redox_from_oxidation_states(reaction: str) -> Dict:
    """Identify redox from oxidation states - backward compatibility function."""
    analyzer = OxidationStateAnalyzer()
    analysis = analyzer.analyze_reaction(reaction)
    
    return {
        "oxidized_elements": [
            {
                "element": elem.element,
                "from_state": elem.from_state,
                "to_state": elem.to_state,
                "change": elem.change,
                "reactant_compound": elem.reactant_compound,
                "product_compound": elem.product_compound
            }
            for elem in analysis.oxidized_elements
        ],
        "reduced_elements": [
            {
                "element": elem.element,
                "from_state": elem.from_state,
                "to_state": elem.to_state,
                "change": elem.change,
                "reactant_compound": elem.reactant_compound,
                "product_compound": elem.product_compound
            }
            for elem in analysis.reduced_elements
        ],
        "oxidation_analysis": analysis.oxidation_analysis,
        "is_redox": analysis.is_redox
    }


def calculate_nernst_equation(e_cell_standard: float, n: int, concentration_ratio: float, temperature: float = 298.15) -> float:
    """Calculate Nernst equation - backward compatibility function."""
    calculator = PotentialCalculator()
    return calculator.calculate_nernst_equation(e_cell_standard, n, concentration_ratio, temperature)


def find_standard_reduction_potential(half_reaction: str) -> Optional[float]:
    """Find standard reduction potential - backward compatibility function."""
    calculator = PotentialCalculator()
    return calculator.find_standard_reduction_potential(half_reaction)


def get_half_reaction_by_element(element: str, state: Optional[str] = None) -> List[Dict]:
    """Find half-reactions by element - backward compatibility function."""
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