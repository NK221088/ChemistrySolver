"""
Refactored Colligative Properties Calculator

This module provides functions to solve colligative property problems:
1. Calculate molecular weight from colligative property data
2. Compare solutions based on colligative properties
3. Calculate colligative property values from concentration data

Supports: freezing point depression, boiling point elevation, osmotic pressure, vapor pressure lowering
"""

from chemistry_solver.molar_mass import calculate_molar_mass
from chemistry_solver.name_to_formula import get_formula_from_name
from .colligative_constants import SOLVENT_CONSTANTS, R_OSMOTIC, GRAMS_TO_KG, CELSIUS_TO_KELVIN


class ColligativePropertyCalculator:
    """Unified calculator for all colligative property problems."""
    
    def __init__(self, solvent="water"):
        """
        Initialize calculator with solvent properties.
        
        Parameters:
        -----------
        solvent : str
            Solvent name (default: "water")
        """
        self.solvent = solvent.lower()
        self.constants = SOLVENT_CONSTANTS.get(self.solvent, SOLVENT_CONSTANTS["water"])
        
    def _get_solvent_molar_mass(self, solvent_formula):
        """Get molar mass of solvent from formula or name."""
        # Convert name to formula if necessary
        if not any(char.isdigit() for char in solvent_formula):
            name_result = get_formula_from_name(solvent_formula)
            if name_result['success']:
                solvent_formula = name_result['formula']
        
        result = calculate_molar_mass(solvent_formula)
        if not result['success']:
            raise ValueError(f"Cannot calculate molar mass for {solvent_formula}: {result['error']}")
        
        return result['molar_mass']
    
    def calculate_property_from_molality(self, molality, property_type, ionization_factor=1):
        """
        Calculate colligative property value from molality.
        
        Parameters:
        -----------
        molality : float
            Molality of solution (mol/kg)
        property_type : str
            "freezing_point_depression", "boiling_point_elevation"
        ionization_factor : float
            van 't Hoff factor (default: 1)
        
        Returns:
        --------
        float
            Property value in °C
        """
        if property_type == "freezing_point_depression":
            return self.constants["Kf"] * molality * ionization_factor
        elif property_type == "boiling_point_elevation":
            return self.constants["Kb"] * molality * ionization_factor
        else:
            raise ValueError(f"Unsupported property type: {property_type}")
    
    def calculate_molecular_weight(self, method, delta_T, solute_mass, solvent_mass, 
                                 ionization_factor=1, answer_choices=None, **kwargs):
        """
        Calculate molecular weight from colligative property data.
        
        Parameters:
        -----------
        method : str
            "freezing_point", "boiling_point", "osmotic_pressure", "vapor_pressure"
        delta_T : float
            Temperature change (°C) or pressure for osmotic/vapor pressure methods
        solute_mass : float
            Mass of solute (g)
        solvent_mass : float
            Mass of solvent (g) - not used for osmotic pressure
        ionization_factor : float
            van 't Hoff factor (default: 1)
        answer_choices : list, optional
            Multiple choice options
        **kwargs : dict
            Method-specific parameters
        
        Returns:
        --------
        dict
            Results with molecular weight and calculation steps
        """
        if method == "freezing_point":
            return self._calculate_mw_freezing_point(
                delta_T, solute_mass, solvent_mass, ionization_factor, answer_choices
            )
        elif method == "boiling_point":
            return self._calculate_mw_boiling_point(
                delta_T, solute_mass, solvent_mass, ionization_factor, answer_choices
            )
        elif method == "osmotic_pressure":
            return self._calculate_mw_osmotic_pressure(
                delta_T, solute_mass, kwargs.get('temperature_c'), 
                kwargs.get('solution_volume_L'), ionization_factor, answer_choices
            )
        elif method == "vapor_pressure":
            return self._calculate_mw_vapor_pressure(
                kwargs.get('P_pure'), kwargs.get('P_solution'), 
                solute_mass, kwargs.get('solvent_formula'), solvent_mass, answer_choices
            )
        else:
            raise ValueError(f"Unsupported method: {method}")
    
    def _calculate_mw_freezing_point(self, delta_T, solute_mass, solvent_mass, 
                                   ionization_factor, answer_choices):
        """Calculate MW from freezing point depression."""
        solvent_mass_kg = solvent_mass / GRAMS_TO_KG
        Kf = self.constants["Kf"]
        
        # Calculate molality: m = ΔT / (Kf × i)
        molality = delta_T / (Kf * ionization_factor)
        
        # Calculate moles of solute
        moles_solute = molality * solvent_mass_kg
        
        # Calculate molecular weight
        molecular_weight = solute_mass / moles_solute
        
        steps = [
            f"Freezing Point Depression Calculation:",
            f"Given: ΔTf = {delta_T}°C, Kf = {Kf}°C/m, i = {ionization_factor}",
            f"Mass of solute = {solute_mass} g, Mass of solvent = {solvent_mass} g",
            f"",
            f"Step 1: Calculate molality",
            f"m = ΔTf / (Kf × i) = {delta_T} / ({Kf} × {ionization_factor}) = {molality:.4f} mol/kg",
            f"",
            f"Step 2: Calculate moles of solute",
            f"moles = m × kg solvent = {molality:.4f} × {solvent_mass_kg} = {moles_solute:.6f} mol",
            f"",
            f"Step 3: Calculate molecular weight",
            f"MW = mass / moles = {solute_mass} / {moles_solute:.6f} = {molecular_weight:.2f} g/mol"
        ]
        
        return self._format_mw_result(molecular_weight, steps, answer_choices, 
                                    molality=molality, moles_solute=moles_solute)
    
    def _calculate_mw_boiling_point(self, delta_T, solute_mass, solvent_mass, 
                                  ionization_factor, answer_choices):
        """Calculate MW from boiling point elevation."""
        solvent_mass_kg = solvent_mass / GRAMS_TO_KG
        Kb = self.constants["Kb"]
        
        # Calculate molality: m = ΔT / (Kb × i)
        molality = delta_T / (Kb * ionization_factor)
        
        # Calculate moles of solute
        moles_solute = molality * solvent_mass_kg
        
        # Calculate molecular weight
        molecular_weight = solute_mass / moles_solute
        
        steps = [
            f"Boiling Point Elevation Calculation:",
            f"Given: ΔTb = {delta_T}°C, Kb = {Kb}°C/m, i = {ionization_factor}",
            f"Mass of solute = {solute_mass} g, Mass of solvent = {solvent_mass} g",
            f"",
            f"Step 1: Calculate molality",
            f"m = ΔTb / (Kb × i) = {delta_T} / ({Kb} × {ionization_factor}) = {molality:.4f} mol/kg",
            f"",
            f"Step 2: Calculate moles of solute",
            f"moles = m × kg solvent = {molality:.4f} × {solvent_mass_kg} = {moles_solute:.6f} mol",
            f"",
            f"Step 3: Calculate molecular weight",
            f"MW = mass / moles = {solute_mass} / {moles_solute:.6f} = {molecular_weight:.2f} g/mol"
        ]
        
        return self._format_mw_result(molecular_weight, steps, answer_choices,
                                    molality=molality, moles_solute=moles_solute)
    
    def _calculate_mw_osmotic_pressure(self, osmotic_pressure_atm, solute_mass, 
                                     temperature_c, solution_volume_L, ionization_factor, answer_choices):
        """Calculate MW from osmotic pressure."""
        temperature_k = temperature_c + CELSIUS_TO_KELVIN
        
        # Calculate moles: π = iMRT → moles = πV/(iRT)
        moles_solute = osmotic_pressure_atm * solution_volume_L / (ionization_factor * R_OSMOTIC * temperature_k)
        
        # Calculate molecular weight
        molecular_weight = solute_mass / moles_solute
        
        steps = [
            f"Osmotic Pressure Calculation:",
            f"Given: π = {osmotic_pressure_atm} atm, T = {temperature_c}°C = {temperature_k} K",
            f"V = {solution_volume_L} L, mass = {solute_mass} g, i = {ionization_factor}",
            f"",
            f"Step 1: Calculate moles using π = iMRT",
            f"moles = πV/(iRT) = ({osmotic_pressure_atm} × {solution_volume_L}) / ({ionization_factor} × {R_OSMOTIC} × {temperature_k})",
            f"moles = {moles_solute:.6f} mol",
            f"",
            f"Step 2: Calculate molecular weight",
            f"MW = mass / moles = {solute_mass} / {moles_solute:.6f} = {molecular_weight:.2f} g/mol"
        ]
        
        return self._format_mw_result(molecular_weight, steps, answer_choices, moles_solute=moles_solute)
    
    def _calculate_mw_vapor_pressure(self, P_pure, P_solution, solute_mass, 
                                   solvent_formula, solvent_mass, answer_choices):
        """Calculate MW from vapor pressure lowering."""
        # Get solvent molar mass
        solvent_molar_mass = self._get_solvent_molar_mass(solvent_formula)
        
        # Calculate vapor pressure lowering and mole fraction
        delta_P = P_pure - P_solution
        mole_fraction_solute = delta_P / P_pure
        
        # Calculate moles of solvent and solute
        moles_solvent = solvent_mass / solvent_molar_mass
        moles_solute = (mole_fraction_solute * moles_solvent) / (1 - mole_fraction_solute)
        
        # Calculate molecular weight
        molecular_weight = solute_mass / moles_solute
        
        steps = [
            f"Vapor Pressure Lowering Calculation:",
            f"Given: P_pure = {P_pure}, P_solution = {P_solution}",
            f"Solvent: {solvent_formula} (MW = {solvent_molar_mass:.2f} g/mol)",
            f"Mass of solute = {solute_mass} g, Mass of solvent = {solvent_mass} g",
            f"",
            f"Step 1: Calculate vapor pressure lowering",
            f"ΔP = {P_pure} - {P_solution} = {delta_P}",
            f"",
            f"Step 2: Calculate mole fraction of solute",
            f"X_solute = ΔP / P_pure = {delta_P} / {P_pure} = {mole_fraction_solute:.6f}",
            f"",
            f"Step 3: Calculate moles of solvent",
            f"n_solvent = {solvent_mass} / {solvent_molar_mass:.2f} = {moles_solvent:.6f} mol",
            f"",
            f"Step 4: Calculate moles of solute",
            f"n_solute = (X_solute × n_solvent) / (1 - X_solute)",
            f"n_solute = ({mole_fraction_solute:.6f} × {moles_solvent:.6f}) / (1 - {mole_fraction_solute:.6f}) = {moles_solute:.6f} mol",
            f"",
            f"Step 5: Calculate molecular weight",
            f"MW = mass / moles = {solute_mass} / {moles_solute:.6f} = {molecular_weight:.2f} g/mol"
        ]
        
        return self._format_mw_result(molecular_weight, steps, answer_choices,
                                    moles_solute=moles_solute, moles_solvent=moles_solvent)
    
    def _format_mw_result(self, molecular_weight, steps, answer_choices, **extra_data):
        """Format molecular weight calculation result."""
        result = {
            'success': True,
            'molecular_weight': molecular_weight,
            'steps': steps,
            **extra_data
        }
        
        if answer_choices:
            closest_answer = min(answer_choices, key=lambda x: abs(x - molecular_weight))
            difference = abs(closest_answer - molecular_weight)
            
            steps.extend([
                f"",
                f"Multiple Choice Analysis:",
                f"Answer choices: {answer_choices}",
                f"Closest match: {closest_answer} g/mol (difference: {difference:.2f})"
            ])
            
            result.update({
                'closest_answer': closest_answer,
                'answer_difference': difference
            })
        
        return result
    
    def compare_solutions(self, solutions, property_type="freezing_point_depression"):
        """
        Compare colligative properties of multiple solutions.
        
        Parameters:
        -----------
        solutions : list of dict
            Each dict contains: 'molality', 'name' (optional), 'ionization_factor' (optional)
        property_type : str
            "freezing_point_depression" or "boiling_point_elevation"
        
        Returns:
        --------
        dict
            Comparison results with rankings
        """
        if property_type not in ["freezing_point_depression", "boiling_point_elevation"]:
            return {'success': False, 'error': f"Unsupported property type: {property_type}"}
        
        results = []
        
        # Calculate property for each solution
        for i, solution in enumerate(solutions):
            molality = solution['molality']
            ionization_factor = solution.get('ionization_factor', 1)
            name = solution.get('name', f"Solution {i+1} ({molality} m)")
            
            property_value = self.calculate_property_from_molality(
                molality, property_type, ionization_factor
            )
            
            results.append({
                'name': name,
                'molality': molality,
                'ionization_factor': ionization_factor,
                'property_value': property_value
            })
        
        # Sort by property value (descending for highest effect)
        results.sort(key=lambda x: x['property_value'], reverse=True)
        
        # Add rankings
        for i, result in enumerate(results):
            result['rank'] = i + 1
        
        # Generate explanation
        constant_name = "Kf" if property_type == "freezing_point_depression" else "Kb"
        constant_value = self.constants["Kf"] if property_type == "freezing_point_depression" else self.constants["Kb"]
        formula = "ΔTf = Kf × m × i" if property_type == "freezing_point_depression" else "ΔTb = Kb × m × i"
        
        steps = [
            f"Comparing {property_type.replace('_', ' ')} for multiple solutions:",
            f"Solvent: {self.solvent} ({constant_name} = {constant_value}°C/m)",
            f"Formula: {formula}",
            f"",
            f"Calculations:"
        ]
        
        for result in results:
            steps.append(f"  {result['name']}: {constant_value} × {result['molality']} × {result['ionization_factor']} = {result['property_value']:.3f}°C")
        
        steps.extend([
            f"",
            f"Ranking (highest to lowest):"
        ])
        
        for result in results:
            steps.append(f"  {result['rank']}. {result['name']}: {result['property_value']:.3f}°C")
        
        steps.extend([
            f"",
            f"Answer: {results[0]['name']} has the highest {property_type.replace('_', ' ')}"
        ])
        
        return {
            'success': True,
            'property_type': property_type,
            'solvent': self.solvent,
            'constant': constant_value,
            'results': results,
            'highest': results[0],
            'steps': steps
        }


# Convenience functions for backward compatibility and ease of use
def calculate_molecular_weight(method, solvent="water", **kwargs):
    """
    Calculate molecular weight from colligative property data.
    
    Parameters:
    -----------
    method : str
        "freezing_point", "boiling_point", "osmotic_pressure", "vapor_pressure"
    solvent : str
        Solvent name (default: "water")
    **kwargs : dict
        Method-specific parameters (see ColligativePropertyCalculator.calculate_molecular_weight)
    
    Returns:
    --------
    dict
        Results with molecular weight and steps
    """
    calculator = ColligativePropertyCalculator(solvent)
    return calculator.calculate_molecular_weight(method, **kwargs)

def compare_colligative_properties(solutions, property_type="freezing_point_depression", solvent="water"):
    """
    Compare colligative properties of multiple solutions.
    
    Parameters:
    -----------
    solutions : list of dict
        Each dict contains: 'molality', 'name' (optional), 'ionization_factor' (optional)
    property_type : str
        "freezing_point_depression" or "boiling_point_elevation"
    solvent : str
        Solvent name (default: "water")
    
    Returns:
    --------
    dict
        Comparison results
    """
    calculator = ColligativePropertyCalculator(solvent)
    return calculator.compare_solutions(solutions, property_type)

def solve_multiple_choice_problem(problem_data):
    """
    Solve various types of colligative property problems.
    
    Parameters:
    -----------
    problem_data : dict
        Problem configuration with keys:
        - 'type': 'molecular_weight' or 'compare_solutions'
        - 'method': calculation method (for molecular weight problems)
        - 'solvent': solvent name (default: "water")
        - Other method-specific parameters
    
    Returns:
    --------
    dict
        Complete solution
    """
    problem_type = problem_data.get('type')
    solvent = problem_data.get('solvent', 'water')
    calculator = ColligativePropertyCalculator(solvent)
    
    if problem_type == 'molecular_weight':
        method = problem_data.get('method')
        if not method:
            return {'success': False, 'error': 'Method is required for molecular weight problems'}
        
        # Remove non-method parameters
        method_params = {k: v for k, v in problem_data.items() 
                        if k not in ['type', 'method', 'solvent']}
        
        return calculator.calculate_molecular_weight(method, **method_params)
    
    elif problem_type == 'compare_solutions':
        solutions = problem_data.get('solutions')
        property_type = problem_data.get('property_type', 'freezing_point_depression')
        
        if not solutions:
            return {'success': False, 'error': 'Solutions list is required'}
        
        return calculator.compare_solutions(solutions, property_type)
    
    else:
        return {'success': False, 'error': f'Unsupported problem type: {problem_type}'}


# Example usage
if __name__ == "__main__":
    # Example 1: Compare solutions for freezing point depression
    print("=== SOLUTION COMPARISON EXAMPLE ===")
    solutions = [
        {'molality': 2.6, 'name': '2.6 m solution'},
        {'molality': 3.3, 'name': '3.3 m solution'},
        {'molality': 1.1, 'name': '1.1 m solution'},
        {'molality': 5.7, 'name': '5.7 m solution'},
        {'molality': 4.4, 'name': '4.4 m solution'}
    ]
    
    result = compare_colligative_properties(solutions, "freezing_point_depression", "water")
    if result['success']:
        for step in result['steps']:
            print(step)
    
    print("\n" + "="*50 + "\n")
    
    # Example 2: Calculate molecular weight from freezing point data
    print("=== MOLECULAR WEIGHT CALCULATION EXAMPLE ===")
    mw_result = calculate_molecular_weight(
        method="freezing_point",
        solvent="water",
        delta_T=1.5,  # 1.5°C depression
        solute_mass=2.0,  # 2.0 g solute
        solvent_mass=100.0,  # 100 g water
        answer_choices=[46, 56, 66, 76, 86]
    )
    
    if mw_result['success']:
        for step in mw_result['steps']:
            print(step)