"""
Ionic Character Calculator Module for Chemistry Problem Solver
"""
import re

class IonicCharacterCalculator:
    """
    Calculator for determining ionic character of chemical bonds and compounds
    """
    
    # Electronegativity values (Pauling scale)
    ELECTRONEGATIVITY = {
        'H': 2.20, 'He': 0.00,
        'Li': 0.98, 'Be': 1.57, 'B': 2.04, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98, 'Ne': 0.00,
        'Na': 0.93, 'Mg': 1.31, 'Al': 1.61, 'Si': 1.90, 'P': 2.19, 'S': 2.58, 'Cl': 3.16, 'Ar': 0.00,
        'K': 0.82, 'Ca': 1.00, 'Sc': 1.36, 'Ti': 1.54, 'V': 1.63, 'Cr': 1.66, 'Mn': 1.55, 'Fe': 1.83,
        'Co': 1.88, 'Ni': 1.91, 'Cu': 1.90, 'Zn': 1.65, 'Ga': 1.81, 'Ge': 2.01, 'As': 2.18, 'Se': 2.55,
        'Br': 2.96, 'Kr': 3.00, 'Rb': 0.82, 'Sr': 0.95, 'Y': 1.22, 'Zr': 1.33, 'Nb': 1.6, 'Mo': 2.16,
        'Tc': 1.9, 'Ru': 2.2, 'Rh': 2.28, 'Pd': 2.20, 'Ag': 1.93, 'Cd': 1.69, 'In': 1.78, 'Sn': 1.96,
        'Sb': 2.05, 'Te': 2.1, 'I': 2.66, 'Xe': 2.60, 'Cs': 0.79, 'Ba': 0.89, 'La': 1.10, 'Ce': 1.12,
        'Pr': 1.13, 'Nd': 1.14, 'Pm': 1.13, 'Sm': 1.17, 'Eu': 1.20, 'Gd': 1.20, 'Tb': 1.10, 'Dy': 1.22,
        'Ho': 1.23, 'Er': 1.24, 'Tm': 1.25, 'Yb': 1.10, 'Lu': 1.27, 'Hf': 1.3, 'Ta': 1.5, 'W': 2.36,
        'Re': 1.9, 'Os': 2.2, 'Ir': 2.20, 'Pt': 2.28, 'Au': 2.54, 'Hg': 2.00, 'Tl': 1.62, 'Pb': 2.33,
        'Bi': 2.02, 'Po': 2.0, 'At': 2.2, 'Rn': 2.2, 'Fr': 0.7, 'Ra': 0.9, 'Ac': 1.1, 'Th': 1.3,
        'Pa': 1.5, 'U': 1.38, 'Np': 1.36, 'Pu': 1.28, 'Am': 1.13, 'Cm': 1.28, 'Bk': 1.3, 'Cf': 1.3,
        'Es': 1.3, 'Fm': 1.3, 'Md': 1.3, 'No': 1.3, 'Lr': 1.3
    }
    
    # Common polyatomic ions
    POLYATOMIC_CATIONS = {
        'NH4': {'name': 'ammonium', 'charge': 1, 'formula': 'NH4+'},
        'H3O': {'name': 'hydronium', 'charge': 1, 'formula': 'H3O+'},
        'NO': {'name': 'nitrosyl', 'charge': 1, 'formula': 'NO+'},
        'NO2': {'name': 'nitronium', 'charge': 1, 'formula': 'NO2+'}
    }
    
    POLYATOMIC_ANIONS = {
        'OH': {'name': 'hydroxide', 'charge': -1, 'formula': 'OH-'},
        'NO3': {'name': 'nitrate', 'charge': -1, 'formula': 'NO3-'},
        'NO2': {'name': 'nitrite', 'charge': -1, 'formula': 'NO2-'},
        'SO4': {'name': 'sulfate', 'charge': -2, 'formula': 'SO4²-'},
        'SO3': {'name': 'sulfite', 'charge': -2, 'formula': 'SO3²-'},
        'PO4': {'name': 'phosphate', 'charge': -3, 'formula': 'PO4³-'},
        'CO3': {'name': 'carbonate', 'charge': -2, 'formula': 'CO3²-'},
        'HCO3': {'name': 'bicarbonate', 'charge': -1, 'formula': 'HCO3-'},
        'ClO4': {'name': 'perchlorate', 'charge': -1, 'formula': 'ClO4-'},
        'ClO3': {'name': 'chlorate', 'charge': -1, 'formula': 'ClO3-'},
        'ClO2': {'name': 'chlorite', 'charge': -1, 'formula': 'ClO2-'},
        'ClO': {'name': 'hypochlorite', 'charge': -1, 'formula': 'ClO-'},
        'MnO4': {'name': 'permanganate', 'charge': -1, 'formula': 'MnO4-'},
        'Cr2O7': {'name': 'dichromate', 'charge': -2, 'formula': 'Cr2O7²-'},
        'CrO4': {'name': 'chromate', 'charge': -2, 'formula': 'CrO4²-'},
        'C2H3O2': {'name': 'acetate', 'charge': -1, 'formula': 'C2H3O2-'},
        'CN': {'name': 'cyanide', 'charge': -1, 'formula': 'CN-'},
        'SCN': {'name': 'thiocyanate', 'charge': -1, 'formula': 'SCN-'},
        'OCN': {'name': 'cyanate', 'charge': -1, 'formula': 'OCN-'}
    }
    
    # Common monoatomic cations and their charges
    MONOATOMIC_CATIONS = {
        'Li': 1, 'Na': 1, 'K': 1, 'Rb': 1, 'Cs': 1, 'Fr': 1,
        'Be': 2, 'Mg': 2, 'Ca': 2, 'Sr': 2, 'Ba': 2, 'Ra': 2,
        'Al': 3, 'Ga': 3, 'In': 3,
        'Zn': 2, 'Cd': 2, 'Hg': 2,
        'Ag': 1, 'Cu': 1, 'Au': 1,
        'Fe': [2, 3], 'Co': [2, 3], 'Ni': 2, 'Mn': [2, 3, 4, 6, 7],
        'Cr': [2, 3, 6], 'V': [2, 3, 4, 5], 'Ti': [2, 3, 4],
        'Pb': [2, 4], 'Sn': [2, 4], 'Bi': 3
    }
    
    # Common monoatomic anions and their charges
    MONOATOMIC_ANIONS = {
        'F': -1, 'Cl': -1, 'Br': -1, 'I': -1, 'At': -1,
        'O': -2, 'S': -2, 'Se': -2, 'Te': -2, 'Po': -2,
        'N': -3, 'P': -3, 'As': -3, 'Sb': -3, 'Bi': -3
    }
    
    def __init__(self):
        pass
    
    def parse_compound(self, compound):
        """
        Parse a compound formula to extract the two main elements
        
        Args:
            compound (str): Chemical formula (e.g., 'HCl', 'NaCl', 'CaF2')
            
        Returns:
            tuple: (element1, element2) or None if parsing fails
        """
        compound = compound.strip()
        
        # Handle simple binary compounds
        elements = []
        i = 0
        while i < len(compound):
            if compound[i].isupper():
                element = compound[i]
                i += 1
                # Check for lowercase letters (part of element symbol)
                while i < len(compound) and compound[i].islower():
                    element += compound[i]
                    i += 1
                # Skip numbers
                while i < len(compound) and compound[i].isdigit():
                    i += 1
                elements.append(element)
            else:
                i += 1
        
        if len(elements) >= 2:
            return elements[0], elements[1]
        return None
    
    def get_electronegativity_difference(self, compound):
        """
        Calculate electronegativity difference for a compound
        
        Args:
            compound (str): Chemical formula
            
        Returns:
            dict: Contains electronegativity difference and other info
        """
        parsed = self.parse_compound(compound)
        if not parsed:
            return {
                'success': False,
                'error': f"Could not parse compound: {compound}"
            }
        
        element1, element2 = parsed
        
        if element1 not in self.ELECTRONEGATIVITY:
            return {
                'success': False,
                'error': f"Electronegativity not found for element: {element1}"
            }
        
        if element2 not in self.ELECTRONEGATIVITY:
            return {
                'success': False,
                'error': f"Electronegativity not found for element: {element2}"
            }
        
        en1 = self.ELECTRONEGATIVITY[element1]
        en2 = self.ELECTRONEGATIVITY[element2]
        en_diff = abs(en1 - en2)
        
        # Determine bond type based on electronegativity difference
        if en_diff < 0.5:
            bond_type = "Nonpolar covalent"
        elif en_diff < 1.7:
            bond_type = "Polar covalent"
        else:
            bond_type = "Ionic"
        
        # Calculate percent ionic character using Pauling's equation
        percent_ionic = 100 * (1 - pow(2.718, -0.25 * en_diff * en_diff))
        
        return {
            'success': True,
            'compound': compound,
            'element1': element1,
            'element2': element2,
            'electronegativity1': en1,
            'electronegativity2': en2,
            'electronegativity_difference': en_diff,
            'bond_type': bond_type,
            'percent_ionic_character': percent_ionic
        }
    
    def analyze_compound_list(self, compounds):
        """
        Analyze a list of compounds and return their ionic character data
        
        Args:
            compounds (list): List of compound formulas
            
        Returns:
            list: List of analysis results for each compound
        """
        results = []
        for compound in compounds:
            result = self.get_electronegativity_difference(compound)
            results.append(result)
        return results
    
    def order_by_ionic_character(self, compounds):
        """
        Order compounds by increasing ionic character
        
        Args:
            compounds (list): List of compound formulas
            
        Returns:
            dict: Contains ordered compounds and analysis
        """
        results = self.analyze_compound_list(compounds)
        
        # Filter successful results
        valid_results = [r for r in results if r['success']]
        failed_results = [r for r in results if not r['success']]
        
        if not valid_results:
            return {
                'success': False,
                'error': "No valid compounds could be analyzed"
            }
        
        # Sort by electronegativity difference (ionic character)
        sorted_results = sorted(valid_results, key=lambda x: x['electronegativity_difference'])
        
        return {
            'success': True,
            'ordered_compounds': [r['compound'] for r in sorted_results],
            'detailed_analysis': sorted_results,
            'failed_compounds': failed_results
        }
    
    def check_order_correctness(self, given_order, correct_order):
        """
        Check if a given order matches the correct order
        
        Args:
            given_order (list): The order to check
            correct_order (list): The correct order
            
        Returns:
            dict: Analysis of correctness
        """
        return {
            'is_correct': given_order == correct_order,
            'given_order': given_order,
            'correct_order': correct_order,
            'differences': [i for i, (g, c) in enumerate(zip(given_order, correct_order)) if g != c]
        }
    
    def solve_multiple_choice_problem(self, choices):
        """
        Solve a multiple choice problem about ionic character ordering
        
        Args:
            choices (list): List of choice lists, each containing compound formulas
            
        Returns:
            dict: Analysis of all choices and the correct answer
        """
        choice_analyses = []
        
        for i, choice in enumerate(choices):
            if isinstance(choice, str) and choice.lower() in ['none of the above', 'none']:
                choice_analyses.append({
                    'choice_number': i + 1,
                    'compounds': choice,
                    'is_none_option': True
                })
                continue
            
            # Analyze this choice
            result = self.order_by_ionic_character(choice)
            if result['success']:
                correct_order = result['ordered_compounds']
                is_correct_order = choice == correct_order
                
                choice_analyses.append({
                    'choice_number': i + 1,
                    'compounds': choice,
                    'correct_order': correct_order,
                    'is_correct_order': is_correct_order,
                    'detailed_analysis': result['detailed_analysis'],
                    'is_none_option': False
                })
            else:
                choice_analyses.append({
                    'choice_number': i + 1,
                    'compounds': choice,
                    'error': result['error'],
                    'is_none_option': False
                })
        
        # Determine if any choice is correct
        correct_choices = [c for c in choice_analyses if c.get('is_correct_order', False)]
        
        if correct_choices:
            answer = correct_choices[0]['choice_number']
        else:
            # Check if there's a "none of the above" option
            none_options = [c for c in choice_analyses if c.get('is_none_option', False)]
            if none_options:
                answer = none_options[0]['choice_number']
            else:
                answer = "None of the given choices are correct"
        
        return {
            'answer': answer,
            'choice_analyses': choice_analyses,
            'explanation': self._generate_explanation(choice_analyses)
        }
    
    def _generate_explanation(self, choice_analyses):
        """Generate explanation for the solution"""
        explanation = []
        explanation.append("Analysis of each choice:")
        
        for analysis in choice_analyses:
            if analysis.get('is_none_option'):
                explanation.append(f"\nChoice {analysis['choice_number']}: {analysis['compounds']}")
                continue
                
            explanation.append(f"\nChoice {analysis['choice_number']}: {', '.join(analysis['compounds'])}")
            
            if 'error' in analysis:
                explanation.append(f"  Error: {analysis['error']}")
                continue
            
            if analysis.get('is_correct_order'):
                explanation.append("  ✓ This order is CORRECT")
            else:
                explanation.append("  ✗ This order is INCORRECT")
                explanation.append(f"  Correct order: {', '.join(analysis['correct_order'])}")
            
            # Add detailed electronegativity analysis
            explanation.append("  Electronegativity differences:")
            for detail in analysis['detailed_analysis']:
                explanation.append(f"    {detail['compound']}: ΔEN = {detail['electronegativity_difference']:.3f}")
        
    def parse_ionic_compound(self, compound):
        """
        Parse an ionic compound to identify cations and anions
        
        Args:
            compound (str): Chemical formula (e.g., 'NH4ClO3', 'K2Cr2O7')
            
        Returns:
            dict: Contains cation and anion information
        """
        compound = compound.strip()
        
        # Try to identify polyatomic ions first
        cations = []
        anions = []
        
        # Check for polyatomic cations
        for poly_cat, info in self.POLYATOMIC_CATIONS.items():
            if poly_cat in compound:
                # Find how many of this cation
                pattern = rf'{poly_cat}(\d*)'
                match = re.search(pattern, compound)
                if match:
                    count = int(match.group(1)) if match.group(1) else 1
                    cations.append({
                        'formula': poly_cat,
                        'name': info['name'],
                        'charge': info['charge'],
                        'count': count,
                        'type': 'polyatomic'
                    })
        
        # Check for polyatomic anions
        for poly_an, info in self.POLYATOMIC_ANIONS.items():
            if poly_an in compound:
                # Find how many of this anion
                pattern = rf'{poly_an}(\d*)'
                match = re.search(pattern, compound)
                if match:
                    count = int(match.group(1)) if match.group(1) else 1
                    anions.append({
                        'formula': poly_an,
                        'name': info['name'],
                        'charge': info['charge'],
                        'count': count,
                        'type': 'polyatomic'
                    })
        
        # If no polyatomic ions found, try to parse as simple ionic compound
        if not cations and not anions:
            # Simple parsing for monoatomic ions
            elements = self._extract_elements(compound)
            
            for element, count in elements:
                if element in self.MONOATOMIC_CATIONS:
                    charge = self.MONOATOMIC_CATIONS[element]
                    if isinstance(charge, list):
                        charge = charge[0]  # Use most common charge
                    cations.append({
                        'formula': element,
                        'name': element.lower(),
                        'charge': charge,
                        'count': count,
                        'type': 'monoatomic'
                    })
                elif element in self.MONOATOMIC_ANIONS:
                    anions.append({
                        'formula': element,
                        'name': element.lower() + 'ide',
                        'charge': self.MONOATOMIC_ANIONS[element],
                        'count': count,
                        'type': 'monoatomic'
                    })
        
        return {
            'compound': compound,
            'cations': cations,
            'anions': anions,
            'has_polyatomic_cation': any(c['type'] == 'polyatomic' for c in cations),
            'has_polyatomic_anion': any(a['type'] == 'polyatomic' for a in anions)
        }
    
    def _extract_elements(self, compound):
        """Extract elements and their counts from a compound formula"""
        pattern = r'([A-Z][a-z]?)(\d*)'
        matches = re.findall(pattern, compound)
        elements = []
        for element, count_str in matches:
            count = int(count_str) if count_str else 1
            elements.append((element, count))
        return elements
    
    def identify_polyatomic_compounds(self, compounds):
        """
        Identify which compounds contain polyatomic cations and anions
        
        Args:
            compounds (list): List of compound formulas
            
        Returns:
            dict: Analysis of each compound
        """
        results = []
        
        for compound in compounds:
            analysis = self.parse_ionic_compound(compound)
            
            # Determine if compound has both polyatomic cation and anion
            has_both_polyatomic = (analysis['has_polyatomic_cation'] and 
                                 analysis['has_polyatomic_anion'])
            
            results.append({
                'compound': compound,
                'analysis': analysis,
                'has_both_polyatomic': has_both_polyatomic,
                'polyatomic_cations': [c for c in analysis['cations'] if c['type'] == 'polyatomic'],
                'polyatomic_anions': [a for a in analysis['anions'] if a['type'] == 'polyatomic'],
                'monoatomic_cations': [c for c in analysis['cations'] if c['type'] == 'monoatomic'],
                'monoatomic_anions': [a for a in analysis['anions'] if a['type'] == 'monoatomic']
            })
        
        return results
    
    def solve_polyatomic_multiple_choice(self, choices):
        """
        Solve multiple choice questions about polyatomic ions
        
        Args:
            choices (list): List of compound choices
            
        Returns:
            dict: Analysis and answer
        """
        choice_analyses = []
        
        for i, choice in enumerate(choices):
            if isinstance(choice, str) and choice.lower() in ['none of the above', 'none']:
                choice_analyses.append({
                    'choice_number': i + 1,
                    'compound': choice,
                    'is_none_option': True
                })
                continue
            
            if isinstance(choice, list):
                # Handle multiple compounds in a choice (shouldn't happen for this type of question)
                choice = choice[0] if choice else ""
            
            analysis = self.parse_ionic_compound(choice)
            
            choice_analyses.append({
                'choice_number': i + 1,
                'compound': choice,
                'analysis': analysis,
                'has_both_polyatomic': (analysis['has_polyatomic_cation'] and 
                                      analysis['has_polyatomic_anion']),
                'polyatomic_cations': [c for c in analysis['cations'] if c['type'] == 'polyatomic'],
                'polyatomic_anions': [a for a in analysis['anions'] if a['type'] == 'polyatomic'],
                'is_none_option': False
            })
        
        # Find the correct answer
        correct_choices = [c for c in choice_analyses 
                          if c.get('has_both_polyatomic', False)]
        
        if correct_choices:
            answer = correct_choices[0]['choice_number']
        else:
            # Check if there's a "none of the above" option
            none_options = [c for c in choice_analyses if c.get('is_none_option', False)]
            if none_options:
                answer = none_options[0]['choice_number']
            else:
                answer = "None of the choices contain both polyatomic cation and anion"
        
        return {
            'answer': answer,
            'choice_analyses': choice_analyses,
            'explanation': self._generate_polyatomic_explanation(choice_analyses)
        }
    
    def _generate_polyatomic_explanation(self, choice_analyses):
        """Generate explanation for polyatomic ion problems"""
        explanation = []
        explanation.append("Analysis of each compound:")
        
        for analysis in choice_analyses:
            if analysis.get('is_none_option'):
                explanation.append(f"\nChoice {analysis['choice_number']}: {analysis['compound']}")
                continue
            
            explanation.append(f"\nChoice {analysis['choice_number']}: {analysis['compound']}")
            
            # List cations
            if analysis.get('polyatomic_cations'):
                for cation in analysis['polyatomic_cations']:
                    explanation.append(f"  ✓ Polyatomic cation: {cation['name']} ({cation['formula']})")
            if analysis.get('analysis', {}).get('cations'):
                mono_cations = [c for c in analysis['analysis']['cations'] if c['type'] == 'monoatomic']
                for cation in mono_cations:
                    explanation.append(f"  - Monoatomic cation: {cation['name']} ({cation['formula']})")
            
            # List anions
            if analysis.get('polyatomic_anions'):
                for anion in analysis['polyatomic_anions']:
                    explanation.append(f"  ✓ Polyatomic anion: {anion['name']} ({anion['formula']})")
            if analysis.get('analysis', {}).get('anions'):
                mono_anions = [a for a in analysis['analysis']['anions'] if a['type'] == 'monoatomic']
                for anion in mono_anions:
                    explanation.append(f"  - Monoatomic anion: {anion['name']} ({anion['formula']})")
            
            # Conclusion for this choice
            if analysis.get('has_both_polyatomic'):
                explanation.append("  ✓ HAS BOTH polyatomic cation AND polyatomic anion")
            else:
                explanation.append("  ✗ Does NOT have both polyatomic cation and anion")
        
        return "\n".join(explanation)