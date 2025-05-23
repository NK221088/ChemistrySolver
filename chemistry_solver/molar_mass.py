# Enhanced Molar Mass Calculator with Avogadro's Number Support
try:
    import molmass
except ImportError:
    molmass = None

def check_molmass():
    if molmass is None:
        raise ImportError("molmass is required. Install it with `pip install molmass`.")

# Constants
AVOGADROS_NUMBER = 6.022e23

def calculate_molar_mass(formula):
    """Calculate molar mass and composition of a chemical formula."""
    try:
        check_molmass()
        f = molmass.Formula(formula)
        comp = [{
            'element': s,
            'count': i.count,
            'atomic_mass': i.mass / i.count,
            'contribution': i.mass
        } for s, i in f.composition().items()]

        return {
            'formula': str(f),
            'hill_notation': f.formula,
            'molar_mass': f.mass,
            'monoisotopic_mass': getattr(f.isotope, 'mass', f.mass),
            'composition': comp,
            'success': True
        }
    except Exception as e:
        return {'success': False, 'error': str(e)}

def count_atoms_in_formula(formula, target_element=None):
    """Count total atoms or specific element atoms in a formula."""
    try:
        check_molmass()
        f = molmass.Formula(formula)
        composition = f.composition()
        
        if target_element:
            # Count specific element
            target_element = target_element.capitalize()
            if target_element in composition:
                return {
                    'element': target_element,
                    'count': composition[target_element].count,
                    'success': True
                }
            else:
                return {
                    'success': False, 
                    'error': f"Element {target_element} not found in {formula}"
                }
        else:
            # Count all atoms
            total_atoms = sum(item.count for item in composition.values())
            atom_breakdown = {element: item.count for element, item in composition.items()}
            return {
                'total_atoms': total_atoms,
                'atom_breakdown': atom_breakdown,
                'success': True
            }
    except Exception as e:
        return {'success': False, 'error': str(e)}

def calculate_atoms_from_mass(formula, mass_grams, target_element=None):
    """Calculate number of atoms from mass in grams."""
    try:
        # Get molar mass
        molar_result = calculate_molar_mass(formula)
        if not molar_result['success']:
            return molar_result
        
        molar_mass = molar_result['molar_mass']
        
        # Calculate moles
        moles = mass_grams / molar_mass
        
        # Get atom count per molecule
        atom_result = count_atoms_in_formula(formula, target_element)
        if not atom_result['success']:
            return atom_result
        
        if target_element:
            atoms_per_molecule = atom_result['count']
            element_name = target_element
        else:
            atoms_per_molecule = atom_result['total_atoms']
            element_name = "total"
        
        # Calculate total atoms
        total_atoms = moles * AVOGADROS_NUMBER * atoms_per_molecule
        
        return {
            'formula': formula,
            'mass_grams': mass_grams,
            'molar_mass': molar_mass,
            'moles': moles,
            'atoms_per_molecule': atoms_per_molecule,
            'total_atoms': total_atoms,
            'element': element_name,
            'success': True
        }
    except Exception as e:
        return {'success': False, 'error': str(e)}

def calculate_mass_from_atoms(formula, num_atoms, target_element=None):
    """Calculate mass in grams from number of atoms."""
    try:
        # Get molar mass
        molar_result = calculate_molar_mass(formula)
        if not molar_result['success']:
            return molar_result
        
        molar_mass = molar_result['molar_mass']
        
        # Get atom count per molecule
        atom_result = count_atoms_in_formula(formula, target_element)
        if not atom_result['success']:
            return atom_result
        
        if target_element:
            atoms_per_molecule = atom_result['count']
            element_name = target_element
        else:
            atoms_per_molecule = atom_result['total_atoms']
            element_name = "total"
        
        # Calculate moles and mass
        molecules = num_atoms / atoms_per_molecule
        moles = molecules / AVOGADROS_NUMBER
        mass_grams = moles * molar_mass
        
        return {
            'formula': formula,
            'num_atoms': num_atoms,
            'element': element_name,
            'atoms_per_molecule': atoms_per_molecule,
            'molecules': molecules,
            'moles': moles,
            'molar_mass': molar_mass,
            'mass_grams': mass_grams,
            'success': True
        }
    except Exception as e:
        return {'success': False, 'error': str(e)}

def solve_isopropanol_problem(mass_grams=36.25):
    """Solve the specific isopropanol problem from the question."""
    # Isopropanol formula: C3H8O
    formula = "C3H8O"
    
    result = calculate_atoms_from_mass(formula, mass_grams, "H")
    
    if result['success']:
        return {
            'problem': f"How many hydrogen atoms in {mass_grams}g of isopropanol?",
            'formula': formula,
            'mass_grams': mass_grams,
            'molar_mass': result['molar_mass'],
            'moles': result['moles'],
            'hydrogen_atoms_per_molecule': result['atoms_per_molecule'],
            'total_hydrogen_atoms': result['total_atoms'],
            'scientific_notation': f"{result['total_atoms']:.2e}",
            'success': True
        }
    else:
        return result

# Test the isopropanol problem
if __name__ == "__main__":
    print("Testing Isopropanol Problem:")
    print("=" * 50)
    
    result = solve_isopropanol_problem()
    
    if result['success']:
        print(f"Problem: {result['problem']}")
        print(f"Formula: {result['formula']}")
        print(f"Mass: {result['mass_grams']} g")
        print(f"Molar Mass: {result['molar_mass']:.4f} g/mol")
        print(f"Moles: {result['moles']:.6f} mol")
        print(f"H atoms per molecule: {result['hydrogen_atoms_per_molecule']}")
        print(f"Total H atoms: {result['total_hydrogen_atoms']:.2e}")
        print(f"Answer: {result['scientific_notation']} H atoms")
        
        # Check which answer choice this matches
        choices = [
            ("A", 2.90e12),
            ("B", 2.40e12), 
            ("C", 2.90e24),
            ("D", 5.80e24),
            ("E", 4.80e12)
        ]
        
        print("\nAnswer choices comparison:")
        for letter, value in choices:
            diff = abs(result['total_hydrogen_atoms'] - value)
            match = "✓" if diff < value * 0.01 else ""  # Within 1%
            print(f"{letter}: {value:.2e} {match}")
    else:
        print(f"Error: {result['error']}")