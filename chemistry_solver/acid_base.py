"""
Complete Acid-Base Chemistry Solver
Handles strong/weak acids/bases with automatic compound identification
"""
import math

def identify_acid_base(compound):
    """
    Identifies if a compound is likely an acid, base, or neutral using Brønsted definition.
    
    Args:
        compound (str): Chemical formula of the compound
    
    Returns:
        dict: Contains classification and explanation
    """
    # Clean the input
    compound = compound.strip()
    
    # Common acids dictionary
    common_acids = {
        "HCl": "Hydrochloric acid - donates H+ in solution",
        "HBr": "Hydrobromic acid - donates H+ in solution",
        "HI": "Hydroiodic acid - donates H+ in solution", 
        "HF": "Hydrofluoric acid - donates H+ in solution",
        "H2SO4": "Sulfuric acid - can donate H+ ions",
        "H2SO3": "Sulfurous acid - can donate H+ ions",
        "HNO3": "Nitric acid - donates H+ in solution",
        "HNO2": "Nitrous acid - donates H+ in solution",
        "H3PO4": "Phosphoric acid - can donate H+ ions",
        "CH3COOH": "Acetic acid - carboxylic acid that donates H+ from -COOH group",
        "HCOOH": "Formic acid - carboxylic acid that donates H+ from -COOH group",
        "H2CO3": "Carbonic acid - can donate H+ ions",
        "HClO4": "Perchloric acid - donates H+ in solution",
        "HClO3": "Chloric acid - donates H+ in solution"
    }
    
    # Common bases dictionary
    common_bases = {
        "NaOH": "Sodium hydroxide - releases OH- which accepts H+",
        "KOH": "Potassium hydroxide - releases OH- which accepts H+",
        "LiOH": "Lithium hydroxide - releases OH- which accepts H+",
        "RbOH": "Rubidium hydroxide - releases OH- which accepts H+",
        "CsOH": "Cesium hydroxide - releases OH- which accepts H+",
        "NH3": "Ammonia - accepts H+ with its lone pair",
        "Ca(OH)2": "Calcium hydroxide - releases OH- which accepts H+",
        "Sr(OH)2": "Strontium hydroxide - releases OH- which accepts H+",
        "Ba(OH)2": "Barium hydroxide - releases OH- which accepts H+",
        "Mg(OH)2": "Magnesium hydroxide - releases OH- which accepts H+",
        "NaH": "Sodium hydride - contains H- which acts as a proton acceptor",
        "KH": "Potassium hydride - contains H- which acts as a proton acceptor",
        "LiH": "Lithium hydride - contains H- which acts as a proton acceptor",
        "NaNH2": "Sodium amide - strong base",
        "NaHCO3": "Sodium bicarbonate - weak base"
    }
    
    # Common neutral compounds
    common_neutral = {
        "CH4": "Methane - C-H bonds are not acidic enough to donate protons under normal conditions",
        "C2H6": "Ethane - hydrocarbon with non-acidic C-H bonds",
        "H2O": "Water - amphoteric (can act as both acid and base)",
        "CO2": "Carbon dioxide - forms carbonic acid in water but molecule itself is neutral",
        "N2": "Nitrogen gas - inert molecule"
    }
    
    # First check if it's a common compound we know
    if compound in common_acids:
        return {"classification": "Acid", "explanation": common_acids[compound]}
    
    if compound in common_bases:
        return {"classification": "Base", "explanation": common_bases[compound]}
    
    if compound in common_neutral:
        return {"classification": "Neutral", "explanation": common_neutral[compound]}
    
    # Rules for identifying acids
    # 1. Inorganic acids often start with H
    if compound.startswith("H") and any(char.isupper() for char in compound[1:]):
        return {
            "classification": "Likely Acid", 
            "explanation": f"Compound {compound} starts with H followed by non-metals, suggesting it may donate H+ ions."
        }
    
    # 2. Carboxylic acids end with COOH
    if "COOH" in compound:
        return {
            "classification": "Acid", 
            "explanation": f"Compound {compound} contains a carboxylic acid group (-COOH) which can donate H+."
        }
    
    # Rules for identifying bases
    # 1. Metal hydroxides
    metal_symbols = ["Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba"]
    if any(compound.startswith(metal) for metal in metal_symbols) and "OH" in compound:
        return {
            "classification": "Base", 
            "explanation": f"Compound {compound} appears to be a metal hydroxide, which releases OH- ions that accept H+."
        }
    
    # 2. Metal hydrides
    if any(compound.startswith(metal) and compound[len(metal):] == "H" for metal in metal_symbols):
        return {
            "classification": "Base", 
            "explanation": f"Compound {compound} appears to be a metal hydride, which contains H- ions that accept H+."
        }
    
    # Alcohols (weak acids, often treated as neutral in intro chemistry)
    if "OH" in compound and not any(compound.startswith(metal) for metal in metal_symbols):
        return {
            "classification": "Very Weak Acid / Practically Neutral", 
            "explanation": f"Compound {compound} appears to be an alcohol with -OH group. These can technically donate H+ but are extremely weak acids."
        }
    
    # Default case - need more information
    return {
        "classification": "Unknown", 
        "explanation": f"Cannot confidently classify {compound} without more information about its structure and properties."
    }

class AcidBaseEquilibrium:
    """
    Complete class for handling all types of acid-base equilibrium calculations
    """
    
    def __init__(self):
        # Dictionary of common weak acids and their Ka values at 25°C
        self.weak_acid_ka = {
            "CH3COOH": 1.8e-5,       # Acetic acid
            "HF": 6.8e-4,            # Hydrofluoric acid
            "HNO2": 4.5e-4,          # Nitrous acid
            "HCOOH": 1.8e-4,         # Formic acid
            "C6H5COOH": 6.5e-5,      # Benzoic acid
            "H2CO3": 4.3e-7,         # Carbonic acid (first dissociation)
            "H2S": 1.0e-7,           # Hydrogen sulfide (first dissociation)
            "HCN": 4.9e-10,          # Hydrocyanic acid
            "HClO": 2.9e-8,          # Hypochlorous acid
            "HBrO": 2.5e-9,          # Hypobromous acid
            "H3BO3": 5.8e-10,        # Boric acid
            "C6H5OH": 1.3e-10,       # Phenol
            "HClO2": 1.1e-2,         # Chlorous acid
            "H3PO4": 7.5e-3,         # Phosphoric acid (first dissociation)
            "H2PO4-": 6.2e-8,        # Phosphoric acid (second dissociation)
            "HPO4-2": 4.8e-13,       # Phosphoric acid (third dissociation)
            "H2SO3": 1.5e-2,         # Sulfurous acid (first dissociation)
            "HSO3-": 6.3e-8          # Sulfurous acid (second dissociation)
        }
        
        # Dictionary of common weak bases and their Kb values at 25°C
        self.weak_base_kb = {
            "NH3": 1.8e-5,           # Ammonia
            "CH3NH2": 4.4e-4,        # Methylamine
            "C2H5NH2": 5.6e-4,       # Ethylamine
            "C6H5NH2": 4.3e-10,      # Aniline
            "C5H5N": 1.7e-9,         # Pyridine
            "NH2OH": 1.1e-8,         # Hydroxylamine
            "N2H4": 1.7e-6,          # Hydrazine
            "(CH3)2NH": 5.4e-4,      # Dimethylamine
            "(CH3)3N": 6.5e-5,       # Trimethylamine
            "NH4OH": 1.8e-5          # Ammonium hydroxide
        }
    
    def get_acid_ka(self, formula):
        """Get the Ka value for a specific acid"""
        if formula in self.weak_acid_ka:
            return self.weak_acid_ka[formula]
        return None
    
    def get_base_kb(self, formula):
        """Get the Kb value for a specific base"""
        if formula in self.weak_base_kb:
            return self.weak_base_kb[formula]
        return None
    
    def solve_strong_acid_equilibrium(self, formula, initial_concentration):
        """
        Calculate pH for a strong acid solution
        
        Args:
            formula (str): Chemical formula of the strong acid
            initial_concentration (float): Initial concentration of the acid in mol/L
            
        Returns:
            dict: Contains calculated pH and concentrations
        """
        # List of common strong acids
        strong_acids = ["HCl", "HBr", "HI", "HNO3", "H2SO4", "HClO4", "HClO3"]
        
        # Check if it's a strong acid
        compound_info = identify_acid_base(formula)
        if "Acid" not in compound_info["classification"]:
            return {
                "error": f"{formula} is not identified as an acid.",
                "classification": compound_info["classification"]
            }
        
        # For strong acids, assume complete dissociation
        # For diprotic acids like H2SO4, need to account for 2 H+ per molecule
        if formula == "H2SO4":
            h_concentration = 2 * initial_concentration
        elif formula == "H3PO4" and formula in strong_acids:  # If treating as strong
            h_concentration = 3 * initial_concentration
        else:
            h_concentration = initial_concentration
        
        # Calculate pH
        if h_concentration <= 0:
            return {"error": "Concentration must be positive"}
        
        ph = -math.log10(h_concentration)
        
        # Calculate pOH and OH- concentration
        poh = 14 - ph
        oh_concentration = 10**(-poh)
        
        return {
            "formula": formula,
            "initial_concentration": initial_concentration,
            "h_concentration": h_concentration,
            "oh_concentration": oh_concentration,
            "ph": ph,
            "poh": poh,
            "acid_type": "strong",
            "complete_dissociation": True
        }
    
    def solve_strong_base_equilibrium(self, formula, initial_concentration):
        """
        Calculate pH for a strong base solution
        
        Args:
            formula (str): Chemical formula of the strong base
            initial_concentration (float): Initial concentration of the base in mol/L
            
        Returns:
            dict: Contains calculated pH and concentrations
        """
        # List of common strong bases
        strong_bases = ["NaOH", "KOH", "LiOH", "RbOH", "CsOH", "Ca(OH)2", "Sr(OH)2", "Ba(OH)2"]
        
        # Check if it's a base
        compound_info = identify_acid_base(formula)
        if "Base" not in compound_info["classification"]:
            return {
                "error": f"{formula} is not identified as a base.",
                "classification": compound_info["classification"]
            }
        
        # For strong bases, assume complete dissociation
        # For bases with multiple OH- groups, need to account for multiple OH- per molecule
        if formula in ["Ca(OH)2", "Sr(OH)2", "Ba(OH)2", "Mg(OH)2"]:
            oh_concentration = 2 * initial_concentration
        else:
            oh_concentration = initial_concentration
        
        # Calculate pOH and pH
        if oh_concentration <= 0:
            return {"error": "Concentration must be positive"}
        
        poh = -math.log10(oh_concentration)
        ph = 14 - poh  # at 25°C, pH + pOH = 14
        
        # Calculate H+ concentration
        h_concentration = 10**(-ph)
        
        return {
            "formula": formula,
            "initial_concentration": initial_concentration,
            "oh_concentration": oh_concentration,
            "h_concentration": h_concentration,
            "ph": ph,
            "poh": poh,
            "base_type": "strong",
            "complete_dissociation": True
        }
    
    def solve_weak_acid_equilibrium(self, formula, initial_concentration, ka=None):
        """
        Calculate equilibrium concentrations and pH for a weak acid solution
        
        Args:
            formula (str): Chemical formula of the weak acid
            initial_concentration (float): Initial concentration of the acid in mol/L
            ka (float, optional): Acid dissociation constant. If None, uses value from database.
            
        Returns:
            dict: Contains calculated equilibrium values
        """
        # Get Ka value if not provided
        if ka is None:
            if formula in self.weak_acid_ka:
                ka = self.weak_acid_ka[formula]
            else:
                return {
                    "error": f"Ka value for {formula} not found in database. Please provide Ka."
                }
        
        # Identify if it's an acid
        compound_info = identify_acid_base(formula)
        if "Acid" not in compound_info["classification"]:
            return {
                "error": f"{formula} is not identified as an acid. Results may not be accurate.",
                "classification": compound_info["classification"]
            }
        
        # Solve the equilibrium problem
        # For a weak acid HA: HA ⇌ H⁺ + A⁻
        # Ka = [H⁺][A⁻]/[HA]
        
        # Using the quadratic formula to solve for x (where x = [H⁺] = [A⁻])
        # x² + Ka*x - Ka*C₀ = 0
        # x = (-Ka + √(Ka² + 4*Ka*C₀))/2
        
        a = 1
        b = ka
        c = -ka * initial_concentration
        
        discriminant = b**2 - 4*a*c
        if discriminant < 0:
            return {"error": "Calculation error: negative discriminant"}
        
        x = (-b + math.sqrt(discriminant)) / (2*a)
        
        # Calculate concentrations at equilibrium
        h_concentration = x  # [H⁺]
        a_concentration = x  # [A⁻]
        ha_concentration = initial_concentration - x  # [HA]
        
        # Calculate pH and pOH
        ph = -math.log10(h_concentration)
        poh = 14 - ph
        oh_concentration = 10**(-poh)
        
        # Calculate percent dissociation
        percent_dissociation = (x / initial_concentration) * 100
        
        return {
            "formula": formula,
            "initial_concentration": initial_concentration,
            "ka": ka,
            "h_concentration": h_concentration,
            "oh_concentration": oh_concentration,
            "conjugate_base_concentration": a_concentration,
            "undissociated_concentration": ha_concentration,
            "ph": ph,
            "poh": poh,
            "percent_dissociation": percent_dissociation,
            "acid_type": "weak",
            "is_approximation_valid": (ka / initial_concentration < 0.05)
        }
    
    def solve_weak_base_equilibrium(self, formula, initial_concentration, kb=None):
        """
        Calculate equilibrium concentrations and pH for a weak base solution
        
        Args:
            formula (str): Chemical formula of the weak base
            initial_concentration (float): Initial concentration of the base in mol/L
            kb (float, optional): Base dissociation constant. If None, uses value from database.
            
        Returns:
            dict: Contains calculated equilibrium values
        """
        # Get Kb value if not provided
        if kb is None:
            if formula in self.weak_base_kb:
                kb = self.weak_base_kb[formula]
            else:
                return {
                    "error": f"Kb value for {formula} not found in database. Please provide Kb."
                }
        
        # Identify if it's a base
        compound_info = identify_acid_base(formula)
        if "Base" not in compound_info["classification"]:
            return {
                "error": f"{formula} is not identified as a base. Results may not be accurate.",
                "classification": compound_info["classification"]
            }
        
        # Solve the equilibrium problem
        # For a weak base B: B + H₂O ⇌ BH⁺ + OH⁻
        # Kb = [BH⁺][OH⁻]/[B]
        
        # Using the quadratic formula to solve for x (where x = [BH⁺] = [OH⁻])
        # x² + Kb*x - Kb*C₀ = 0
        # x = (-Kb + √(Kb² + 4*Kb*C₀))/2
        
        a = 1
        b = kb
        c = -kb * initial_concentration
        
        discriminant = b**2 - 4*a*c
        if discriminant < 0:
            return {"error": "Calculation error: negative discriminant"}
        
        x = (-b + math.sqrt(discriminant)) / (2*a)
        
        # Calculate concentrations at equilibrium
        oh_concentration = x  # [OH⁻]
        bh_concentration = x  # [BH⁺]
        b_concentration = initial_concentration - x  # [B]
        
        # Calculate pOH and pH
        poh = -math.log10(oh_concentration)
        ph = 14 - poh  # at 25°C, pH + pOH = 14
        h_concentration = 10**(-ph)
        
        # Calculate percent ionization
        percent_ionization = (x / initial_concentration) * 100
        
        return {
            "formula": formula,
            "initial_concentration": initial_concentration,
            "kb": kb,
            "oh_concentration": oh_concentration,
            "h_concentration": h_concentration,
            "conjugate_acid_concentration": bh_concentration,
            "undissociated_concentration": b_concentration,
            "ph": ph,
            "poh": poh,
            "percent_ionization": percent_ionization,
            "base_type": "weak",
            "is_approximation_valid": (kb / initial_concentration < 0.05)
        }
    
    def solve_ph_problem(self, formula, concentration, acid_base_type=None, ka=None, kb=None):
        """
        General pH solver that automatically determines if compound is strong/weak acid/base
        
        Args:
            formula (str): Chemical formula
            concentration (float): Concentration in mol/L
            acid_base_type (str, optional): Force specific calculation type
                                          Options: 'strong_acid', 'weak_acid', 'strong_base', 'weak_base'
            ka (float, optional): Acid dissociation constant for weak acids
            kb (float, optional): Base dissociation constant for weak bases
            
        Returns:
            dict: Complete pH calculation results
        """
        # Identify the compound
        compound_info = identify_acid_base(formula)
        
        # Lists of strong acids and bases
        strong_acids = ["HCl", "HBr", "HI", "HNO3", "H2SO4", "HClO4", "HClO3"]
        strong_bases = ["NaOH", "KOH", "LiOH", "RbOH", "CsOH", "Ca(OH)2", "Sr(OH)2", "Ba(OH)2"]
        
        # Determine calculation type if not specified
        if acid_base_type is None:
            if "Acid" in compound_info["classification"]:
                if formula in strong_acids:
                    acid_base_type = "strong_acid"
                else:
                    acid_base_type = "weak_acid"
            elif "Base" in compound_info["classification"]:
                if formula in strong_bases:
                    acid_base_type = "strong_base"
                else:
                    acid_base_type = "weak_base"
            else:
                return {
                    "error": f"Cannot determine if {formula} is an acid or base.",
                    "classification": compound_info["classification"]
                }
        
        # Perform appropriate calculation
        if acid_base_type == "strong_acid":
            return self.solve_strong_acid_equilibrium(formula, concentration)
        elif acid_base_type == "weak_acid":
            return self.solve_weak_acid_equilibrium(formula, concentration, ka)
        elif acid_base_type == "strong_base":
            return self.solve_strong_base_equilibrium(formula, concentration)
        elif acid_base_type == "weak_base":
            return self.solve_weak_base_equilibrium(formula, concentration, kb)
        else:
            return {"error": f"Invalid acid_base_type: {acid_base_type}"}
    
    def analyze_weak_acid_question(self, options, initial_concentration=1.0, formula="HX", ka=1e-5):
        """
        Analyzes a multiple choice question about weak acids
        
        Args:
            options (list): List of statement options to verify
            initial_concentration (float): Initial concentration of the acid
            formula (str): Chemical formula of the acid
            ka (float): Acid dissociation constant
            
        Returns:
            dict: Analysis of each option and the correct answer
        """
        # Get equilibrium data
        equilibrium = self.solve_weak_acid_equilibrium(formula, initial_concentration, ka)
        
        if "error" in equilibrium:
            return {"error": equilibrium["error"]}
        
        # Extract key values
        h_concentration = equilibrium["h_concentration"]
        x_concentration = equilibrium["conjugate_base_concentration"]
        hx_concentration = equilibrium["undissociated_concentration"]
        ph = equilibrium["ph"]
        
        # Analyze each option
        analysis = {}
        for i, option in enumerate(options):
            is_correct = False
            explanation = ""
            
            # Common option patterns
            if "pH = 0" in option:
                is_correct = abs(ph) < 0.01
                explanation = f"The pH is actually {ph:.4f}, not 0. This would only be true for a 1M strong acid."
                
            elif "[X-] = 1 M" in option or "[X⁻] = 1 M" in option:
                is_correct = abs(x_concentration - initial_concentration) < 0.0001
                explanation = f"[X⁻] = {x_concentration:.4e} M, which is less than 1 M because the weak acid only partially dissociates."
                
            elif "[HX] > [H+]" in option or "[HX] > [H⁺]" in option:
                is_correct = hx_concentration > h_concentration
                explanation = f"[HX] = {hx_concentration:.4e} M, which is {'greater' if is_correct else 'not greater'} than [H⁺] = {h_concentration:.4e} M."
                
            elif "[H+] = 1 M" in option or "[H⁺] = 1 M" in option:
                is_correct = abs(h_concentration - initial_concentration) < 0.0001
                explanation = f"[H⁺] = {h_concentration:.4e} M, which is much less than 1 M because the weak acid only partially dissociates."
                
            elif "both" in option.lower():
                explanation = "This option depends on which individual options are correct."
            
            analysis[f"option_{i+1}"] = {
                "statement": option,
                "is_correct": is_correct,
                "explanation": explanation
            }
        
        # Find the correct answer
        correct_options = [i+1 for i, data in enumerate(analysis.values()) if data["is_correct"]]
        
        return {
            "equilibrium_data": equilibrium,
            "options_analysis": analysis,
            "correct_options": correct_options,
            "summary": f"For a {initial_concentration} M solution of weak acid {formula} with Ka = {ka}, "
                     f"the equilibrium concentrations are: [H⁺] = {h_concentration:.4e} M, "
                     f"[X⁻] = {x_concentration:.4e} M, [HX] = {hx_concentration:.4e} M, pH = {ph:.4f}"
        }