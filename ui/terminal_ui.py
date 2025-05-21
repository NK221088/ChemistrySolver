"""
Common Terminal UI Functions for Chemistry Problem Solver
"""
import os

def clear_screen():
    """Clear the terminal screen."""
    os.system('cls' if os.name == 'nt' else 'clear')

def display_title(title):
    """Display a centered title for a module."""
    clear_screen()
    print("\n" + "="*50)
    print(title.center(50))
    print("="*50 + "\n")

def display_main_menu():
    """Display the main menu options."""
    print("\nPlease select a module:")
    print("  [1] Electron Configuration")
    print("  [2] Chemical Kinetics")
    print("  [3] Molar Mass Calculator")
    print("  [4] Chemical Equation Balancer")
    print("  [5] Thermodynamics Calculator")
    print("  [6] Acid-Base Analyzer")
    print("  [7] Name to Formula Converter")  # New option
    print("  [0] Exit")

def get_user_choice():
    """Get user choice from the main menu."""
    return input("\nEnter your choice: ").strip()

def display_results_header():
    """Display a header for results section."""
    print("\n" + "-"*50)
    print("RESULTS".center(50))
    print("-"*50)

def display_steps(steps):
    """Display solution steps."""
    print("\n" + "-"*50)
    print("SOLUTION STEPS".center(50))
    print("-"*50)
    
    for step in steps:
        print(step)

def wait_for_user():
    """Wait for user before continuing."""
    input("\nPress Enter to continue...")
