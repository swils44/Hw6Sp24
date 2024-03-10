from Rankine_stem import rankine
from Steam_stem import steam
def analyze_rankine_cycle(p_high, p_low, t_high=None, name='Rankine Cycle'):
    cycle = rankine(p_low=p_low, p_high=p_high, t_high=t_high, name=name)
    efficiency = cycle.calc_efficiency()
    print(f"\n{name} Efficiency: {efficiency:.2f}%")
    cycle.print_summary()

def main():
    # Rankine cycle with saturated vapor entering the turbine
    print("Analyzing Rankine Cycle with Saturated Vapor:")
    analyze_rankine_cycle(p_high=8000, p_low=8, name='Rankine Cycle (Saturated Vapor)')

    # Rankine cycle with superheated steam entering the turbine
    # Assuming you have a way to calculate the saturation temperature at the given pressure (8000 kPa) using Steam.py
    saturated_steam= steam(pressure=8000, x=1)  # This is just an example; you would use your actual Steam class
    saturated_steam.calc()
    Tsat = saturated_steam.T
    T_superheated = 1.7 * Tsat

    print("\nAnalyzing Rankine Cycle with Superheated Steam:")
    analyze_rankine_cycle(p_high=8000, p_low=8, t_high=T_superheated, name='Rankine Cycle (Superheated Steam)')

if __name__ == "__main__":
    main()
