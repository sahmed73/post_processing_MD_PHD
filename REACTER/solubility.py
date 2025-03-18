import magnolia.log_parser_SA as lfp
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('default')
plt.rcParams['font.size'] = 18

fig, ax = plt.subplots(dpi=300)
E_values = []
n_bulk = 100  # Number of molecules in the bulk simulation
molecule = 'Chloroethane'

for mode in ['Gas', 'Bulk']:
    print('-' * 50)
    print(mode)
    
    sim = 'Sim-1'
    # sim = 'Sim-1' if mode == 'Gas' else 'Sim-1_NVT-NPT'
    
    dirr = rf"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\Solubility\Chloroethane_OLD_FINE\PCFF-Gastiger\Eq\{mode}_phase\{sim}"
    filename = r'\log.lammps'
    
    logfile = dirr + filename
    thermo = lfp.thermo_panda(logfile, serial=-1, zero_ref='Time')
    
    time, pe = thermo['Time'], thermo['PotEng']
    start = -2000
    end = None

    if mode == 'Bulk':
        pe = pe / n_bulk
        V_bulk = thermo['Volume'].iloc[start:end].mean()  # Average last 1000 values
        eq_density = thermo['Density'].iloc[start:end].mean()  # Average density
        
    time, pe = time.iloc[start:end] - time.iloc[start:end].min(), pe.iloc[start:end]
    E_values.append(pe.mean())
    
    ax.plot(time, [pe.mean()] * time.size, color='k')
    ax.scatter(time, pe, s=10, label=mode)
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('P.E. per molecule (Kcal/mol)')

ax.legend(fontsize=16)

E_gas, E_bulk = E_values  # kcal/mol
V_m = (V_bulk / n_bulk) * (6.022e23 * 1e-24)  # cm³ per mole

conversion_factor = 64.6838  # (MPa^0.5) / (kcal/mol)^0.5 / cm^(3/2)
delta_E_vap = (E_gas - E_bulk) 
delta = np.sqrt(delta_E_vap / V_m) * conversion_factor

print('-' * 50)
print(f"E_gas   = {E_gas:.4f} kcal/mol")
print(f"E_bulk  = {E_bulk:.4f} kcal/mol")
print(f"E_vap   = {E_gas - E_bulk:.4f} kcal/mol")

print(f"V_bulk  = {V_bulk:.4f} Å³")
print(f"V_m     = {V_m:.4f} cm³/mol")
print(f"Density = {eq_density:.2f} g/cc")

print(f"Solubility Parameter (δ): {delta:.2f} MPa^0.5")
