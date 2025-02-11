import pandas as pd
import magnolia.bondfile_parser as bfp
import networkx as nx
import multiprocessing as mp

# Define your directory and other parameters
dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\11_ARLOI_Slow_Ramp_1KpPS_300-500K\A0001\Production\300-500K_TRR=1Kpps"
H, C, O = 1, 2, 3  # Atom types
Lmin, Lmax = 30, 30  # Chain length limits
timestep, N_PAOr = 0.25, 20

# Function to process each simulation
def process_sim(sim):
    bondfile = rf"{dirr}\Sim-{1+sim}\bonds.out"
    atominfo = bfp.parsebondfile(bondfile, mtypes=True)
    
    neighbours, atypes, mtypes = atominfo["neighbours"], atominfo["atypes"], atominfo["mtypes"]
    SR = {}  # Count of scavenged radicals
    
    for step, neigh in neighbours.items():
        molecules = bfp.get_molecules(neigh)
        G = nx.Graph(neigh)
        time = step * timestep / 1000  # Convert to ps
        SR[time] = 0
        for molecule in molecules:
            g = G.subgraph(molecule)
            L = [atypes[x] for x in molecule].count(C)  # Count carbons only
            if L < Lmin or L > Lmax:
                continue
            for parent in g.nodes():
                c_types = sorted([atypes[x] for x in g[parent]])
                if atypes[parent] == O and mtypes[parent] <= N_PAOr and c_types == [H, C]:
                    SR[time] += 1
    return f'Sim-{sim+1}', SR

# Multi-processing to process simulations in parallel
if __name__ == "__main__":
    with mp.Pool(processes=mp.cpu_count()) as pool:  # Use all available cores
        results = pool.map(process_sim, range(3))  # Run 3 simulations in parallel

    # Convert results to DataFrame and save
    result_dict = dict(results)
    df = pd.DataFrame(result_dict).dropna(how='any')
    df.to_csv(dirr + "\\number_of_SR_vs_time.csv")
