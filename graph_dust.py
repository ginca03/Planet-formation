import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm  # For logscale in density
import argparse  #for getting input from bash

def plot_density_evolution(filename):
 
    # open file root
    with uproot.open(filename) as file:
    
        tree_name = list(file.keys())[0]
        tree = file[tree_name]

        t = tree["t"].array(library="np")
        z = tree["z"].array(library="np")
        rho = tree["rho"].array(library="np")

 
    t_unique = np.unique(t)
    z_unique = np.unique(z)

  
    matrix = rho.reshape(len(t_unique), len(z_unique))

    # Plot
    plt.figure(figsize=(10, 6))
    plt.imshow(
        matrix,
        aspect="auto",
        origin="lower",
        extent=[z.min(), z.max(), t.min(), t.max()],
        cmap="viridis",
        norm=LogNorm(vmin=matrix[matrix > 0].min(), vmax=matrix.max())
    )
    plt.colorbar(label="Density (⍴_0)")
    plt.xlabel("z [H(R)]")
    plt.ylabel("t (years)")
    plt.title("Dust Density Evolution (Log Scale)")

    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot dust evolution from file .root")
    parser.add_argument("filename", type=str, help="name of .root file to plot")
    args = parser.parse_args()

    plot_density_evolution(args.filename)
