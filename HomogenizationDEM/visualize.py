#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Rescaling factors (adjust as needed)
FORCE_VECTOR_SCALE = 50        # Larger value → shorter arrows.
PARTICLE_SIZE_FACTOR = 1e4       # Scales the area of particle circles.

##########################################
# File-loading Functions
##########################################

def load_dem_data(index):
    """Load DEM data from DEM_data/DEMdemo_output_<index>.csv."""
    filename = f"DEM_data/DEMdemo_output_{index}.csv"
    return pd.read_csv(filename)

def load_contact_data(index):
    """Load contact data from DEM_data/Contact_pairs_<index>.csv."""
    filename = f"DEM_data/Contact_pairs_{index}.csv"
    try:
        return pd.read_csv(filename)
    except Exception as e:
        print(f"Error loading contact file: {e}")
        return pd.DataFrame()  # Return an empty DataFrame if file not found

def load_output_data(index):
    """Load homogenization output data from output_data/output_data_<index>.csv."""
    filename = f"output_data/output_data_{index}.csv"
    return pd.read_csv(filename)

##########################################
# Net Force Computation
##########################################

def compute_particle_forces(dem_df, contact_df):
    """
    Compute the net contact force for each particle.
    The contact file is assumed to have columns:
      Type, a_idx, b_idx, f_x, f_y, f_z, x_c, y_c, z_c
    For each contact of type "SS", add the force (f_x, f_y, f_z) to particle a_idx
    and subtract it from particle b_idx.
    (Adjust indices here if your C code uses an offset.)
    """
    # Initialize a DataFrame of zeros for forces (one row per particle)
    forces = pd.DataFrame(0, index=dem_df.index, columns=["fx", "fy", "fz"], dtype=float)
    if contact_df.empty:
        return forces
    for _, row in contact_df.iterrows():
        if row["contact_type"].strip() != "SS":
            break
        try:
            # Adjust these indices if your C code uses an offset (e.g., subtract 4)
            a_idx = int(row["A"])
            b_idx = int(row["B"])
            fx = float(row["f_x"])
            fy = float(row["f_y"])
            fz = float(row["f_z"])
        except Exception as e:
            print(f"Error processing contact row: {e}")
            continue
        if a_idx in forces.index:
            forces.loc[a_idx, "fx"] += fx
            forces.loc[a_idx, "fy"] += fy
            forces.loc[a_idx, "fz"] += fz
        if b_idx in forces.index:
            forces.loc[b_idx, "fx"] -= fx
            forces.loc[b_idx, "fy"] -= fy
            forces.loc[b_idx, "fz"] -= fz
    return forces

##########################################
# Visualization Functions
##########################################

def plot_particles_forces_and_grid(dem_df, contact_df, output_df):
    """
    Plot the grid nodes from output_df, particle positions from dem_df,
    and overlay net contact force vectors computed from contact_df.
    """
    # Compute net forces from contact data
    net_forces = compute_particle_forces(dem_df, contact_df)
    dem_df = dem_df.copy()
    dem_df = dem_df.join(net_forces)
    
    plt.figure(figsize=(8,8))
    
    # 1) Plot grid nodes as small gray dots.
    plt.scatter(
        output_df['X'], output_df['Y'],
        s=10, c='gray', alpha=0.5, marker='.',
        label='Grid Nodes'
    )
    
    # 2) Plot particles as larger circles.
    plt.scatter(
        dem_df['X'], dem_df['Y'],
        s=dem_df['r'] * PARTICLE_SIZE_FACTOR,
        c='black', label='Particles'
    )
    
    # 3) Overlay net contact force vectors.
    plt.quiver(
        dem_df['X'], dem_df['Y'],
        dem_df['fx'], dem_df['fy'],
        color='green', angles='xy', scale_units='xy', scale=FORCE_VECTOR_SCALE,
        label='Net Contact Force'
    )
    
    # 4) Set axis limits to show the entire grid (based on output_df).
    x_min, x_max = output_df['X'].min(), output_df['X'].max()
    y_min, y_max = output_df['Y'].min(), output_df['Y'].max()
    margin = 0.01
    plt.xlim(x_min - margin, x_max + margin)
    plt.ylim(y_min - margin, y_max + margin)
    
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Particles, Net Contact Forces, and Grid Nodes")
    plt.legend()
    plt.axis('equal')
    plt.grid(True)
    plt.show()

def plot_density_field(output_df):
    """
    Plot a density field using the grid node positions (X, Y) and Density values.
    """
    pivot = output_df.pivot_table(index='Y', columns='X', values='Density')
    pivot = pivot.sort_index(ascending=True)
    plt.figure(figsize=(8,6))
    plt.imshow(
        pivot,
        extent=[pivot.columns.min(), pivot.columns.max(),
                pivot.index.min(), pivot.index.max()],
        origin='lower', aspect='auto', cmap='viridis'
    )
    plt.colorbar(label='Density')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Density Field")
    plt.show()

def plot_velocity_field(output_df):
    """
    Plot the homogenized velocity field as arrows using Vx and Vy.
    """
    plt.figure(figsize=(8,6))
    plt.quiver(
        output_df['X'], output_df['Y'],
        output_df['Vx'], output_df['Vy'],
        color='blue', angles='xy', scale_units='xy', scale=FORCE_VECTOR_SCALE
    )
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Homogenized Velocity Field")
    plt.axis('equal')
    plt.grid(True)
    plt.show()

def plot_stf_component(output_df, component):
    """
    Plot one component of the stress tensor field as a color-mapped image.
    """
    pivot = output_df.pivot_table(index='Y', columns='X', values=component)
    pivot = pivot.sort_index(ascending=True)
    plt.figure(figsize=(8,6))
    plt.imshow(
        pivot,
        extent=[pivot.columns.min(), pivot.columns.max(),
                pivot.index.min(), pivot.index.max()],
        origin='lower', aspect='auto', cmap='inferno'
    )
    plt.colorbar(label=component)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(f"{component} Field")
    plt.show()

##########################################
# Main Routine
##########################################

def main():
    # Usage: python visualize_forces.py <index> [<step>] [<STF_component>]
    if len(sys.argv) < 2:
        print("Usage: python visualize_forces.py <index> [<step>] [<STF_component>]")
        sys.exit(1)
    index = sys.argv[1]
    step = float(sys.argv[2]) if len(sys.argv) >= 3 else 0.2
    stf_component = sys.argv[3] if len(sys.argv) >= 4 else "STF_xx"

    print(f"Loading DEM data from DEM_data/DEMdemo_output_{index}.csv ...")
    dem_df = load_dem_data(index)
    print(f"Loading contact data from DEM_data/Contact_pairs_{index}.csv ...")
    contact_df = load_contact_data(index)
    print(f"Loading homogenization output from output_data/output_data_{index}.csv ...")
    output_df = load_output_data(index)
    
    print("Plotting particles, net contact forces, and grid nodes...")
    plot_particles_forces_and_grid(dem_df, contact_df, output_df)
    
    print("Plotting density field...")
    plot_density_field(output_df)
    
    print("Plotting homogenized velocity field...")
    plot_velocity_field(output_df)
    
    print(f"Plotting {stf_component} field...")
    plot_stf_component(output_df, component=stf_component)

if __name__ == '__main__':
    main()
