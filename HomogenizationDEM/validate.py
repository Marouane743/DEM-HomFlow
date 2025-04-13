#!/usr/bin/env python3
import csv
import math
import sys
import numpy as np
from scipy.integrate import simpson

# Tolerance for numerical comparisons
TOLERANCE = 1e-9

def almost_equal(a, b, tol=TOLERANCE):
    return abs(a - b) < tol

##########################################
# File-loading and Parsing Functions
##########################################

def load_header(header_filename):
    """
    Parse the header file.
    Expects lines like:
       @ Np
       2
       @ RadiusMax
       0.01
       @ Density
       2500
       @ Dimension
       2
       @ xsmin
       0.2
       @ xsmax
       0.8
       @ ysmin
       0.2
       @ ysmax
       0.8
       @ zsmin
       0.0
       @ zsmax
       0.0
    Returns a dictionary with keys:
       Np, RadiusMax, Density, Dimension, xsmin, xsmax, ysmin, ysmax, zsmin, zsmax.
    """
    header = {}
    with open(header_filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip() != ""]
    i = 0
    while i < len(lines):
        if lines[i].startswith("@"):
            parts = lines[i].split()
            if len(parts) >= 2:
                key = parts[1]
                if i+1 < len(lines):
                    value_line = lines[i+1]
                    try:
                        if '.' in value_line or 'e' in value_line.lower():
                            header[key] = float(value_line)
                        else:
                            header[key] = int(value_line)
                    except:
                        header[key] = value_line
                    i += 2
                    continue
        i += 1
    # For 2D, ensure zsmin and zsmax are set to 0 if missing.
    if header.get("Dimension", 2) == 2:
        header.setdefault("zsmin", 0.0)
        header.setdefault("zsmax", 0.0)
    return header

def load_dem(dem_filename):
    """
    Reads the DEM CSV file.
    Expects a header line like:
       x,y,z,r,vx,vy,vz   or   X,Y,Z,r,v_x,v_y,v_z
    Returns a list of particle dictionaries.
    """
    particles = []
    with open(dem_filename, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                # Try lower-case keys; if not present, try uppercase (and for velocities, also check for underscore)
                x_val = row.get("x", row.get("X"))
                y_val = row.get("y", row.get("Y"))
                z_val = row.get("z", row.get("Z"))
                r_val = row.get("r")  # assuming 'r' is the same in both cases
                vx_val = row.get("vx", row.get("v_x"))
                vy_val = row.get("vy", row.get("v_y"))
                vz_val = row.get("vz", row.get("v_z"))
                
                # Convert strings to floats
                particle = {
                    "x": float(x_val),
                    "y": float(y_val),
                    "z": float(z_val),
                    "r": float(r_val),
                    "vx": float(vx_val),
                    "vy": float(vy_val),
                    "vz": float(vz_val)
                }
            except Exception as e:
                print(f"Error converting DEM row: {e}")
                continue
            particles.append(particle)
    return particles


def read_output_csv(csv_filename):
    """
    Reads the output CSV (produced by your C project) and returns a list of dictionaries.
    """
    with open(csv_filename, newline='') as f:
        reader = csv.DictReader(f)
        return list(reader)

##########################################
# Core Mathematical Functions
##########################################

def compute_mass(r, density):
    return (4.0/3.0) * math.pi * (r**3) * density

def phi(q, s3):
    if 0 <= q <= 1:
        return s3 * (1 - 1.5 * q * q * (1 - q/2.0))
    elif 1 < q <= 2:
        return (s3 / 4.0) * ((2 - q) ** 3)
    else:
        return 0.0

def distance3d(p1, p2):
    return math.sqrt((p1['x'] - p2['x'])**2 +
                     (p1['y'] - p2['y'])**2 +
                     (p1['z'] - p2['z'])**2)

##########################################
# Integration for Collision Contributions
##########################################

def norm_func(s, x, dj, theta, RC):
    term1 = (1 - s)**2 * x**2
    term2 = 2 * (1 - s) * s * dj * math.cos(theta) * x
    term3 = s**2 * dj**2
    inside = term1 + term2 + term3
    return math.sqrt(inside) / RC

def integrand(s, x, dj, theta, s3, RC):
    q = 2.0 * norm_func(s, x, dj, theta, RC)
    return phi(q, s3)

def integral_phi(x, dj, theta, s3, RC, n_points=101):
    s_vals = np.linspace(0, 1, n_points)
    vals = [integrand(s, x, dj, theta, s3, RC) for s in s_vals]
    return simpson(vals, s_vals)

##########################################
# Expected Value Computation (Kinetic Part)
##########################################

def compute_expected(header, particles, RC, step):
    """
    Re‑implements the non‑collision part of processSingleSet():
      - Builds a grid with the given step.
      - For each particle (inside the grid), computes contributions to density,
        momentum, and the kinetic part of the stress tensor.
    Returns a list of dictionaries (one per grid node) with keys:
      NodeID, X, Y, Z, Density, Vx, Vy, Vz, STF_xx, STF_xy, STF_xz,
      STF_yx, STF_yy, STF_yz, STF_zx, STF_zy, STF_zz.
    Also returns grid parameters (xgmin, ygmin, step, Cx, Cy, Cz) for later use.
    """
    RadiusMax = header["RadiusMax"]
    xsmin = header["xsmin"]
    xsmax = header["xsmax"]
    ysmin = header["ysmin"]
    ysmax = header["ysmax"]
    D = header["Dimension"]
    if D == 2:
        zgmin = 0.0
        zgmax = 0.0
    else:
        zgmin = 2 * RC
        zgmax = header["zsmax"] - 2 * RC

    xgmin = xsmin + 2 * RC
    xgmax = xsmax - 2 * RC
    ygmin = ysmin + 2 * RC
    ygmax = ysmax - 2 * RC

    Cx = int((xgmax - xgmin) / step)
    Cy = int((ygmax - ygmin) / step)
    Cz = 0 if D == 2 else int((zgmax - zgmin) / step)
    total_nodes = (Cx + 1) * (Cy + 1) * (Cz + 1)

    # Build grid nodes (ordered as in C: for k, for j, for i)
    grid_nodes = []
    for k in range(Cz+1):
        for j in range(Cy+1):
            for i in range(Cx+1):
                node = {"x": xgmin + i * step,
                        "y": ygmin + j * step,
                        "z": zgmin + k * step}
                grid_nodes.append(node)

    density_field = [0.0] * total_nodes
    momentum = [{"x": 0.0, "y": 0.0, "z": 0.0} for _ in range(total_nodes)]
    stf = [{"xx": 0.0, "xy": 0.0, "xz": 0.0,
            "yx": 0.0, "yy": 0.0, "yz": 0.0,
            "zx": 0.0, "zy": 0.0, "zz": 0.0} for _ in range(total_nodes)]

    s3 = 1.0 / (math.pi * ((RC/2.0) ** 3))
    # Loop over all particles and accumulate contributions
    for particle in particles:
        # Skip particles outside grid bounding box
        if (particle["x"] < xgmin or particle["x"] > xgmax or
            particle["y"] < ygmin or particle["y"] > ygmax or
            particle["z"] < zgmin or particle["z"] > zgmax):
            continue
        mass = compute_mass(particle["r"], header["Density"])
        w = int(RC / step)
        ix = int((particle["x"] - xgmin) / step)
        iy = int((particle["y"] - ygmin) / step)
        iz = int((particle["z"] - zgmin) / step) if D == 3 else 0
        nxmin = max(ix - w, 0)
        nymin = max(iy - w, 0)
        nzmin = max(iz - w, 0)
        nxmax = min(ix + w + 1, Cx)
        nymax = min(iy + w + 1, Cy)
        nzmax = min(iz + w + 1, Cz) if D == 3 else 0

        for k in range(nzmin, (nzmax+1) if D==3 else 1):
            for j in range(nymin, nymax+1):
                for i in range(nxmin, nxmax+1):
                    idx = i + j*(Cx+1) + (k*(Cx+1)*(Cy+1) if D==3 else 0)
                    node = grid_nodes[idx]
                    d = math.sqrt((node["x"] - particle["x"])**2 +
                                  (node["y"] - particle["y"])**2 +
                                  (node["z"] - particle["z"])**2)
                    q = 2.0 * d / RC
                    if q >= 2.0:
                        continue
                    ph_val = phi(q, s3)
                    density_field[idx] += mass * ph_val
                    momentum[idx]["x"] += mass * particle["vx"] * ph_val
                    momentum[idx]["y"] += mass * particle["vy"] * ph_val
                    momentum[idx]["z"] += mass * particle["vz"] * ph_val
                    stf[idx]["xx"] -= mass * (particle["vx"]**2) * ph_val
                    stf[idx]["xy"] -= mass * (particle["vx"] * particle["vy"]) * ph_val
                    stf[idx]["xz"] -= mass * (particle["vx"] * particle["vz"]) * ph_val
                    stf[idx]["yy"] -= mass * (particle["vy"]**2) * ph_val
                    stf[idx]["yz"] -= mass * (particle["vy"] * particle["vz"]) * ph_val
                    stf[idx]["zz"] -= mass * (particle["vz"]**2) * ph_val
                    stf[idx]["yx"] = stf[idx]["xy"]
                    stf[idx]["zx"] = stf[idx]["xz"]
                    stf[idx]["zy"] = stf[idx]["yz"]

    # Compute velocity field and subtract kinetic part from STF.
    velocity_field = [{"x": 0.0, "y": 0.0, "z": 0.0} for _ in range(total_nodes)]
    EPSILON = 1e-8
    for idx in range(total_nodes):
        if density_field[idx] > EPSILON:
            inv_rho = 1.0 / density_field[idx]
            velocity_field[idx]["x"] = momentum[idx]["x"] * inv_rho
            velocity_field[idx]["y"] = momentum[idx]["y"] * inv_rho
            velocity_field[idx]["z"] = momentum[idx]["z"] * inv_rho
            stf[idx]["xx"] += density_field[idx] * (velocity_field[idx]["x"]**2)
            stf[idx]["xy"] += density_field[idx] * (velocity_field[idx]["x"] * velocity_field[idx]["y"])
            stf[idx]["xz"] += density_field[idx] * (velocity_field[idx]["x"] * velocity_field[idx]["z"])
            stf[idx]["yy"] += density_field[idx] * (velocity_field[idx]["y"]**2)
            stf[idx]["yz"] += density_field[idx] * (velocity_field[idx]["y"] * velocity_field[idx]["z"])
            stf[idx]["zz"] += density_field[idx] * (velocity_field[idx]["z"]**2)
            stf[idx]["yx"] = stf[idx]["xy"]
            stf[idx]["zx"] = stf[idx]["xz"]
            stf[idx]["zy"] = stf[idx]["yz"]

    results = []
    for idx, node in enumerate(grid_nodes):
        result = {
            "NodeID": idx,
            "X": node["x"],
            "Y": node["y"],
            "Z": node["z"],
            "Density": density_field[idx],
            "Vx": velocity_field[idx]["x"],
            "Vy": velocity_field[idx]["y"],
            "Vz": velocity_field[idx]["z"],
            "STF_xx": stf[idx]["xx"],
            "STF_xy": stf[idx]["xy"],
            "STF_xz": stf[idx]["xz"],
            "STF_yx": stf[idx]["yx"],
            "STF_yy": stf[idx]["yy"],
            "STF_yz": stf[idx]["yz"],
            "STF_zx": stf[idx]["zx"],
            "STF_zy": stf[idx]["zy"],
            "STF_zz": stf[idx]["zz"]
        }
        results.append(result)
    grid_params = {"xgmin": xgmin, "ygmin": ygmin, "step": step, "Cx": Cx, "Cy": Cy, "Cz": Cz}
    return results, grid_params

##########################################
# Collision (Contact) Contributions
##########################################

def apply_collision_contributions(index_str, results, grid_params, header, RC, particles):
    """
    Reads the contact file from DEM_data/Contact_pairs_<index_str>.csv.
    For each contact of type "SS", computes the collision integral using SciPy's Simpson integration
    and updates the stress tensor (STF) in the computed results.
    No extra scaling is applied.
    """
    step = grid_params["step"]
    RadiusMax = header["RadiusMax"]
    xsmin = header["xsmin"]
    ysmin = header["ysmin"]
    D = header["Dimension"]
    if D == 2:
        zgmin = 0.0
    else:
        zgmin = 4 * RadiusMax
    xgmin = xsmin + 4 * RadiusMax
    ygmin = ysmin + 4 * RadiusMax
    Cx = grid_params["Cx"]
    Cy = grid_params["Cy"]
    Cz = grid_params["Cz"]

    s3 = 1.0/(math.pi*((RC/2.0)**3))
    contact_filename = f"DEM_data/Contact_pairs_{index_str}.csv"
    try:
        with open(contact_filename, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Contact file {contact_filename} not found; skipping collision contributions.")
        return results

    if len(lines) < 2:
        print(f"Contact file {contact_filename} is empty; no collision contributions.")
        return results

    # Process each contact line (skip header)
    for line in lines[1:]:
        parts = line.strip().split(',')
        if len(parts) < 9:
            continue
        contact_type = parts[0].strip()
        if contact_type != "SS":
            break
        try:
            a_idx = int(parts[1])
            b_idx = int(parts[2])
            F = {"x": float(parts[3]), "y": float(parts[4]), "z": float(parts[5])}
        except:
            continue
        # Adjust indices as in C: A = a_idx - 4, B = b_idx - 4.
        A = a_idx - 4
        B = b_idx - 4
        if A < 0 or A >= len(particles) or B < 0 or B >= len(particles):
            continue
        rpA = {"x": particles[A]["x"], "y": particles[A]["y"], "z": particles[A]["z"]}
        rpB = {"x": particles[B]["x"], "y": particles[B]["y"], "z": particles[B]["z"]}
        if (rpA["x"]==0 and rpA["y"]==0 and rpA["z"]==0 and 
            rpB["x"]==0 and rpB["y"]==0 and rpB["z"]==0):
            continue
        rpAB = {"x": rpA["x"] - rpB["x"],
                "y": rpA["y"] - rpB["y"],
                "z": rpA["z"] - rpB["z"]}
        x_min = min(rpA["x"], rpB["x"]) - RC
        x_max = max(rpA["x"], rpB["x"]) + RC
        y_min = min(rpA["y"], rpB["y"]) - RC
        y_max = max(rpA["y"], rpB["y"]) + RC
        z_min = min(rpA["z"], rpB["z"]) - RC
        z_max = max(rpA["z"], rpB["z"]) + RC

        nxmin = max(int((x_min - xgmin)/step), 0)
        nymin = max(int((y_min - ygmin)/step), 0)
        nzmin = 0 if D==2 else max(int((z_min - zgmin)/step), 0)
        nxmax = min(int((x_max - xgmin)/step), Cx)
        nymax = min(int((y_max - ygmin)/step), Cy)
        nzmax = 0 if D==2 else min(int((z_max - zgmin)/step), Cz)

        for k in range(nzmin, (nzmax+1) if D==3 else 1):
            for j in range(nymin, nymax+1):
                for i in range(nxmin, nxmax+1):
                    idx = i + j*(Cx+1) + (k*(Cx+1)*(Cy+1) if D==3 else 0)
                    node = {"x": xgmin + i*step,
                            "y": ygmin + j*step,
                            "z": zgmin + (k*step if D==3 else 0.0)}
                    xA = math.sqrt((node["x"] - rpA["x"])**2 +
                                   (node["y"] - rpA["y"])**2 +
                                   (node["z"] - rpA["z"])**2)
                    xB = math.sqrt((node["x"] - rpB["x"])**2 +
                                   (node["y"] - rpB["y"])**2 +
                                   (node["z"] - rpB["z"])**2)
                    vA = (rpA["x"] - node["x"], rpA["y"] - node["y"], rpA["z"] - node["z"])
                    vB = (rpB["x"] - node["x"], rpB["y"] - node["y"], rpB["z"] - node["z"])
                    lenA = xA
                    lenB = xB
                    dotAB = vA[0]*vB[0] + vA[1]*vB[1] + vA[2]*vB[2]
                    cosT = 0.0
                    if lenA > 1e-8 and lenB > 1e-8:
                        cosT = dotAB / (lenA * lenB)
                    cosT = max(min(cosT, 1.0), -1.0)
                    theta = math.acos(cosT)
                    integral_value = integral_phi(xA, xB, theta, s3, RC)
                    if abs(integral_value) < 1e-12:
                        continue
                    results[idx]["STF_xx"] -= F["x"] * rpAB["x"] * integral_value
                    results[idx]["STF_xy"] -= F["x"] * rpAB["y"] * integral_value
                    results[idx]["STF_xz"] -= F["x"] * rpAB["z"] * integral_value
                    results[idx]["STF_yx"] -= F["y"] * rpAB["x"] * integral_value
                    results[idx]["STF_yy"] -= F["y"] * rpAB["y"] * integral_value
                    results[idx]["STF_yz"] -= F["y"] * rpAB["z"] * integral_value
                    results[idx]["STF_zx"] -= F["z"] * rpAB["x"] * integral_value
                    results[idx]["STF_zy"] -= F["z"] * rpAB["y"] * integral_value
                    results[idx]["STF_zz"] -= F["z"] * rpAB["z"] * integral_value
    return results

##########################################
# Comparison Function
##########################################

def compare_results(computed, output):
    """
    Compare each numeric field of the computed results with the output CSV.
    Returns a list of error messages.
    """
    errors = []
    if len(computed) != len(output):
        errors.append(f"Number of nodes mismatch: computed {len(computed)} vs output {len(output)}")
    n = min(len(computed), len(output))
    for i in range(n):
        comp = computed[i]
        out = output[i]
        for key in comp:
            try:
                if key == "NodeID":
                    comp_val = int(comp[key])
                    out_val = int(float(out[key]))
                else:
                    comp_val = float(comp[key])
                    out_val = float(out[key])
            except Exception as e:
                errors.append(f"Row {i} field {key}: conversion error ({e})")
                continue
            if not almost_equal(comp_val, out_val):
                errors.append(f"Row {i} field {key}: computed {comp_val} vs output {out_val}")
    return errors

##########################################
# Main Routine
##########################################

def main():
    # Usage: python validate.py <RC> [<index>] [<step>]
    if len(sys.argv) < 2:
        print("Usage: python validate.py <RC> [<index>] [<step>]")
        sys.exit(1)
    try:
        RC = float(sys.argv[1])
    except:
        print("RC must be a positive number")
        sys.exit(1)
    index_str = sys.argv[2] if len(sys.argv) >= 3 else "0000"
    step = float(sys.argv[3]) if len(sys.argv) >= 4 else 0.2

    # Build file paths automatically
    header_file = f"DEM_data/header{index_str}.txt"
    dem_file = f"DEM_data/DEMdemo_output_{index_str}.csv"
    output_csv = f"output_data/output_data_{index_str}.csv"

    print(f"Using header file: {header_file}")
    print(f"Using DEM file: {dem_file}")
    print(f"Using output CSV file: {output_csv}")
    print(f"Using grid step: {step}")

    header = load_header(header_file)
    particles = load_dem(dem_file)
    if not particles:
        print("No particle data found in DEM file.")
        sys.exit(1)

    computed_results, grid_params = compute_expected(header, particles, RC, step)
    computed_results = apply_collision_contributions(index_str, computed_results, grid_params, header, RC, particles)
    output_results = read_output_csv(output_csv)
    errors = compare_results(computed_results, output_results)
    if errors:
        print("Validation FAILED with the following errors:")
        for err in errors:
            print("  " + err)
    else:
        print("All tests passed!")

if __name__ == "__main__":
    main()
