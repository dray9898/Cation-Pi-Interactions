import os
import subprocess

# Parameters
orca_executable = "/home/dr1015/orca_6_0_1/orca"  # Full path to ORCA executable
input_template = "test.inp"  # Template input file
output_folder = "orca_runs"            # Folder to store all output files
start_distance = 3.0                   # Starting z-offset for the group of 5 atoms
end_distance = 1.0                     # Ending z-offset for the group of 5 atoms
step_size = -0.25                       # Step size for z-offset
atom_indices = [14, 13, 16, 15, 17]    # Indices for the last 5 atoms (1-based)

# Function to modify z-coordinates for the group of atoms
def modify_geometry(template, z_offset, atom_indices):
    """
    Modify the z-coordinates of specified atoms in the input file.
    """
    with open(template, "r") as file:
        lines = file.readlines()

    # Locate the geometry block
    for i, line in enumerate(lines):
        if line.startswith("* xyz") or line.startswith("*"):
            geometry_start = i + 1
            break

    geometry = lines[geometry_start:]
    for atom_idx in atom_indices:
        atom_line = geometry[atom_idx - 1].split()  # Extract atom line
        atom_line[3] = str(float(atom_line[3]) + z_offset)  # Adjust z-coordinate
        geometry[atom_idx - 1] = "  ".join(atom_line) + "\n"  # Update geometry

    lines[geometry_start:] = geometry
    return lines

# Function to extract E(CCSD) energy from ORCA output
def extract_energy(output_file):
    """
    Extract the E(CCSD) energy from ORCA output file.
    """
    energy = None
    with open(output_file, "r") as file:
        for line in file:
            if "E(CCSD)" in line:  # Look for the E(CCSD) line
                energy = float(line.split()[-1])  # Extract the energy value
                break
    if energy is None:
        print(f"Warning: E(CCSD) energy not found in {output_file}.")
    return energy

# Main script
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

results = []
current_offset = start_distance

while current_offset >= end_distance:
    # Generate input file with modified geometry
    input_file = os.path.join(output_folder, f"input_{current_offset:.2f}.inp")
    modified_input = modify_geometry(input_template, current_offset - start_distance, atom_indices)
    with open(input_file, "w") as file:
        file.writelines(modified_input)

    # Run ORCA
    output_file = os.path.join(output_folder, f"output_{current_offset:.2f}.out")
    command = [orca_executable, input_file]
    print(f"Running ORCA for z-offset {current_offset:.2f} Å...")
    try:
        subprocess.run(command, stdout=open(output_file, "w"), stderr=subprocess.STDOUT, check=True)
        # Extract energy
        energy = extract_energy(output_file)
        if energy is not None:
            results.append((current_offset, energy))
            print(f"  Distance: {current_offset:.2f} Å, Energy: {energy:.6f} Eh")
        else:
            print(f"  Warning: No energy found for distance {current_offset:.2f} Å. Check output.")
    except subprocess.CalledProcessError as e:
        print(f"  Error: ORCA calculation failed for z-offset {current_offset:.2f} Å.")

    # Update z-offset for the next step
    current_offset += step_size

# Save results to a file
results_file = os.path.join(output_folder, "results.txt")
with open(results_file, "w") as file:
    file.write("Z-offset (Å)\tEnergy (Eh)\n")
    for z_offset, energy in results:
        file.write(f"{z_offset:.2f}\t{energy:.6f}\n")

print(f"Completed. Results saved in {results_file}.")
