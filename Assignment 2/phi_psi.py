import numpy as np
import matplotlib.pyplot as plt
import math

def read_pdb(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    atom_lines = [line for line in lines if line.startswith("ATOM")]
    return atom_lines

def parse_atoms(atom_lines):
    atoms = []
    for line in atom_lines:
        atom_info = {
            'atom_number': int(line[6:11].strip()),
            'atom_name': line[12:16].strip(),
            'residue_name': line[17:20].strip(),
            'chain': line[21],
            'residue_number': int(line[22:26].strip()),
            'x': float(line[30:38].strip()),
            'y': float(line[38:46].strip()),
            'z': float(line[46:54].strip())
        }
        atoms.append(atom_info)
    return atoms

def group_atoms_by_residue(atoms):
    residues = {}
    for atom in atoms:
        key = (atom['chain'], atom['residue_number'])
        if key not in residues:
            residues[key] = []
        residues[key].append(atom)
    return residues

def calculate_dihedral(p1, p2, p3, p4):
    p1, p2, p3, p4 = np.array(p1), np.array(p2), np.array(p3), np.array(p4)
    v1 = p2 - p1
    v2 = p3 - p2
    v3 = p4 - p3
    n1 = np.cross(v1, v2)
    n2 = np.cross(v2, v3)
    x = np.dot(n1, n2) # essentially gives the cos(theta) or the angle, but
    y = np.dot(np.cross(n1, v2), n2) # used this to find the sign/ direction of the angle
    angle = -math.atan2(y, x)
    return round(np.degrees(angle), 3)

def calculate_phi_psi(residues):
    phi_psi_angles = []
    for i, ((chain1, res_num1), atoms1) in enumerate(residues.items()):
        if i == 0 or i == len(residues) - 1:
            continue

        (chain0, res_num0), atoms0 = list(residues.items())[i-1]
        (chain2, res_num2), atoms2 = list(residues.items())[i+1]

        if chain0 != chain1 or chain1 != chain2:
            continue

        c0_list  = [a for a in atoms0 if a['atom_name'] == 'C']
        n1_list  = [a for a in atoms1 if a['atom_name'] == 'N']
        ca1_list = [a for a in atoms1 if a['atom_name'] == 'CA']
        c1_list  = [a for a in atoms1 if a['atom_name'] == 'C']
        n2_list  = [a for a in atoms2 if a['atom_name'] == 'N']
        if len(c0_list) and len(n1_list) and len(ca1_list) and len(c1_list) and len(n2_list):
            c0  = c0_list[0]
            n1  = n1_list[0]
            ca1 = ca1_list[0]
            c1  = c1_list[0]
            n2  = n2_list[0]
        else:
            continue  # Skip this iteration and proceed to the next
        
        phi = calculate_dihedral([c0['x'], c0['y'], c0['z']],
                                 [n1['x'], n1['y'], n1['z']],
                                 [ca1['x'], ca1['y'], ca1['z']],
                                 [c1['x'], c1['y'], c1['z']])
        
        psi = calculate_dihedral([n1['x'], n1['y'], n1['z']],
                                 [ca1['x'], ca1['y'], ca1['z']],
                                 [c1['x'], c1['y'], c1['z']],
                                 [n2['x'], n2['y'], n2['z']])
        
        phi_psi_angles.append((phi, psi))
    return phi_psi_angles

def plot_ramachandran(angles):
    phis, psis = zip(*angles)
    plt.scatter(phis, psis)
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.xlabel("Phi angles (degrees)")
    plt.ylabel("Psi angles (degrees)")
    plt.title("Ramachandran Plot")
    plt.grid(True)
    plt.savefig("Ramachandran_Plot.png")
    plt.show()

atom_lines = read_pdb("1ASY.PDB")
atoms = parse_atoms(atom_lines)
residues = group_atoms_by_residue(atoms)
phi_psi_angles = calculate_phi_psi(residues)
plot_ramachandran(phi_psi_angles)

# Calculate Phi and Psi for the 5th residue
phi_5th, psi_5th = phi_psi_angles[4]
print(f"Phi and Psi angles for the 5th residue are {phi_5th} and {psi_5th}, respectively.")
with open('output.txt', 'w') as file:
    file.write(f"Phi and Psi angles for the 5th residue are {phi_5th} and {psi_5th}, respectively.\n")
