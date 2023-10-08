import numpy as np
# ======================= DEFINING FUNCTIONS (BEGIN) ============================= #
def read_asa_file(filename):
    asa_data = []
    protein_asa_data = []
    rna_asa_data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                parts = line.split()
                # asa_value = float(line[54:62].strip())
                asa_value = float(parts[-2])
                atom_info = {
                    'line': line.strip(),
                    'atom_number': int(line[6:11].strip()),
                    'atom_type': line[12:16].strip(),
                    'class' : 'protein' if len(line[17:20].strip()) == 3 else 'rna',
                    'chain': line[21],
                    'asa': asa_value
                }
                if atom_info['class'] == 'protein':
                    protein_asa_data.append(atom_info)
                else:
                    rna_asa_data.append(atom_info)
                asa_data.append(atom_info)
    return asa_data, protein_asa_data, rna_asa_data

def write_int_file(filename, asa_bound, asa_unbound):
    with open(filename, 'w') as f:
        for bound, unbound in zip(asa_bound, asa_unbound):
            diff = bound['asa'] - unbound['asa']
            f.write(f"{bound['line']}\t{bound['asa']}\t{unbound['asa']}\t{diff}")

def calculate_sasa(asa_data):
    return np.sum([atom['asa'] for atom in asa_data])

def filter_atoms(asa_data, atom_types):
    return [atom for atom in asa_data if atom['atom_type'] in atom_types]

def interface_atoms(asa_bound, asa_unbound):
    atoms = [{**b, 'difference_in_asa' : u['asa'] - b['asa']} for (b, u) in zip(asa_bound, asa_unbound) if u['asa'] != b['asa']]
    return atoms

def count_boundary(interface_asa_data, threshold):
    boundary_atoms = [atom for atom in interface_asa_data if atom['difference_in_asa'] < threshold and atom['difference_in_asa'] > 0.0]
    return len(boundary_atoms)

def count_core(asa_data_protein, asa_data_rna, threshold):
    core_atoms = [atom for atom in asa_data_protein + asa_data_rna if abs(atom['asa']) == threshold]
    return len(core_atoms)

def count_surface_atoms(asa_data_protein, asa_data_rna, threshold):
    return len([atom for atom in asa_data_protein + asa_data_rna if atom['asa'] > threshold])

def percentage(numerator, denominator):
    return 100 * numerator / denominator if denominator else 0

# ======================= DEFINING FUNCTIONS (END) =============================== #


# ========================= IMPLEMENTATION (BEGIN) =============================== #

# configuring variables to segeregate and store information required to compute SASA etc
asa_complex, asa_complex_protein, asa_complex_rna = read_asa_file('1ASY_C.asa')
___________,    asa_protein     , _______________ = read_asa_file('1ASY_P.asa')
___________, ___________________,     asa_rna     = read_asa_file('1ASY_R.asa')
# backbone atoms - N, C, alpha C
backbone_atoms_protein = filter_atoms(asa_protein, ['N', 'C', 'CA'])
sidechain_atoms_protein = filter_atoms(asa_protein, set([atom['atom_type'] for atom in asa_protein]) - {'N', 'C', 'CA'}) 
# (set data structure stores only distinct elements, atoms in this case)


# SASA Values
total_sasa_complex      = round(calculate_sasa(asa_complex), 3)
total_sasa_protein      = round(calculate_sasa(asa_protein), 3)
total_sasa_rna          = round(calculate_sasa(asa_rna), 3)
total_sasa_of_interface = round(total_sasa_protein + total_sasa_rna - total_sasa_complex, 3)
total_sasa_backbone     = round(calculate_sasa(backbone_atoms_protein), 3)
total_sasa_sidechain    = round(calculate_sasa(sidechain_atoms_protein), 3)

# Counting number of atoms
protein_atoms_on_interface     = interface_atoms(asa_complex_protein, asa_protein)
rna_atoms_on_interface         = interface_atoms(asa_complex_rna, asa_rna)
atoms_on_interface             = protein_atoms_on_interface + rna_atoms_on_interface
num_atoms_interface            = len(atoms_on_interface)
boundary_atoms_count           = count_boundary(atoms_on_interface, 1.0)
core_atoms_count               = count_core(asa_protein, asa_rna, 0.0)
# core_atoms_count               = count_core(asa_complex, 0.0)
surface_atoms_count            = count_surface_atoms(asa_protein, asa_rna, 0.5)
# surface_atoms_count            = count_surface_atoms(asa_complex, 0.5)

percent_sasa_sidechain_overall = round(percentage(total_sasa_sidechain, total_sasa_protein), 3)
percent_sasa_sidechain_protein = round(percentage(total_sasa_sidechain, calculate_sasa(protein_atoms_on_interface)), 3)
# ========================== IMPLEMENTATION (END) ================================ #


# ============================= OUTPUT (BEGIN) =================================== #
print(f"number of sidechain of protein atoms: {len(sidechain_atoms_protein)}")
print(f"number of protein atoms on interface: {len(protein_atoms_on_interface)}")
print(f"Total SASA of the interface in the complex: {total_sasa_of_interface}")
print(f"Total SASA of the protein: {total_sasa_protein}")
print(f"Total SASA of the RNA: {total_sasa_rna}")
print(f"Total SASA of backbone atoms of protein: {total_sasa_backbone}")
print(f"Total SASA of side chain atoms of the protein: {total_sasa_sidechain}")
print(f"Number of atoms in the interface: {num_atoms_interface}")
print(f"Number of atoms in the boundary of the interface: {boundary_atoms_count}")
print(f"Number of atoms in the core: {core_atoms_count}")
print(f"Number of atoms on the surface: {surface_atoms_count}")
print(f"% Total SASA contribution by side chain protein atoms w.r.t all protein atoms: {percent_sasa_sidechain_overall}")
print(f"% Total SASA contribution by side chain protein atoms w.r.t interface protein atoms: {percent_sasa_sidechain_protein}")

write_int_file('1ASY_P.INT', asa_complex, asa_protein)
write_int_file('1ASY_R.INT', asa_complex, asa_rna)

with open('output.txt', 'w') as f:
    f.write(f"Total SASA of the interface in the complex: {total_sasa_of_interface}\n")
    f.write(f"Total SASA of the protein: {total_sasa_protein}\n")
    f.write(f"Total SASA of the RNA: {total_sasa_rna}\n")
    f.write(f"Total SASA of backbone atoms of protein: {total_sasa_backbone}\n")
    f.write(f"Total SASA of side chain atoms of the protein: {total_sasa_sidechain}\n")
    f.write(f"Number of atoms in the interface: {num_atoms_interface}\n")
    f.write(f"Number of atoms in the boundary of the interface: {boundary_atoms_count}\n")
    f.write(f"Number of atoms in the core: {core_atoms_count}\n")
    f.write(f"Number of atoms on the surface: {surface_atoms_count}\n")
    f.write(f"% Total SASA contribution by side chain protein atoms overall: {percent_sasa_sidechain_overall}\n")
# ============================= OUTPUT (END) ===================================== #
