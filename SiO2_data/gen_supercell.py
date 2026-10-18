from pymatgen.core import Structure

# Load original POSCAR
in_poscar_file_name="POSCAR_alpha_quartz"
struct = Structure.from_file(in_poscar_file_name)

# Multiply the cell (e.g., 2x2x2 supercell)
n0=3
n1=3
n2=3
struct.make_supercell([n0, n1, n2])

# Write out the new supercell POSCAR
struct.to(filename=in_poscar_file_name+f"{n0}{n1}{n2}", fmt="poscar")