from pymatgen.core import Structure

# Load original POSCAR
struct = Structure.from_file("POSCAR")

# Multiply the cell (e.g., 2x2x2 supercell)
struct.make_supercell([3, 3, 3])

# Write out the new supercell POSCAR
struct.to(filename="POSCAR_333", fmt="poscar")