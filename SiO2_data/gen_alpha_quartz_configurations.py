# this script generates configurations based on alpha-quartz

from pymatgen.core import Structure
in_file_name="POSCAR_alpha_quartz333"
import numpy as np

# 1. Load the structure
structure = Structure.from_file(in_file_name)
# -----------------------------------------------------------------
# 2. Get the 3 basis vectors (Lattice Vectors)
# -----------------------------------------------------------------
# lattice.matrix is a 3x3 numpy array where each row is a vector: [a, b, c]
lattice_matrix = structure.lattice.matrix

a_vec = lattice_matrix[0]  # vector a
b_vec = lattice_matrix[1]  # vector b
c_vec = lattice_matrix[2]  # vector c


# -----------------------------------------------------------------
# 3. Get all atomic coordinates (as full Nx3 arrays)
# -----------------------------------------------------------------
# Fractional coordinates (Direct, dimensionless [0, 1)):

frac_coords = structure.frac_coords

def perturb_lattice(structure: Structure, eps_max: float = 0.1) -> Structure:
    # 1. Generate 3x3 matrix of random epsilons
    eps = np.random.uniform(-eps_max, eps_max, size=(3, 3))
    print(f"eps={eps}")
    # 2. Construct ones matrix and add eps to each element: (1 + eps)
    scale_matrix = np.ones((3, 3)) + eps
    # 3. Hadamard product (element-wise multiplication '*') with original lattice
    new_lattice = structure.lattice.matrix * scale_matrix
    print(f"new_lattice={new_lattice}")
    # 4. Return new Structure keeping fractional coordinates invariant
    return Structure(
        lattice=new_lattice,
        species=structure.species,
        coords=structure.frac_coords,
        coords_are_cartesian=False,
    )

st1=perturb_lattice(structure,0.1)