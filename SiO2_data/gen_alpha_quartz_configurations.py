# this script generates configurations based on alpha-quartz

from pymatgen.core import Structure

import numpy as np
from pymatgen.core import Lattice
from pathlib import Path





# def perturb_lattice(structure: Structure, eps_max: float = 0.1) -> Structure:
#     # 1. Generate 3x3 matrix of random epsilons
#     eps = np.random.uniform(-eps_max, eps_max, size=(3, 3))
#     print(f"eps={eps}")
#     # 2. Construct ones matrix and add eps to each element: (1 + eps)
#     scale_matrix = np.ones((3, 3)) + eps
#     # 3. Hadamard product (element-wise multiplication '*') with original lattice
#     new_lattice = structure.lattice.matrix * scale_matrix
#     print(f"new_lattice={new_lattice}")
#     # 4. Return new Structure keeping fractional coordinates invariant
#     return Structure(
#         lattice=new_lattice,
#         species=structure.species,
#         coords=structure.frac_coords,
#         coords_are_cartesian=False,
#     )



def increase_length(
        structure: Structure,
        length_eps_max: float = 0.1,
        frac_eps_max: float = 0.02,
        min_atom_dist: float = 1.0,
) -> Structure:
    """Expands lattice vector lengths (a, b, c) by up to +length_eps_max and

    randomly perturbs atomic positions while preventing collisions (d >= min_atom_dist).
    """
    # 1. Expand lattice parameters
    a, b, c, alpha, beta, gamma = structure.lattice.parameters
    length_eps_vec = np.random.uniform(0, length_eps_max, size=3)
    a_new = a * (1.0 + length_eps_vec[0])
    b_new = b * (1.0 + length_eps_vec[1])
    c_new = c * (1.0 + length_eps_vec[2])

    # 2. Build the new expanded Lattice object
    new_lattice = Lattice.from_parameters(
        a_new, b_new, c_new, alpha, beta, gamma
    )

    # 3. Resample fractional coordinates until no atoms collide
    frac_coords = structure.frac_coords
    while True:
        frac_displacement = np.random.uniform(
            -frac_eps_max, frac_eps_max, size=frac_coords.shape
        )
        frac_coords_new = (frac_coords + frac_displacement) % 1.0

        perturbed_struct = Structure(
            lattice=new_lattice,
            species=structure.species,
            coords=frac_coords_new,
            coords_are_cartesian=False,
            to_unit_cell=True,
        )

        dists = perturbed_struct.distance_matrix[
            np.triu_indices(len(perturbed_struct), k=1)
        ]
        if np.min(dists) >= min_atom_dist:
            return perturbed_struct


def decrease_length(
        structure: Structure,
        length_eps_max: float = 0.1,
        frac_eps_max: float = 0.02,
        min_atom_dist: float = 1.0,
) -> Structure:
    """Compresses lattice vector lengths (a, b, c) by up to -length_eps_max and

    randomly perturbs atomic positions while preventing collisions (d >= min_atom_dist).
    """
    # 1. Compress lattice parameters
    a, b, c, alpha, beta, gamma = structure.lattice.parameters
    length_eps_vec = np.random.uniform(-length_eps_max, 0, size=3)
    a_new = a * (1.0 + length_eps_vec[0])
    b_new = b * (1.0 + length_eps_vec[1])
    c_new = c * (1.0 + length_eps_vec[2])

    # 2. Build the new compressed Lattice object
    new_lattice = Lattice.from_parameters(
        a_new, b_new, c_new, alpha, beta, gamma
    )

    # 3. Resample fractional coordinates until no atoms collide
    frac_coords = structure.frac_coords
    while True:
        frac_displacement = np.random.uniform(
            -frac_eps_max, frac_eps_max, size=frac_coords.shape
        )
        frac_coords_new = (frac_coords + frac_displacement) % 1.0

        perturbed_struct = Structure(
            lattice=new_lattice,
            species=structure.species,
            coords=frac_coords_new,
            coords_are_cartesian=False,
            to_unit_cell=True,
        )

        dists = perturbed_struct.distance_matrix[
            np.triu_indices(len(perturbed_struct), k=1)
        ]
        if np.min(dists) >= min_atom_dist:
            return perturbed_struct

def perturb_all(
        structure: Structure,
        length_eps_max: float = 0.1,
        frac_eps_max: float = 0.02,  # Recommended: 0.01 - 0.02 for ~0.1 - 0.3 Å
        angle_eps_max: float = 0.05,
        min_atom_dist: float = 1.0,  # Minimum allowed interatomic distance in Å
) -> Structure:
    """Perturbs cell lengths, cell angles, and atomic positions without allowing

    any two atoms to collide (d >= min_atom_dist).
    """
    a, b, c, alpha, beta, gamma = structure.lattice.parameters

    # 1. Perturb cell lengths
    length_eps_vec = np.random.uniform(-length_eps_max, length_eps_max, size=3)
    a_new = a * (1.0 + length_eps_vec[0])
    b_new = b * (1.0 + length_eps_vec[1])
    c_new = c * (1.0 + length_eps_vec[2])

    # 2. Resample angles until a valid 3D lattice is obtained
    new_lattice = None
    while True:
        angle_eps_vec = np.random.uniform(-angle_eps_max, angle_eps_max, size=3)
        alpha_new = alpha * (1.0 + angle_eps_vec[0])
        beta_new = beta * (1.0 + angle_eps_vec[1])
        gamma_new = gamma * (1.0 + angle_eps_vec[2])

        if not (
                0.0 < alpha_new < 180.0
                and 0.0 < beta_new < 180.0
                and 0.0 < gamma_new < 180.0
        ):
            continue

        try:
            cand_lattice = Lattice.from_parameters(
                a_new, b_new, c_new, alpha_new, beta_new, gamma_new
            )
            if cand_lattice.volume > 1e-3 and not np.isnan(cand_lattice.volume):
                new_lattice = cand_lattice
                break
        except (ValueError, ZeroDivisionError):
            continue

    # 3. Resample fractional coordinates until no atoms overlap
    frac_coords = structure.frac_coords
    while True:
        frac_displacement = np.random.uniform(
            -frac_eps_max, frac_eps_max, size=frac_coords.shape
        )
        frac_coords_new = (frac_coords + frac_displacement) % 1.0

        perturbed_struct = Structure(
            lattice=new_lattice,
            species=structure.species,
            coords=frac_coords_new,
            coords_are_cartesian=False,
            to_unit_cell=True,
        )

        # Check periodic pairwise distance: distance_matrix diagonal is 0,
        # so we take only the upper triangle (k=1)
        dists = perturbed_struct.distance_matrix[
            np.triu_indices(len(perturbed_struct), k=1)
        ]
        if np.min(dists) >= min_atom_dist:
            return perturbed_struct


poscar_file_name="POSCAR_alpha_quartz333"
structure = Structure.from_file(poscar_file_name)

incar_file_name="INCAR_alpha_quartz_333"
potcar_file_name="POTCAR"
kpoints_file_name="KPOINTS_alpha_quartz_333"

N_increase=5
N_decrease=5
N_perturb_all=5
out_dir_base="alpha_quartz"


with open(incar_file_name, "r") as f:
    incar_content = f.read()

with open(potcar_file_name,"r") as f:
    potcar_content=f.read()

with open(kpoints_file_name,"r") as f:
    kpoints_content=f.read()



def gen_dir_increase_length(incar_content,
                            potcar_content,
                            kpoints_content,
                            out_dir_base,
                            ind,
                            structure: Structure,
                            length_eps_max: float = 0.1,
                            frac_eps_max: float = 0.02,
                            min_atom_dist: float = 1.0,):
    """

    Args:
        incar_content:
        potcar_content:
        kpoints_content:
        out_dir_base:
        ind:

    Returns:

    """
    out_dir=out_dir_base+f"/increase_length_{ind}"
    Path(out_dir).mkdir(exist_ok=True,parents=True)

    out_incar_file_name=out_dir+"/INCAR"
    out_potcar_file_name=out_dir+"/POTCAR"
    out_kpoints_file_name=out_dir+"/KPOINTS"
    out_poscar_file_name=out_dir+"/POSCAR"

    new_struct = increase_length(
        structure,
        length_eps_max=length_eps_max,
        frac_eps_max=frac_eps_max,
        min_atom_dist=min_atom_dist,
    )

    new_struct.to(filename=out_poscar_file_name, fmt="poscar")

    Path(out_incar_file_name).write_text(incar_content)
    Path(out_potcar_file_name).write_text(potcar_content)
    Path(out_kpoints_file_name).write_text(kpoints_content)

    print(f"Generated: {out_dir}")
    return out_dir



def gen_dir_decrease_length(
        incar_content: str,
        potcar_content: str,
        kpoints_content: str,
        out_dir_base: str,
        ind: int,
        structure: Structure,
        length_eps_max: float = 0.1,
        frac_eps_max: float = 0.02,
        min_atom_dist: float = 1.0,
):
    out_dir = f"{out_dir_base}/decrease_length_{ind}"
    Path(out_dir).mkdir(exist_ok=True, parents=True)

    out_incar_file_name = f"{out_dir}/INCAR"
    out_potcar_file_name = f"{out_dir}/POTCAR"
    out_kpoints_file_name = f"{out_dir}/KPOINTS"
    out_poscar_file_name = f"{out_dir}/POSCAR"

    new_struct = decrease_length(
        structure,
        length_eps_max=length_eps_max,
        frac_eps_max=frac_eps_max,
        min_atom_dist=min_atom_dist,
    )

    new_struct.to(filename=out_poscar_file_name, fmt="poscar")

    Path(out_incar_file_name).write_text(incar_content)
    Path(out_potcar_file_name).write_text(potcar_content)
    Path(out_kpoints_file_name).write_text(kpoints_content)

    print(f"Generated: {out_dir}")
    return out_dir


def gen_dir_perturb_all(
        incar_content: str,
        potcar_content: str,
        kpoints_content: str,
        out_dir_base: str,
        ind: int,
        structure: Structure,
        length_eps_max: float = 0.1,
        frac_eps_max: float = 0.02,
        angle_eps_max: float = 0.05,
        min_atom_dist: float = 1.0,
):
    out_dir = f"{out_dir_base}/perturb_all_{ind}"
    Path(out_dir).mkdir(exist_ok=True, parents=True)

    out_incar_file_name = f"{out_dir}/INCAR"
    out_potcar_file_name = f"{out_dir}/POTCAR"
    out_kpoints_file_name = f"{out_dir}/KPOINTS"
    out_poscar_file_name = f"{out_dir}/POSCAR"

    new_struct = perturb_all(
        structure,
        length_eps_max=length_eps_max,
        frac_eps_max=frac_eps_max,
        angle_eps_max=angle_eps_max,
        min_atom_dist=min_atom_dist,
    )

    new_struct.to(filename=out_poscar_file_name, fmt="poscar")

    Path(out_incar_file_name).write_text(incar_content)
    Path(out_potcar_file_name).write_text(potcar_content)
    Path(out_kpoints_file_name).write_text(kpoints_content)

    print(f"Generated: {out_dir}")
    return out_dir

# ==============================================================================
# Run generation loops
# ==============================================================================

# 1. Generate increase_length configurations
for i in range(N_increase):
    gen_dir_increase_length(
        incar_content=incar_content,
        potcar_content=potcar_content,
        kpoints_content=kpoints_content,
        out_dir_base=out_dir_base,
        ind=i,
        structure=structure,
    )


# 2. Generate decrease_length configurations
for i in range(N_decrease):
    gen_dir_decrease_length(
        incar_content=incar_content,
        potcar_content=potcar_content,
        kpoints_content=kpoints_content,
        out_dir_base=out_dir_base,
        ind=i,
        structure=structure,
    )


# 3. Generate perturb_all configurations
for i in range(N_perturb_all):
    gen_dir_perturb_all(
        incar_content=incar_content,
        potcar_content=potcar_content,
        kpoints_content=kpoints_content,
        out_dir_base=out_dir_base,
        ind=i,
        structure=structure,
    )