1. cat ~/Documents/cppCode/vasp640/POT/potpaw_PBE.54/Si/POTCAR >> POTCAR  
   cat ~/Documents/cppCode/vasp640/POT/potpaw_PBE.54/O/POTCAR >> POTCAR
    
2. generating supercell
   pip install pymatgen
   python gen_supercell.py
    mv POSCAR_XXX to working directory