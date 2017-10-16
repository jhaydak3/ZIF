This measures the distrubution of linker angles relative to some reference structure, presumably an energy minimized structure.

1. run createAtomList.py (for Im and BzIm linkers) or createAtomListmethyl.py (for mIm linkers) to generate a list of atom ndx's from a single snapshot xyz file of structure being looked at. This ignores linkers near boundary because during the simulation they may cross the boundary and that makes it hard to calculate the angle. 

2. Run getSwingAnglewf.py and change the parameters in the main() section. Requires an xyz file from an MD simulation of the structure being considered and atom list from step 1. This will output a file with distribution of angles.
