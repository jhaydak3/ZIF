This program works in two parts. A short description of what it does: it takes the structure (single snapshot, xyz file) of a ZIF framework and a file containing information on the trajectories of adsorbants and produces an overlay of trajectory path on the framework in two dimensions (xy,yz,xz). It is a density plot - the framework is split into an N by N grid and for each snapshot of the MD simulation each grid point is samples to see if each adsorbant is within a certain cutoff range of the grid point. For every time there is an absorbant close to the gridpoint, the value of that gridpoint is increased by 1. Please see sample output.

Each file contains various settings and inputs. Inputs are hard-coded into the file in the main() section usually.

The way to generate a 2d density plot: 

-have xyz file of empty ZIF structure
-run make2dplot.py on a MSD trajectory file, select which atom to track, produces a new output trajectory file
-run makedensityplot.py on new trajectory file. In the main() section, set inputs (xyz file, traj file, cell parameters, gridpoints, axis, etc). This will generate several .txt files that get fed to a matlab script. (This function imports a fortran module so it runs faster, takes very long without this)
-run plotZn.m. Set inputs in first few lines of code.






Here is a sample of a trajectory file format for the adsorbants:
COM stands for center of mass. The absorbant laballed #1 is the entire framework itself.

%%%%%
# Time-averaged data for fix 1
# TimeStep Number-of-rows
# Row c_COMsorbate[1] c_COMsorbate[2] c_COMsorbate[3]
200000 7
1 10.9507 19.4689 15.936
2 -6.09933 29.5416 16.0792
3 22.1313 25.8135 21.4737
4 7.90769 21.1487 22.2817
5 -0.703876 0.166057 23.786
6 -0.309408 33.0587 26.9071
7 -10.8316 22.724 25.6012
200010 7
1 10.952 19.4689 15.9359
2 -6.10276 29.486 16.1056
3 22.1047 25.8446 21.4484
4 7.91615 21.1673 22.3583
5 -0.742162 0.158077 23.7567
6 -0.323982 33.0469 26.953
7 -10.8481 22.6856 25.5851
200020 7
1 10.9534 19.469 15.9358
2 -6.10468 29.4327 16.1285
3 22.0807 25.8777 21.4238
4 7.92506 21.1859 22.4343
5 -0.780356 0.149702 23.7279
6 -0.342113 33.0341 26.9966
7 -10.866 22.6493 25.5675
%%%
