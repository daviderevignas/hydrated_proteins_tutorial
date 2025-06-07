#######################################################
TO COMPILE THE PROGRAM

gfortran -O3 PROTEIN_FOLDING_APRIL_2021.f90 -o folding.exe

if you have problems with the flag -O3 use

gfortran -O2 PROTEIN_FOLDING_APRIL_2021.f90 -o folding.exe

on the cluster you can use

ifort -O3 PROTEIN_FOLDING_APRIL_2021.f90 -o folding.exe

#######################################################
TO RUN THE PROGRAM
./folding.exe

to run in background (in this way you can keep using the terminal on your laptop, but if you close the terminal the program stops)

./folding.exe >& out &

#######################################################
INFO PROGRAM

This program runs Monte Carlo moves of lattice proteins in a coarse grained water.
You can choose if the protein is fixed or mobile. The fixed protein can be used to mimic a membrane that crosses the simulation box.
The program specifies the number of proteins (through the variable Number_Proteins) and the number of protein species (through the variable Number_Species). In particular, if two proteins share the same sequence and native structure they belong to the same species

#######################################################
INPUT FILES

1) aapot_water.dat  :  is the interaction matrix including the water-amino acid; water-membrane; amino acid-membrane terms.
The file INTERACTION_MATRIX.info contains the information on any term.
In particular the first 20 number represents the water-amino acid interaction terms, and defines if the amino acid is hydrophilic (the term is negative) or hydrophobic (the term is null or positive)
In particular, all the hydrophilic amino acid have a negative interaction term, while the hydrophobic amino acids are set to 0.
Note that the amino acids are written according to the fasta code (https://en.wikipedia.org/wiki/FASTA_format)


2) input_data_folding : contains the setup of the simulation
Box_Size                is the size of the simulation box, and should be larger than the length of any mobile protein. If one protein is fixed and stretched -- i.e. will act as a membrane -- the box size must equal to the membrane length
Temperature            
Pressure                
Number_Proteins        
Equilibration_Steps     Monte Carlo steps to equilibrate the system. No histogram is sampled during these steps. Use at least 10^6 steps (up to 10^7 should be fine)
Running_Steps           Monte Carlo steps after the equilibration. Histograms are sampled during these steps. USe at least 10^6 steps (up to 10^7 should be fine)
Sampling_Data_Time      Timesteps between two samples of the Histograms and other thermodynamic variables. Fix it such that Running_Steps / Sampling_Data_Time >=10^4 to have enough statistics
Sampling_Conf_Time      Timesteps between two samples off the protein conformations. Fix it such that Running_Steps / Sampling_Con_Time ~100 to have some conformations, but not to much as the file occupies a lot of memory
Seed                    Seed to initialize the Montecarlo move. Use a different number for any replica of the system
Number_Species          Number of different proteins
Flux                    This is a flag to establish if the protein move are isotropic (set 0) or if the protein is shifted always in the same direction (set 1)

3) input_sequences.dat : contains the sequence(s) of the protein(s). The file has two columns: the first is the ordered position of the amino acid (i.e. a number going from 1 to N where N is the length of the protein). The second is the id of the aminoacid. As if we use a negative number going from -1 to -21 according to the following scheme:

-1 A
-2 C
-3 D
-4 E
-5 F
-6 G
-7 H
-8 I
-9 K
-10 L
-11 M
-12 N
-13 P
-14 Q
-15 R
-16 S
-17 T
-18 V
-19 W
-20 Y

If you simulate more proteins, just append one sequence after the other. For example, if you simulate three proteins, the first two identical (same sequence and same native conformation) of 10 amino acids and the third of 15 amino acids, the input_sequences.dat looks like
1 -18
2 -11
3 -6
4 -11
5 -14
6 -15
7 -8
8 -6
9 -15
10 -2
1 -18
2 -11
3 -6
4 -11
5 -14
6 -15
7 -8
8 -6
9 -15
10 -2
1 -17
2 -5
3 -1
4 -13
5 -8
6 -12
7 -14
8 -1
9 -16
10 -16
11 -7
12 -14
13 -7
14 -7
15 -11

4) protein_lengths.dat  : contains the length of each protein. Following to the previous example, the protein_lengths.dat looks like
10
10
15

5) protein_species.dat : contains the "species" of proteins. Always start from 1. If two proteins are identical (same sequence and same native conformation), they belong to the same species. Following the previous example, the file looks like
1
1
2

6) protein_moves.dat : contains a flag that tell if the protein is mobile (set 1) or is fixed (set 0). When the protein is fixed, it can simulate an obstacle like a membrane or something else. In such a case, you can even choose to change the sequence of such a membrane putting all the amino acids id to -21 (it is not mandatory, in the aapot_water.dat you can tune the membrane-amino acid interaction term for any amino acid).
In the previous example, if all the proteins are mobile the protein_moves.dat looks like
1
1
1
If the first does not move, the file is like
0
1
1

7) target_structures.dat : contains the native conformations of any simulated protein. 
The file has two columns, indicating the X and Y position of the amino acids in the native state. 
If you simulate more proteins, just append one after the other all the native conformations. 
It does not matter if the coordinates of two proteins overlap, as the program analyzes separately any native conformation.

8) water_parameters : contains the water-water interaction terms. Look at the water model to understand the terms

9) Usually the program starts from a random initial configuration of proteins (always stretched) and water. The starting configuration is ALWAYS saved in the files Initial_Protein_Conformation.dat (containing the coordinates of the proteins), Initial_Configuration.dat (containing all the variables sigma of water and proteins).
If the files Final_Configuration.dat and Final_Protein_Conformation.dat are found, that are the final configuration output of a previous simulation, the program will continue from such a configuration (in such a case better to set to 0 the Equilibration_Steps of the input_data_folding).
If you want to start from a specific protein conformation and random water, just create the file Starting_Protein_Conformation.dat that contains the X and Y coordinates of all the proteins.

10) Acceptance.dat contains the acceptance ratios of the various protein moves.

REMEMBER: The order of the files is important. The first sequence must be the one of the first target structure, first species placed, first flag_move etc etc.


#################################################################
OUTPUT OF THE SIMULATION

1) Data_folding.dat contains 10 columns:
- sampling time
- energy of the system
- protein energy
- system volume
- number of water-water hydrogen bonds (HB) in the bulk
- number of water-water HB at the protein hydrophobic interfaces
- number of water-water HB at the protein hydrophilic interface
- number of water-water HB at the protein mixed (hydrophobic + hydrophilic) interface
- total number of cooperative bonds (J_sigma)
- total number of cooperative bonds at the hydrophobic interface

2) Initial_Configuration.dat, Initial_Protein_Conformation.dat, Final_Configuration.dat, Final_Protein_Conformation.dat.
The program will always write the initial and final configuration of all the system and of the proteins.
When you run the first simulation, the initial configuration of the system is generated randomly, and the proteins are placed completely stretched and possibly at equal distance from each other.
Once the simulation is finished, you can continue it if you have the files "Final_Configuration.dat" and "Final_Protein_Conformation.dat". Those files are read and used as a starting point. The program will copy them into the files "Initial_Configuration.dat" and "Initial_Protein_Conformation.dat" respectively, and start the simulation.
If you cancel the files Final_Configuration.dat and Final_Protein_Conformation.dat, the program will restart from scratch.

3) Protein_Configurations.dat
This file contains the saved configurations of the proteins.
The columns of the file are in order:
1) time
2) x of the monomer
3) y of the monomer
4) identity of the monomer (i.e. which of the 20 amino acid it is)
5) id of the monomer (i.e. which bead along the chain it is)

5) Acceptance.dat
It contains the acceptance ratio of any monte carlo move during the simulation. You can use it to have an idea of how much you are sampling the conformational space. If a ratio for a particular move is very low it means that, at the thermodynamic conditions and/or with the interaction parameters you are using, that trial move is almost always rejected.

#################################################################
VISUALIZE THE PROTEIN TRAJECTORY

After the simulation, you can plot a "movie" of the protein conformations with the script PLOTTING_PROTEIN.sh

Run it in the simulation folder where you have the file Protein_Configurations.dat with the command

bash PLOTTING_PROTEIN.sh 60 100 10000 100

where the 4 arguments of the script are
1) the size of the simulation box
2) the initial timestep from which you want to visualize the proteins. The movie will start from here.
3) the last timestep you want to visualize the proteins. The movie will finish here.
4) the timesteps between one conformation and the next (you can use Sampling_Conf_Time on the input file or multiple of it to skip some frames)
