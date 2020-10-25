# optimizeOERmechanism

This repository enables one to enumerate all possible electrocatalytic pathways underlying the oxygen evolution reaction (OER) on a given electrocatalyst surface. It uses the Gurobi API in Python to accomplish this task. If you use our code, please cite: 

Govind Rajan, A.; Carter, E. A. Discovering Competing Electrocatalytic Mechanisms and their Overpotentials: Automated Enumeration of Oxygen Evolution Pathways. J. Phys. Chem. C, 2020. DOI: 10.1021/acs.jpcc.0c08120

The code is written by: Ananth Govind Rajan and Garrett R. Dowdy

The file is to be executed as: python optimizeOERmechanism.py (name-of-species-library-file) (number-of-reactions-to-consider), e.g.,

python optimizeOERmechanism.py OER_NiOOH_0001_Oxo_HSE.txt 6

There are 4 different species library files already available, these are:
OER_NiOOH_0001_Oxo_PBE.txt,
OER_NiOOH-Fe_-12-11_Hydroxo_HSE.txt,
OER_NiOOH_0001_MultiOxo_PBE.txt, and
OER_NiOOH_0001_Oxo_HSE.txt

The crystallographic facet, active site, and level of DFT theory are indicated in each filename.

