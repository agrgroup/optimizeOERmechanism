# Mechanism enumeration for the oxygen evolution reaction using Gurobi and Python
# Ref: Govind Rajan, A.; Carter, E. A. Discovering Competing Electrocatalytic Mechanisms and their Overpotentials: Automated Enumeration
# of Oxygen Evolution Pathways. J. Phys. Chem. C, 2020. DOI: 10.1021/acs.jpcc.0c08120
# Code written by: Ananth Govind Rajan, Garrett R. Dowdy
# The file is to be executed as: python optimizeOERmechanism.py <name-of-species-library-file> <number-of-reactions-to-consider>, e.g.,
# python optimizeOERmechanism.py OER_NiOOH_0001_Oxo_HSE.txt 6

# Start timing the code
import timeit
tic = timeit.default_timer();

# Define the species class
class Species:
    def __init__(self, Name, n_H, n_O, G, c):

        self.Name = Name
        self.n_H = n_H  # number of H atoms
        self.n_O = n_O  # number of O atoms
        self.G = G   # Gibbs free energy (eV)
        self.c = c  # Coefficient in the overall OER reaction
        self.InvolvesSite = ('*' in Name)

import sys

# Retrieve the species dictionary from the text file name provided as the first command line argument 
execfile(sys.argv[1])

# Obtain the number of reactions from the second command line argument
n_R = int(sys.argv[2])

# Specify the maximum number of mechanisms you want to see, currently set to an arbitrarily high value
NumSolutionsToFind = 100000

# Initialize the model
from gurobipy import *
m = Model()

# Specify the parameters related to the "Solution Pool"
m.setParam('OutputFlag', True)
m.setParam('PoolSearchMode',2)  # See the gurobi documentation.  A value of 2 indicates that Gurobi will do a systematic search for the n best solutions
m.setParam('PoolSolutions', NumSolutionsToFind)  # Specifies the number of solutions you want the solver to look for
m.setParam('MIPGap',0) # Specifies the tolerance for declaring two objective function values to be the same

# Declare the sets
S = Species.keys()
S_site = [s for s in Species.keys() if Species[s].InvolvesSite == True]
S_molecule = [s for s in Species.keys() if ( (Species[s].InvolvesSite == False) and (Species[s].Name != '(H++e-)') )]   
R = range(1, n_R+1)
R_subset =  range(1, n_R)

# Declare the decision variables
g = m.addVar(vtype = GRB.CONTINUOUS, name = 'g')

v_products = m.addVars(S, R, vtype = GRB.BINARY, name = 'ProductCoefficients')
v_reactants = m.addVars(S, R, vtype = GRB.BINARY, name = 'ReactantCoefficients')

# Declare the constraints
m.addConstrs( (sum([(v_products[s,r] - v_reactants[s,r]) * Species[s].n_H for s in S]) == 0 for r in R), "H_balances"  )

m.addConstrs( (sum([(v_products[s,r] - v_reactants[s,r]) * Species[s].n_O for s in S]) == 0 for r in R), "O_balances"  )

m.addConstrs( (sum([v_products[s,r] - v_reactants[s,r] for s in S_site]) == 0 for r in R), "Site_balances")

m.addConstrs( (sum([(v_products[s,r] - v_reactants[s,r]) for r in R]) == Species[s].c for s in S), "OverallCoefficientBalances" )

G_Values = [s.G for s in Species.values()]
M = 2*(max(G_Values) - min(G_Values))
m.addConstrs( ((g >= sum([(v_products[s,r] - v_reactants[s,r]) * Species[s].G for s in S])-M*(1-v_products['(H++e-)',r])) for r in R), "MaxInequalities" )

m.addConstrs( (v_reactants[s,r] + v_products[s,r] <= 1 for s in S for r in R), "ReactanctXorProduct" )

m.addConstrs( (sum([v_reactants[s,r] for s in S]) <= 2 for r in R), "AtMost2Reactants" )

m.addConstrs( (sum([v_products[s,r] for s in S])-v_products['O2',r] <= 2 for r in R), "AtMost2ProductsExceptIfO2Released" )

m.addConstrs( (sum([v_reactants[s,r] for s in S_site]) <=1 for r in R), "AtMost1SiteInReactants" ) 

m.addConstrs( (sum([v_products[s,r] for s in S_site]) <=1 for r in R), "AtMost1SiteInProducts" )

m.addConstrs(  (v_reactants['(H++e-)',r] == 0 for r in R), "ElectronCannotBeReactantinOxidation" )

m.addConstrs(  (v_reactants['O2',r] == 0 for r in R), "O2CannotBeReactantinOER" )

m.addConstr( (v_reactants[InitialStateSurface, 1] == 1), "ReactionOrderingFirstSpecies")
m.addConstrs( ((v_reactants[s,r+1] >=  v_products[s,r] ) for s in S_site for r in R_subset), "ReactionOrderingIntermediateReactions")

m.addConstrs( (((sum([(v_products[s,r] - v_reactants[s,r]) * Species[s].G for s in S])) <= MaxGibbsNonElectroactive+M*v_products['(H++e-)',r]) for r in R), "ExcludeHighEnergyNonElectroactiveStep")

# Appear only once in cycle as reactant/product to avoid loops
m.addConstrs( (sum([v_reactants[s,r] for r in R])<=1 for s in S_site), "StateOnlyOnceReactant" )
m.addConstrs( (sum([v_products[s,r] for r in R])<=1 for s in S_site), "StateOnlyOnceProduct" )

# Declare the objective
m.setObjective(g, GRB.MINIMIZE)

# Solve the model
m.update()
m.optimize()

#toc = timeit.default_timer();
#print 'Time elapsed: ', round(toc-tic,1), 'seconds'
#exit();

if (m.SolCount==0):
	toc = timeit.default_timer();
	print 'Time elapsed: ', round(toc-tic,1), 'seconds'
	exit();

# Import Numpy
import numpy as np
numExcludedMechanisms=0
n_S = sum(1 for s in S)

# Print out the results
# Loop over all the solutions found
etaList = np.zeros((m.SolCount,1))
gHNEGSList = np.zeros((m.SolCount,1))
indices = np.zeros((m.SolCount,1))
for SolInd in range(0, m.SolCount):
	# Display an update to inform the user which solution we are looking at
	print('\nSolution ' + str(SolInd))

        # Set the Gurobi parameter to indicate which solution we are looking at
        m.setParam('SolutionNumber', SolInd)

	gmaxNonElectro=-np.inf
        ExcludeMechanism = False
        for r in R:
		DelG = sum([(v_products[s,r].xn - v_reactants[s,r].xn) * Species[s].G for s in S])
		
		if (DelG>gmaxNonElectro and round(v_products['(H++e-)',r].xn)==0):
			gmaxNonElectro=DelG


	gHNEGSList[SolInd]=round(gmaxNonElectro,2)
	gmaxElectro = g.xn
        eta = gmaxElectro - DelGTot
        etaList[SolInd] = round(eta,2)
	indices[SolInd]=SolInd

etaGHNEGSList=np.hstack((indices,etaList,gHNEGSList))

from operator import itemgetter
sortedEtaGHNEGSList = np.array(sorted(etaGHNEGSList, key=itemgetter(1,2)))

indices = sortedEtaGHNEGSList[:,0]

for SolInd in range(0, m.SolCount):
        # Display an update to inform the user which solution we are looking at
        print('\nSolution ' + str(SolInd))

	ActualSolInd = int(indices[SolInd])

        # Set the Gurobi parameter to indicate which solution we are looking at
        m.setParam('SolutionNumber', ActualSolInd)

        if (ExcludeMechanism == True):
                numExcludedMechanisms=numExcludedMechanisms+1
 		print "Exclude mechanism"
	else:
		# Loop over the reactions
		for r in R:
			print('Reaction ' + str(r))
		
			# Build up the reactant string
			ReactantString = ''
			FirstSpecies = True
			for s in S:
				if round(v_reactants[s,r].xn) == 1:
					if FirstSpecies == False:
						ReactantString = ReactantString + ' + '
					else:
						FirstSpecies = False
					ReactantString = ReactantString + s
		
			# Build up the product string
			ProductString = ''
			FirstSpecies = True
			for s in S:
				if round(v_products[s,r].xn) == 1:
					if FirstSpecies == False:
						ProductString = ProductString + ' + '
					else:
						FirstSpecies = False
					ProductString = ProductString + s
			
			RxnString = '\t' + ReactantString + ' --> ' + ProductString
			DelG = sum([(v_products[s,r].xn - v_reactants[s,r].xn) * Species[s].G for s in S])
			DelGString = 'DelG = ' + str(round(DelG,2)) + ' eV'	
		
			if (round(DelG,2)==round(g.xn,2) and (round(v_reactants['(H++e-)',r].xn)==1 or round(v_products['(H++e-)',r].xn)==1)):
				PDSstring = ' (PDS)'
			else:
				PDSstring = ''
		
			if (round(DelG,2)==round(gHNEGSList[ActualSolInd],2) and (round(v_reactants['(H++e-)',r].xn)==0 and round(v_products['(H++e-)',r].xn)==0)):
		                nonElectroString = ' (HNEGS)'
		        else:
		                nonElectroString = ''
		
				
			# Print the full string
			print(RxnString.ljust(45) + DelGString + PDSstring + nonElectroString)

		# Calculate the maximum g value
		gmaxOverall = max([sum([(v_products[s,r].xn - v_reactants[s,r].xn) * Species[s].G for s in S]) for r in R])
		gmaxElectro = g.xn
		eta = gmaxElectro - DelGTot
		etaList[ActualSolInd] = eta
		
		print('\nMaximum DelG value = ' + str(round(gmaxOverall,2)) + ' eV')    
		print('Maximum Electroactive DelG value = ' + str(round(gmaxElectro,2)) + ' eV')
		print('Maximum Non-Electroactive DelG value = ' + str(round(gHNEGSList[ActualSolInd],2)) + ' eV')
		print('Overpotential = ' + str(round(eta,2)) + ' V\n')

print('Number of mechanisms requested: ' + str(NumSolutionsToFind))
print('Total number of mechanisms found: ' + str(m.SolCount))
#print('Number of unique mechanisms (i.e., not found amongst mechanisms with fewer number of steps): ' + str(m.SolCount-numExcludedMechanisms))
print('Lowest overpotential found: ' + str(round(g.x-DelGTot,2)) + ' V')

fileName = sys.argv[1]
fileName=fileName.replace(".txt","_eta.txt")
fileObj = open(fileName,"a")
etaGHNEGSList = np.hstack((etaList, gHNEGSList))
for eta,gHNEGS in etaGHNEGSList:
	eta=str(round(eta,2))
	eta=eta.replace("[","")
	eta=eta.replace("]","")

        gHNEGS=str(round(gHNEGS,2))
        gHNEGS=gHNEGS.replace("[","")
        gHNEGS=gHNEGS.replace("]","")

	fileObj.write(str(eta) + "\t" + str(gHNEGS))
        fileObj.write("\n")	

toc = timeit.default_timer();
print 'Time elapsed: ', round(toc-tic,1), 'seconds'
