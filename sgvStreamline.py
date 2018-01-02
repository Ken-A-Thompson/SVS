#Authors: Matthew Osmond <mmosmond@zoology.ubc.ca> & Ken A. Thompson <ken.thompson@zoology.ubc.ca>
#Description: The role of standing genetic variance in hybrid load

import numpy as np
import time
import matplotlib.pyplot as plt

######################################################################
##UNIVERSAL PARAMETERS##
######################################################################

n = 2 #phenotypic dimensions (positive integer)

######################################################################
##PARAMETERS TO MAKE ANCESTOR##
######################################################################

K = 1000 #number of individuals (positive integer)
n_muts = 10 #number of mutations (positive integer)
p_mut = 0.1 #probability of having mutation at any one locus (0<p<1) #set this to zero for de novo only
alpha = 0.1 #mutational sd (positive real number)

######################################################################
##MAKE ANCESTOR##
######################################################################

pop = np.random.binomial(1, p_mut, (K, n_muts)) #p_mut chance of having each of n_muts mutations, for all K individuals
mut = np.random.normal(0, alpha, (n_muts, n)) #create n_muts mutations, each with a random normal phenotypic effect in each n dimension with mean 0 and sd alpha

######################################################################
##PARAMETERS FOR ADAPTING POPULATIONS##
######################################################################

K_adapt = K #number of individuals (positive integer)
alpha_adapt = alpha #mutational sd (positive real number)
B = 2 #number of offspring per generation per parent (positive integer)
u = 0.001 #mutation probability per generation per genome (0<u<1)

theta1 = np.array([0.5] * n) #optimum phenotype for population 1
theta2 = np.array([-0.5] * n) #optimum phenotype for population 2

maxgen = 1000 #total number of generations populations adapt for
outputFreq = 10 #generation interval between plotting

remove_lost = True #If true, remove mutations that are lost (0 for all individuals)
remove = 'derived' #.. any derived (not from ancestor) mutation that is lost 

######################################################################
##FOUND ADAPTING POPULATIONS FROM ANCESTOR##
######################################################################

#population 1
whofounds = np.random.choice(K, size = K_adapt, replace = False) #random choice of nfounders from ancestral population 
popfound1 = pop[whofounds] #list of mutations held by each founding individual
if remove_lost and remove == 'any': #if removing ancestral mutations when lost
	keep = popfound1.any(axis=0)
	mutfound1 = mut[keep]
	popfound1 = pop[:, keep]
else:
	mutfound1 = mut

#population 2
whofounds = np.random.choice(K, size = K_adapt, replace = False) #random choice of nfounders from ancestral population 
popfound2 = pop[whofounds] #list of mutations held by each founding individual
if remove_lost and remove == 'any': #if removing ancestral mutations when lost
	keep = popfound2.any(axis=0)
	mutfound2 = mut[keep]
	popfound2 = pop[:, keep]
else:
	mutfound2 = mut

######################################################################
##PARAMETERS FOR HYBRIDS##
######################################################################

nHybrids = 10

######################################################################
##FUNCTION FOR POPULATIONS TO ADAPT##
######################################################################

def main():

	#initialize adapting populations
	pop1 = popfound1
	mut1 = mutfound1

	pop2 = popfound2
	mut2 = mutfound2

	hyload_old = 0

	gen = 0 #generation
	while gen < maxgen + 1:

		# make hybrids every outputFreq generations
		if gen % outputFreq == 0:

			offphenos = dict()
			offpheno = []
			for k in range(nHybrids):
			    # random parents
			    randpar1 = pop1[np.random.choice(len(pop1))] 
			    randpar2 = pop2[np.random.choice(len(pop2))]
			    # random parent phenotypes
			    phenpar1 = np.dot(randpar1, mut1) 
			    phenpar2 = np.dot(randpar2, mut2)
			    # mutations held by random parents
			    mutpar1 = mut1 * randpar1[:, None]
			    mutpar2 = mut2 * randpar2[:, None]
			    setA = set(tuple(x) for x in mutpar1)
			    setB = set(tuple(x) for x in mutpar2)
			    # mutations shared by two parents (all in offspring)
			    sharedmuts = np.array([x for x in setA & setB])
			    # mutations not shared by two parents
			    unsharedmuts = np.array([x for x in setA ^ setB])
			    # which unshared mutations in offspring (free recombination between all loci, therefore gets each with 0.5 probability)
			    randmuts = np.random.randint(2, size = (len(unsharedmuts)))
			    offmuts = unsharedmuts * randmuts[:, None]
			    if len(offmuts) < 1:
			        offmuts = np.array([[0] * n]) #give something in case empty
			    # offspring phenotype
			    offpheno.append(sum(np.append(sharedmuts, offmuts, axis = 0)))

			offpheno = np.array(offpheno) #reformat correctly
			dist = np.linalg.norm(offpheno - np.mean(offpheno, axis=0), axis=1) #phenotypic distance from mean hybrid
			w = np.exp(-dist**2) #probability of survival
			hyload = np.log(1*B) - np.mean(np.log(w*B)) #hybrid load as defined by Chevin et al 2014
			
			#print some data (gen, pheno var in each dimension, hybrid load)
			print(gen, np.var(offpheno, axis=0), hyload)
			
			#plot hybrid load
			plt.axis([0, maxgen, 0, 0.2])
			plt.ion()

			plt.plot([gen-outputFreq, gen], [hyload_old,hyload], '-o', color='k')
			plt.pause(0.01)
			hyload_old = hyload

			# while True:
			# 	plt.pause(0.05)

		# genotype to phenotype for adapting populations
		phenos1 = np.dot(pop1, mut1) #sum mutations held by each individual
		phenos2 = np.dot(pop2, mut2) #sum mutations held by each individual

		# viability selection for adapting populations
		phenos = phenos1
		theta = theta1
		pop = pop1
		dist = np.linalg.norm(phenos - theta, axis=1) #phenotypic distance from optimum
		w = np.exp(-dist**2) #probability of survival
		rand = np.random.uniform(size = len(pop)) #random uniform number in [0,1] for each individual
		surv = pop[rand < w] #survivors
		if len(surv) > K_adapt:
			surv = surv[np.random.randint(len(surv1), size = K_adapt)] #randomly choose K_adapt individuals if more than K_adapt
		surv1 = surv

		phenos = phenos2
		theta = theta2
		pop = pop2
		dist = np.linalg.norm(phenos - theta, axis=1) #phenotypic distance from optimum
		w = np.exp(-dist**2) #probability of survival
		rand = np.random.uniform(size = len(pop)) #random uniform number in [0,1] for each individual
		surv = pop[rand < w] #survivors
		if len(surv) > K_adapt:
			surv = surv[np.random.randint(len(surv), size = K_adapt)] #randomly choose K_adapt individuals if more than K_adapt
		surv2 = surv
				
		#end simulation if either adapting population extinct        
		if len(surv1) == 0 or len(surv2) == 0: 
			print("Extinct")              
			break 
			
		# birth for adapting populations	
		surv = surv1
		pairs = np.resize(np.random.choice(len(surv), size=len(surv), replace=False), (int(len(surv)/2), 2)) #random mate pairs (each mates at most once and not with self)
		rand2 = np.random.randint(2, size=(len(pairs), len(surv[0]))) #from which parent each offspring inherits each allele (free recombination, fair transmission)
		rec = np.resize(np.append(rand2, 1-rand2, axis=1),(len(rand2), 2, len(rand2[0]))) #reshape
		off_1 = np.sum(surv[pairs] * rec, axis=1) #one product of meiosis
		off_2 = np.sum(surv[pairs] * (1-rec), axis=1) #other product of meiosis
		off = np.repeat(np.append(off_1, off_2, axis=0), B, axis=0) #each product of meiosis produced B times
		off1 = off

		surv = surv2
		pairs = np.resize(np.random.choice(len(surv), size=len(surv), replace=False), (int(len(surv)/2), 2)) #random mate pairs (each mates at most once and not with self)
		rand2 = np.random.randint(2, size=(len(pairs), len(surv[0]))) #from which parent each offspring inherits each allele (free recombination, fair transmission)
		rec = np.resize(np.append(rand2, 1-rand2, axis=1),(len(rand2), 2, len(rand2[0]))) #reshape
		off_1 = np.sum(surv[pairs] * rec, axis=1) #one product of meiosis
		off_2 = np.sum(surv[pairs] * (1-rec), axis=1) #other product of meiosis
		off = np.repeat(np.append(off_1, off_2, axis=0), B, axis=0) #each product of meiosis produced B times
		off2 = off

		# mutation for adapting populations
		off = off1
		rand3 = np.random.uniform(size = len(off)) #random uniform number in [0,1] for each offspring
		nmuts = sum(rand3 < u) # mutate if random number is below mutation rate; returns number of new mutations
		whomuts1 = np.where(rand3 < u) #indices of mutants
		newmuts = np.random.normal(0, alpha, (nmuts, n)) #phenotypic effect of new mutations
		newmuts1 = newmuts

		off = off2
		rand3 = np.random.uniform(size = len(off)) #random uniform number in [0,1] for each offspring
		nmuts = sum(rand3 < u) # mutate if random number is below mutation rate; returns number of new mutations
		whomuts2 = np.where(rand3 < u) #indices of mutants
		newmuts = np.random.normal(0, alpha, (nmuts, n)) #phenotypic effect of new mutations
		newmuts2 = newmuts

		# update adapting populations
		pop1 = np.append(off1, np.transpose(np.identity(len(off1), dtype=int)[whomuts1[0]]), axis=1) #add new loci and identify mutants
		mut1 = np.append(mut1, newmuts1, axis=0) #append effect of new mutations to mutation list
		
		pop2 = np.append(off2, np.transpose(np.identity(len(off2), dtype=int)[whomuts2[0]]), axis=1) #add new loci and identify mutants
		mut2 = np.append(mut2, newmuts2, axis=0) #append effect of new mutations to mutation list

		# remove lost mutations (all zero columns in pop)
		if remove_lost:
			if remove == 'any':
				keep = pop1.any(axis=0)
				mut1 = mut1[keep]
				pop1 = pop1[:, keep]
			if remove == 'derived':
				segregating = pop1.any(axis=0)
				ancestral = np.array(range(len(mut1))) < len(mutfound1)
				keep = np.add(segregating, ancestral)
				mut1 = mut1[keep]
				pop1 = pop1[:, keep]

		if remove_lost:
			if remove == 'any':
				keep = pop2.any(axis=0)
				mut2 = mut2[keep]
				pop2 = pop2[:, keep]
			if remove == 'derived':
				segregating = pop2.any(axis=0)
				ancestral = np.array(range(len(mut2))) < len(mutfound2)
				keep = np.add(segregating,ancestral)
				mut2 = mut2[keep]
				pop2 = pop2[:, keep]
		
		# go to next generation
		gen += 1

######################################################################
##RUNNING ADAPTATION FUNCTION##
######################################################################    
	
start = time.time()
main()
end = time.time()
print('this took %.2f seconds to complete' %(end-start))
