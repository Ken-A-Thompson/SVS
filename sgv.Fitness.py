#Authors: Matthew Osmond <mmosmond@zoology.ubc.ca> & Ken A. Thompson <ken.thompson@zoology.ubc.ca>
#Description: The role of standing genetic variance in speciation

import numpy as np
import time
import matplotlib.pyplot as plt
import csv
import math

######################################################################
##FUNCTIONS##
######################################################################

def open_output_files(n, K, p_mut, alpha, B, u, data_dir):
	"""
	This function opens the output files and returns file
	handles to each.
	"""
	sim_id = 'n%d_K%d_pmut%.1f_alpha%.1f_B%d_u%.3f' %(n, K, p_mut, alpha, B, u)
	outfile_A = open("%s/fitness_%s.csv" %(data_dir, sim_id), "w")
	return outfile_A

def write_data_to_output(fileHandles, data):
	"""
	This function writes a (time, data) pair to the
	corresponding output file. We write densities
	not abundances.
	"""
	writer = csv.writer(fileHandles)
	writer.writerow(data)

def close_output_files(fileHandles):
	"""
	This function closes all output files.
	"""
	fileHandles.close()

def found(K, K_adapt, pop, mut, remove_lost, remove):
	"""
	This function creates a founding population from an ancestral one
	"""
	whofounds = np.random.choice(K, size = K_adapt, replace = False) #random choice of K_adapt founders from ancestral population 
	popfound = pop[whofounds] #list of mutations held by each founding individual
	if remove_lost and remove == 'any': #if removing ancestral mutations when lost
		keep = popfound.any(axis = 0) #loci that have mutation in at least one founding individual
		mutfound = mut[keep] #keep those mutations
		popfound = pop[:, keep] #keep those loci
	else:
		mutfound = mut #else just keep all mutations (and all loci)
	return [popfound, mutfound]

def survival(dist):
	"""
	This function gives the probability of survival
	"""
	return np.exp(-0.5 * dist**2) #probability of survival

def viability(phenos, theta, pop, K_adapt):
	"""
	This function determines which individuals survive viability selection
	"""
	dist = np.linalg.norm(phenos - theta, axis=1) #phenotypic distance from optimum
	w = survival(dist) #probability of survival
	rand = np.random.uniform(size = len(pop)) #random uniform number in [0,1] for each individual
	surv = pop[rand < w] #survivors
	if len(surv) > K_adapt:
		surv = surv[np.random.randint(len(surv), size = K_adapt)] #randomly choose K_adapt individuals if more than K_adapt
	return surv

def recomb(surv, B):
	"""
	This function creates offspring through pairing of parents (haploid) and recombination (i.e, meiosis)
	"""
	pairs = np.resize(np.random.choice(len(surv), size=len(surv), replace=False), (int(len(surv)/2), 2)) #random mate pairs (each mates at most once and not with self)
	rand2 = np.random.randint(2, size=(len(pairs), len(surv[0]))) #from which parent each offspring inherits each allele (free recombination, fair transmission)
	rec = np.resize(np.append(rand2, 1-rand2, axis=1),(len(rand2), 2, len(rand2[0]))) #reshape
	off_1 = np.sum(surv[pairs] * rec, axis=1) #one product of meiosis
	off_2 = np.sum(surv[pairs] * (1-rec), axis=1) #other product of meiosis
	off = np.repeat(np.append(off_1, off_2, axis=0), B, axis=0) #each product of meiosis produced B times
	return off

def mutate(off, u, alpha, n, mut):
	"""
	This function creates mutations and updates population
	"""
	rand3 = np.random.uniform(size = len(off)) #random uniform number in [0,1] for each offspring
	nmuts = sum(rand3 < u) # mutate if random number is below mutation rate; returns number of new mutations
	whomuts = np.where(rand3 < u) #indices of mutants
	newmuts = np.random.normal(0, alpha, (nmuts, n)) #phenotypic effect of new mutations
	pop = np.append(off, np.transpose(np.identity(len(off), dtype=int)[whomuts[0]]), axis=1) #add new loci and identify mutants
	mut = np.append(mut, newmuts, axis=0) #append effect of new mutations to mutation list
	return [pop, mut]

def remove_muts(remove, remove_lost, pop, mut, mutfound):
	"""
	This function creates mutations and updates population
	"""
	if remove_lost:
		if remove == 'any':
			keep = pop.any(axis=0)
			mut = mut[keep]
			pop = pop[:, keep]
		elif remove == 'derived':
			segregating = pop.any(axis=0)
			ancestral = np.array(range(len(mut))) < len(mutfound)
			keep = np.add(segregating, ancestral)
			mut = mut[keep]
			pop = pop[:, keep]
	return [pop, mut]

######################################################################
##UNIVERSAL PARAMETERS##
######################################################################

nreps = 10 #number of replicates for each set of parameters
n = 2 #phenotypic dimensions (positive integer >=1); actually need it greater or equal to 2 for angle code below to work (otherwise just 0 or 180 valid) 
data_dir = 'data'

######################################################################
##PARAMETERS TO MAKE ANCESTOR##
######################################################################

K = 1000 #number of individuals (positive integer >=1)
n_mut_list = list(np.arange(0, 31, 30)) #start, end+1, interval number >=0)
p_mut = 0.1 #probability of having mutation at any one locus (0<=p<=1)
alpha = 0.1 #mutational sd (positive real number)

######################################################################
##PARAMETERS FOR ADAPTING POPULATIONS##
######################################################################

K_adapt = K #number of individuals (positive integer)
alpha_adapt = alpha #mutational sd (positive real number)
B = 2 #number of offspring per generation per parent (positive integer)
u = 0.001 #mutation probability per generation per genome (0<u<1)

# opt_angles = [0] #angles of 0, 60 and 180 degrees, which give distances between optima of 0, d, and 2d
opt_angles = [0, np.pi/3, np.pi] #angles of 0, 60 and 180 degrees, which give distances between optima of 0, d, and 2d
opt_dists = list(np.arange(0, 0.801, 0.05)) #distances from origin to optima (d)

maxgen = 2000 #total number of generations populations adapt for

remove_lost = True #If true, remove mutations that are lost (0 for all individuals)
remove = 'derived' #.. any derived (not from ancestor) mutation that is lost 

######################################################################
##PARAMETERS FOR HYBRIDS##
######################################################################

nHybrids = 100 #number of hybrids to make at end of each replicate

######################################################################
##FUNCTION FOR POPULATIONS TO ADAPT##
######################################################################

def main():

	# open output files
	fileHandles = open_output_files(n, K, p_mut, alpha, B, u, data_dir) 

	#loop over angle between optima
	k = 0
	while k < len(opt_angles):

		#loop over adaptive distances
		j = 0
		while j < len(opt_dists):
			
			#set optima
			theta1 = np.append(opt_dists[j],[0]*(n-1)) #set one optima to be on axis 1
			theta2 = np.append([opt_dists[j] * math.cos(opt_angles[k]), opt_dists[j] * math.sin(opt_angles[k])], [0]*(n-2)) #set other to be on plane of 1st and 2nd axes (note this requires n>=2) 

			#loop over all n_muts values
			i = 0
			while i < len(n_mut_list):

				n_muts = n_mut_list[i] #set number of mutations in ancestor (ie how much SGV)

				#loop over all replicates
				rep = 0
				while rep < nreps:

					#make ancestor
					if n_muts > 0:
						pop = np.random.binomial(1, p_mut, (K, n_muts)) #p_mut chance of having each of n_muts mutations, for all K individuals
						mut = np.random.normal(0, alpha, (n_muts, n)) #create n_muts mutations, each with a random normal phenotypic effect in each n dimension with mean 0 and sd alpha
					else: #de novo only, even if p_mut>0
						pop = np.array([[1]] * K)
						mut = np.array([[0] * n])

					#found adapting populations
					[popfound1, mutfound1] = found(K, K_adapt, pop, mut, remove_lost, remove)
					[popfound2, mutfound2] = found(K, K_adapt, pop, mut, remove_lost, remove)

					#initialize adapting populations
					[pop1, mut1] = [popfound1, mutfound1]
					[pop2, mut2] = [popfound2, mutfound2]

					#intitialize generation counter
					gen = 0

					#run until maxgen
					while gen < maxgen + 1:

						# genotype to phenotype
						phenos1 = np.dot(pop1, mut1) #sum mutations held by each individual
						phenos2 = np.dot(pop2, mut2) #sum mutations held by each individual

						# viability selection
						surv1 = viability(phenos1, theta1, pop1, K_adapt)
						surv2 = viability(phenos2, theta2, pop2, K_adapt)

						#end simulation if either population extinct (or unable to produce offspring)        
						if len(surv1) < 2 or len(surv2) < 2: 
							print("Extinct")              
							break 
							
						# meiosis
						off1 = recomb(surv1, B)
						off2 = recomb(surv2, B)

						# mutation and population update
						[pop1, mut1] = mutate(off1, u, alpha, n, mut1)
						[pop2, mut2] = mutate(off2, u, alpha, n, mut2)

						# remove lost mutations (all zero columns in pop)
						[pop1, mut1] = remove_muts(remove, remove_lost, pop1, mut1, mutfound1)
						[pop2, mut2] = remove_muts(remove, remove_lost, pop2, mut2, mutfound2)

						# go to next generation
						gen += 1

					#parent fitness and load (use parent 1, but could be either)	
					parents = np.random.randint(len(pop1), size = nHybrids)
					parent_phenos = np.dot(pop1[parents],mut1)	
					dist = np.linalg.norm(parent_phenos - np.mean(parent_phenos, axis=0), axis=1) #phenotypic distance from mean
					fitness = survival(dist) #viability
					pfit = np.mean(fitness) #mean fitness
					growth = np.log(fitness*B) #continuous time growth rates
					pload = np.log(1*B) - np.mean(growth) #segregation load
					
					#make variables to hold hybrid phenotypes
					offphenos = dict()
					offpheno = []

					#make each of nHybrids hybrids
					for m in range(nHybrids):
					    # choose random parents
						randpar1 = pop1[np.random.choice(len(pop1))] 
						randpar2 = pop2[np.random.choice(len(pop2))]
						# get random parent phenotypes
						phenpar1 = np.dot(randpar1, mut1) 
						phenpar2 = np.dot(randpar2, mut2)
						# get mutations held by random parents
						mutpar1 = mut1 * randpar1[:, None]
						mutpar2 = mut2 * randpar2[:, None]
						setA = set(tuple(x) for x in mutpar1)
						setB = set(tuple(x) for x in mutpar2)
						# find mutations shared by two parents (all in offspring)
						sharedmuts = np.array([x for x in setA & setB])
						if len(sharedmuts) < 1:
							sharedmuts = np.array([[0] * n]) #give something in case empty
						# find mutations not shared by two parents
						unsharedmuts = np.array([x for x in setA ^ setB])
						# which unshared mutations in offspring (free recombination between all loci, therefore gets each with 0.5 probability)
						randmuts = np.random.randint(2, size = (len(unsharedmuts)))	
						unsharedoffmuts = unsharedmuts * randmuts[:, None]
						if len(unsharedoffmuts) < 1:
						    unsharedoffmuts = np.array([[0] * n]) #give something in case empty
						# offspring phenotype is collection of shared and random unshared mutations
						offpheno.append(sum(np.append(sharedmuts, unsharedoffmuts, axis = 0)))

					#hybrid load
					offpheno = np.array(offpheno) #reformat correctly
					dist = np.linalg.norm(offpheno - np.mean(offpheno, axis=0), axis=1) #phenotypic distance from mean hybrid
					hyload = np.log(1*B) - np.mean(np.log(survival(dist)*B)) #segregation load of hybrids

					#mean hybrid fitness
					dist1 = np.linalg.norm(offpheno - theta1, axis=1) #phenotypic distance from parental 1 optimum
					dist2 = np.linalg.norm(offpheno - theta2, axis=1) #phenotypic distance from parental 2 optimum
					dist = np.minimum(dist1, dist2) #distance to closest optimum
					fitness = survival(dist) #viability
					hyfit = np.mean(fitness) #mean fitness
					rhyfit = hyfit/pfit #relative fitness

					#print an update
					print('rep=%d, n_muts=%d, dist=%.1f, angle=%.1f, parent load=%.3f, hybrid load=%.3f, parent fitness=%.3f, relative hybrid fitness=%.3f' %(rep+1, n_mut_list[i], opt_dists[j], opt_angles[k]*180/np.pi, pload, hyload, pfit, rhyfit)) 
					
					#save data
					write_data_to_output(fileHandles, [rep+1, n_mut_list[i], opt_dists[j], opt_angles[k]*180/np.pi, pload, hyload, pfit, rhyfit])

					# go to next rep
					rep += 1

				#go to next n_muts value
				i += 1

			#go to next optima
			j += 1

		#go to next type of selection
		k += 1

	# cleanup
	close_output_files(fileHandles)

######################################################################
##RUNNING ADAPTATION FUNCTION##
######################################################################    
	
start = time.time()
main()
end = time.time()
print('this took %.2f seconds to complete' %(end-start))
