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

def open_output_files(n, K, alpha, B, u, data_dir):
	"""
	This function opens the output files and returns file
	handles to each.
	"""
	sim_id = 'n%d_K%d_alpha%.1f_B%d_u%.3f' %(n, K, alpha, B, u)
	outfile_A = open("%s/Movie_data_%s.csv" %(data_dir, sim_id), "w")
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

def Wrv(slist):
	"""
	Frequency of deleterious alleles randomly chosen from Wright's distribution without mutation given list of selection coefficients (uses inverse transform sampling)
	"""
	x = np.random.uniform(size = len(slist)) #random uniform numbers in [0,1]
	Ne = 2 * B / ( 2 * B - 1 ) * K #effective population size (eqn 13 in Burger & Lynch 1995 Evol)
	return [math.log( 1 - x[i] + np.exp(4 * Ne * slist[i]) * x[i] ) / (4 * Ne * slist[i] ) for i in range(len(x))]

######################################################################
##UNIVERSAL PARAMETERS##
######################################################################

nreps = 1 #number of replicates for each set of parameters
n = 2 #phenotypic dimensions (positive integer >=1)
data_dir = 'data'

######################################################################
##PARAMETERS TO MAKE ANCESTOR##
######################################################################

K = 50 #number of individuals (positive integer >=1)
n_mut_list = list([0]) #number of mutations in ancestor
alpha = 0.01 #mutational sd (positive real number)

######################################################################
##PARAMETERS FOR ADAPTING POPULATION##
######################################################################

K_adapt = K #number of individuals (positive integer)
alpha_adapt = alpha #mutational sd (positive real number)
B = 2 #number of offspring per generation per parent (positive integer)
u = 1 #mutation probability per generation per genome (0<u<1)

theta_list = np.array([[0.5,0.5]]) #optimum phenotypes for population

maxgen = 1000 #total number of generations population adapts for
gen_rec = 25 #print and save after this many generations

remove_lost = True #If true, remove mutations that are lost (0 for all individuals)
remove = 'derived' #.. any derived (not from ancestor) mutation that is lost 

######################################################################
##FUNCTION FOR POPULATION TO ADAPT##
######################################################################

def main():

	# open output files
	fileHandles = open_output_files(n, K, alpha, B, u, data_dir) 

	#loop over optima
	j = 0
	while j < len(theta_list):
		
		#set optimum
		theta = theta_list[j]
			
		#loop over all n_muts values
		i = 0
		while i < len(n_mut_list):

			n_muts = n_mut_list[i] #set number of mutations in ancestor (ie how much SGV)

			#loop over all replicates
			rep = 0
			while rep < nreps:

				#make ancestor
				if n_muts > 0:
					mut = np.random.normal(0, alpha, (n_muts, n)) #create n_muts mutations, each with a random normal phenotypic effect in each n dimension with mean 0 and sd alpha
					slist = survival(np.linalg.norm(mut, axis=1)) - 1 #selection coefficients for each mutant (in isolation in a background at the optimum)
					p_mut = Wrv(slist) #frequency of each mutation (randomly chosen from Wright's distribution, without mutation)
					pop = np.random.binomial(1, p_mut, (K, n_muts)) #randomly assign each of n_muts mutations to each of K individuals, weighted by random freqeucny from Wright's distribution
				else: #de novo only, even if p_mut>0
					pop = np.array([[1]] * K) #start all individuals with one "mutation"
					mut = np.array([[0] * n]) #but let the "mutation" do nothing (ie stay at origin)

				#found adapting population
				[popfound, mutfound] = found(K, K_adapt, pop, mut, remove_lost, remove)

				#initialize adapting population
				[pop, mut] = [popfound, mutfound]

				#intitialize generation counter
				gen = 0

				#run until maxgen
				while gen < maxgen + 1:

					# genotype to phenotype
					phenos = np.dot(pop, mut) #sum mutations held by each individual

					# viability selection
					surv = viability(phenos, theta, pop, K_adapt)

					#end simulation if population extinct (or unable to produce offspring)        
					if len(surv) < 2: 
						print("Extinct")              
						break 
						
					# meiosis
					off = recomb(surv, B)

					# mutation and population update
					[pop, mut] = mutate(off, u, alpha, n, mut)

					# remove lost mutations (all zero columns in pop)
					[pop, mut] = remove_muts(remove, remove_lost, pop, mut, mutfound)

					if gen % gen_rec == 0:
						
						p = np.dot(pop, mut) #phenotypes
						w = survival(np.linalg.norm(p - theta, axis=1)) #fitnesses


						#print update
						print('gen=%d, mean fitness=%.2f' %(gen, np.mean(w))) 
						
						#save data
						for x in range(len(p)):
							write_data_to_output(fileHandles, [gen] + p[x].tolist() + [w[x]]) #gen, phenotype in each dimension, fitness
						
					# go to next generation
					gen += 1

				# go to next rep
				rep += 1

			#go to next n_muts value
			i += 1

		#go to next optima
		j += 1

	# cleanup
	close_output_files(fileHandles)

######################################################################
##RUNNING ADAPTATION FUNCTION##
######################################################################    
	
start = time.time()
main()
end = time.time()
print('this took %.2f seconds to complete' %(end-start))
