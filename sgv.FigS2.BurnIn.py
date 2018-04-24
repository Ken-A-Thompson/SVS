#Authors: Matthew Osmond <mmosmond@zoology.ubc.ca> & Ken A. Thompson <ken.thompson@zoology.ubc.ca>
#Description: The role of standing genetic variance in speciation

import numpy as np
import time
import matplotlib.pyplot as plt
import csv
import math
import ast

######################################################################
##FUNCTIONS##
######################################################################

def open_output_files(n, K, alpha, B, u, data_dir):
	"""
	This function opens the output files and returns file
	handles to each.
	"""
	sim_id = 'n%d_K%d_alpha%.1f_B%d_u%.3f' %(n, K, alpha, B, u)
	outfile_A = open("%s/segregation_load_%s.csv" %(data_dir, sim_id), "w")
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

def found(n_muts, nmuts_max, ancestor_muts, ancestor_freqs, K, n):
	"""
	This function creates a founding population from an ancestral one
	"""

	#make ancestor
	if n_muts > 0:
		mut_choice = np.random.choice(nmuts_max, size=n_muts, replace=False) #indices of mutations to take from ancestor
		mutfound = ancestor_muts[mut_choice] #mutational effects
		p_mut = ancestor_freqs[mut_choice] #expected frequency of these mutations
		popfound = np.random.binomial(1, p_mut, (K, n_muts)) #p_mut chance of having each of n_muts mutations, for all K individuals
	else: #de novo only, even if p_mut>0
		popfound = np.array([[1]] * K)
		mutfound = np.array([[0] * n])
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
n = 2 #phenotypic dimensions (positive integer >=1)
data_dir = 'data'

######################################################################
##PARAMETERS OF ANCESTOR##
######################################################################

n_reps = 10 #number of reps of ancestor that exist
K = 10000 #number of individuals (positive integer >=1)
alpha = 0.1 #mutational sd (positive real number)
B = 2 #number of offspring per generation per parent (positive integer)
u = 0.001 #mutation probability per generation per genome (0<u<1)
sigma = 0.1 #selection strength
burn_dir = 'data/burnins'
rrep = np.random.choice(n_reps, nreps, replace=False) #randomly assign each rep an ancestor

######################################################################
##PARAMETERS FOR ADAPTING POPULATION##
######################################################################

K_adapt = K #number of individuals (positive integer)
alpha_adapt = alpha #mutational sd (positive real number)
B = B #number of offspring per generation per parent (positive integer)
u = u #mutation probability per generation per genome (0<u<1)

opt_dists = list(np.arange(0.2, 1.01, 0.4)) #distances to optima

n_mut_list = list(np.arange(1, 2, 1))

maxgen = 200 #total number of generations population adapts for
gen_rec = 100 #print and save after this many generations

remove_lost = True #If true, remove mutations that are lost (0 for all individuals)
remove = 'derived' #.. any derived (not from ancestor) mutation that is lost 

nHybrids = 100 #number of offspring to make to calculate load from

######################################################################
##FUNCTION FOR POPULATION TO ADAPT##
######################################################################

def main():

	# open output files
	fileHandles = open_output_files(n, K, alpha, B, u, data_dir) 

	#loop over optima
	j = 0
	while j < len(opt_dists):
		
		#set optima
		theta = np.append(opt_dists[j],[0]*(n-1)) #set optimum
		
		#loop over all n_muts values
		i = 0
		while i < len(n_mut_list):

			n_muts = n_mut_list[i] #set number of mutations in ancestor (ie how much SGV)

			#loop over all replicates
			rep = 0
			while rep < nreps:

				burn_id = 'n%d_K%d_alpha%.1f_B%d_u%.4f_sigma%.1f_rep%d' %(n, K, alpha, B, u, sigma, rrep[rep]+1)

				with open('%s/BurnIn_%s.csv' %(burn_dir,burn_id)) as csvfile:
				    readCSV = csv.reader(csvfile, delimiter=',')
				    ancestor_muts = []
				    ancestor_freqs = []
				    for row in readCSV:
				        mut = ast.literal_eval(row[0].replace("[ ","[").replace("  "," ").replace(" ",",").replace(",,",",").replace(",,",","))
				        freq = float(row[1])
				        ancestor_muts.append(mut)
				        ancestor_freqs.append(freq)

				ancestor_muts = np.array(ancestor_muts)
				ancestor_freqs = np.array(ancestor_freqs)
				nmuts_max = len(ancestor_freqs) #number of mutations in ancestor

				#found adapting population
				[popfound, mutfound] = found(n_muts, nmuts_max, ancestor_muts, ancestor_freqs, K, n)

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
						
						#make variables to hold offspring phenotypes
						offphenos = dict()
						offpheno = []

						#make each of nHybrids hybrids
						for k in range(nHybrids):
						    # choose random parents
							randpar1 = pop[np.random.choice(len(pop))] 
							randpar2 = pop[np.random.choice(len(pop))]
							# get random parent phenotypes
							phenpar1 = np.dot(randpar1, mut) 
							phenpar2 = np.dot(randpar2, mut)
							# get mutations held by random parents
							mutpar1 = mut * randpar1[:, None]
							mutpar2 = mut * randpar2[:, None]
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

						offpheno = np.array(offpheno) #reformat correctly
						dist = np.linalg.norm(offpheno - np.mean(offpheno, axis=0), axis=1) #phenotypic distance from mean hybrid
						hyload = np.log(1*B) - np.mean(np.log(survival(dist)*B)) #hybrid load as defined by Chevin et al 2014
		
						#print update
						print('opt_dist=%.3f, n_muts=%d, rep=%d, gen=%d, segregation load=%.3f' %(opt_dists[j], n_muts,  rep+1, gen, hyload)) 
						
						#save data
						write_data_to_output(fileHandles, [opt_dists[j], n_muts, rep+1,  gen, hyload])
						
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