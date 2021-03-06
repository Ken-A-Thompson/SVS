#Authors: Matthew Osmond <mmosmond@zoology.ubc.ca> & Ken A. Thompson <ken.thompson@zoology.ubc.ca>
#Description: The role of standing genetic variance in speciation

import numpy as np
import time
import csv
import random
# import matplotlib.pyplot as plt


######################################################################
##FUNCTIONS##
######################################################################

def open_output_files(n, K, alpha, B, u, data_dir):
	"""
	This function opens the output files and returns file
	handles to each.
	"""
	sim_id = 'n%d_K%d_alpha%.1f_B%d_u%.3f' %(n, K_adapt, alpha, B, u)
	outfile_A = open("%s/hybrid_loads_%s.csv" %(data_dir, sim_id), "w")
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

def found(n_muts, ancestor_muts, ancestor_freqs, K_adapt, n):
	"""
	This function creates a founding population from an ancestral one
	"""

	#make ancestor
	if n_muts > 0:
		seg_id = [ancestor_freqs < 1] #indices for segregating mutations in ancestor
		nmuts_max = np.sum(seg_id) #number of segregating mutations in ancestor
		probs = ancestor_freqs[seg_id]/ancestor_freqs[seg_id].sum() #probability of choosing each mutation in sgv (make pdf)
		mut_choice = np.random.choice(nmuts_max, size=n_muts, replace=False, p=probs) #indices of mutations to take from ancestor
		mutfound = (ancestor_muts[seg_id])[mut_choice] #mutational effects
		p_mut = (ancestor_freqs[seg_id])[mut_choice] #expected frequency of these mutations
		popfound = np.random.binomial(1, p_mut, (K_adapt, n_muts)) #p_mut chance of having each of n_muts mutations, for all K_adapt individuals
		fix_id = [ancestor_freqs == 1] #indices for fixed mutations in ancestor
		mutfound = np.append(mutfound, ancestor_muts[fix_id], axis=0) #add fixed mutations to founding mutation matrix
		addpop = np.array([1]*K_adapt*np.sum(fix_id)).reshape(K_adapt,np.sum(fix_id)) #matrix of 1s for fixed mutations
		popfound = np.append(popfound, addpop, axis=1) #add fixed mutations to founding pop matrix
	else: #de novo only, even if p_mut>0
		popfound = np.array([[1]] * K_adapt)
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
sigma = 0.01 #selection strength

burn_dir = 'data/burnins_jun25'
rrep = np.random.choice(n_reps, nreps, replace=False) #randomly assign each rep an ancestor, without or with replacement (i.e., unique ancestor for each sim or not)

######################################################################
##PARAMETERS FOR ADAPTING POPULATIONS##
######################################################################

n_mut_list = list(np.arange(0, 101, 5)) #starting nmuts, final n_muts, interval

K_adapt = 1000 #number of individuals (positive integer)
alpha_adapt = alpha #mutational sd (positive real number)
B_adapt = B #number of offspring per generation per parent (positive integer)
u_adapt = u #mutation probability per generation per genome (0<u<1)

opt_dists = list(np.arange(1, 1.01, 0.1)) #distances to optima

# selection = 'divergent' #divergent selection (angle = 180 deg)
# selection = 'parallel' #parallel selection (angle = 0)
selection = 'both' #both divergent and parallel selection

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
	fileHandles = open_output_files(n, K_adapt, alpha_adapt, B_adapt, u_adapt, data_dir) 

	if selection == 'both':
		k = 0
		kmax = 1
	elif selection == 'parallel':
		k = 0
		kmax = 0
	elif selection == 'divergent':
		k = 1
		kmax = 1

	#loop of selection styles
	while k < kmax + 1:

		#loop over optima
		j = 0
		while j < len(opt_dists):
			
			#set optima
			theta1 = np.append(opt_dists[j],[0]*(n-1)) #set one optima

			if k == 0: #parallel
				theta2 = theta1
			elif k == 1: #divergent
				theta2 = np.append(-opt_dists[j],[0]*(n-1))
				
			# #set up plot of hybrid load versus number of ancestral mutations (n_muts)
			# plt.axis([0, max(n_mut_list)+1, 0, 0.1])
			# plt.ylabel('hybrid load at generation %d (mean $\pm$ SD of %d replicates)' %(maxgen,nreps))
			# plt.xlabel('number of ancestral mutations')
			# plt.ion()

			#loop over all n_muts values
			i = 0
			while i < len(n_mut_list):

				n_muts = n_mut_list[i] #set number of mutations in ancestor (ie how much SGV)

				# hyloads = [0] * nreps #initialize vector to store hybrid loads in from each replicate

				#loop over all replicates
				rep = 0
				while rep < nreps:

					#load ancestor
					burn_id = 'n%d_K%d_alpha%.1f_B%d_u%.4f_sigma%.2f_rep%d' %(n, K, alpha, B, u, sigma, rrep[rep]+1)

					filename = "%s/Muts_%s.npy" %(burn_dir, burn_id)
					ancestor_muts = np.load(filename) #load mutations

					filename = "%s/Freqs_%s.npy" %(burn_dir, burn_id)
					ancestor_freqs = np.load(filename) #load frequencies

					#found adapting populations
					# [popfound1, mutfound1] = found(n_muts, nmuts_max, ancestor_muts, ancestor_freqs, K, n)
					# [popfound2, mutfound2] = found(n_muts, nmuts_max, ancestor_muts, ancestor_freqs, K, n)

					#initialize adapting populations
					# [pop1, mut1] = [popfound1, mutfound1]
					# [pop2, mut2] = [popfound2, mutfound2]

					#found identical populations
					[popfound, mutfound] = found(n_muts, ancestor_muts, ancestor_freqs, K_adapt, n)
					[pop1, mut1] = [popfound, mutfound]
					[pop2, mut2] = [popfound, mutfound]

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
						[pop1, mut1] = remove_muts(remove, remove_lost, pop1, mut1, mutfound)
						[pop2, mut2] = remove_muts(remove, remove_lost, pop2, mut2, mutfound)

						# go to next generation
						gen += 1

					#make variables to hold offspring phenotypes
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

					offpheno = np.array(offpheno) #reformat correctly
					dist = np.linalg.norm(offpheno - np.mean(offpheno, axis=0), axis=1) #phenotypic distance from mean hybrid
					# hyload = np.log(1*B) - np.mean(np.log(survival(dist)*B)) #hybrid load as defined by Chevin et al 2014
					segvar = np.mean(np.var(offpheno, axis = 0))
					
					#print an update
					print('opt1=%r, opt2=%r, rep=%d, n_muts=%d, segregation variance=%.3f, distance=%.3f, selection=%r' %(theta1, theta2, rep+1, n_mut_list[i], segvar, opt_dists[j], ['parallel','divergent'][k])) 
					
					#save data
					write_data_to_output(fileHandles, [theta1, theta2, rep+1, n_mut_list[i], segvar, opt_dists[j], ['parallel','divergent'][k]])

					# hyloads[rep] = hyload #save hybrid load for this replicate

					# go to next rep
					rep += 1

				# #plot mean and SD hybrid load over all replicates for this n_muts value
				# plt.errorbar(n_mut_list[i], np.mean(hyloads), yerr=np.var(hyloads)**0.5, fmt='o', color='k')
				# plt.pause(0.01)

				#go to next n_muts value
				i += 1

			# plt.pause(1) #pause on finished plot for a second
			# plt.savefig('Figs/HLvsNMUT.png') #save finished plot

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
