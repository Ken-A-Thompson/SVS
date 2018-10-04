#Authors: Matthew Osmond <mmosmond@zoology.ubc.ca> & Ken A. Thompson <ken.thompson@zoology.ubc.ca>
#Description: Dynamics of parent populations

import numpy as np
import time
import matplotlib.pyplot as plt
import csv

######################################################################
##FUNCTIONS##
######################################################################

def open_output_files(n, N, alpha, u, sigma, data_dir):
	"""
	This function opens the output files and returns file
	handles to each.
	"""
	sim_id = 'n%d_N%d_alpha%.4f_u%.4f_sigma%.4f' %(n, N, alpha, u, sigma)
	outfile_A = open("%s/FigS2_%s.csv" %(data_dir, sim_id), "w")
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

def found(n_muts, ancestor_muts, ancestor_freqs, N_adapt, n):
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
		popfound = np.random.binomial(1, p_mut, (N_adapt, n_muts)) #p_mut chance of having each of n_muts mutations, for all K_adapt individuals
		fix_id = [ancestor_freqs == 1] #indices for fixed mutations in ancestor
		mutfound = np.append(mutfound, ancestor_muts[fix_id], axis=0) #add fixed mutations to founding mutation matrix
		addpop = np.array([1]*N_adapt*np.sum(fix_id)).reshape(N_adapt,np.sum(fix_id)) #matrix of 1s for fixed mutations
		popfound = np.append(popfound, addpop, axis=1) #add fixed mutations to founding pop matrix
	else: #de novo only, even if p_mut>0
		popfound = np.array([[1]] * N_adapt)
		mutfound = np.array([[0] * n])
	return [popfound, mutfound]

def fitness(phenos, theta, sigma):
	"""
	This function determines relative fitness
	"""
	dist = np.linalg.norm(phenos - theta, axis=1) #phenotypic distance from optimum
	w = np.exp(-0.5 * sigma * dist**2) #fitness
	return w

def recomb(surv):
	"""
	This function creates offspring through pairing of parents (haploid) and recombination (i.e, meiosis)
	"""
	pairs = np.resize(np.random.choice(len(surv), size=len(surv), replace=False), (int(len(surv)/2), 2)) #random mate pairs (each mates at most once and not with self)
	rand2 = np.random.randint(2, size=(len(pairs), len(surv[0]))) #from which parent each offspring inherits each allele (free recombination, fair transmission)
	rec = np.resize(np.append(rand2, 1-rand2, axis=1),(len(rand2), 2, len(rand2[0]))) #reshape
	off_1 = np.sum(surv[pairs] * rec, axis=1) #one product of meiosis
	off_2 = np.sum(surv[pairs] * (1-rec), axis=1) #other product of meiosis
	off = np.append(off_1, off_2, axis=0) #each product of meiosis
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

nreps = 2 #number of replicates for each set of parameters
ns = [2, 5, 10] #phenotypic dimensions (positive integer >=1)
data_dir = 'data'

######################################################################
##PARAMETERS OF ANCESTOR##
######################################################################

n_reps = 2 #number of reps of ancestor that exist
N = 10**4 #number of individuals (positive integer >=1)
alpha = 2*10**(-1) #mutational sd (positive real number)
u = 10**(-3) #mutation probability per generation per genome (0<u<1)
sigma = 10**(-2) #selection strength
burn_dir = 'data'
rrep = np.random.choice(n_reps, nreps, replace=True) #randomly assign each rep an ancestor, without or with replacement (i.e., unique ancestor for each sim or not)

######################################################################
##PARAMETERS FOR ADAPTING POPULATION##
######################################################################

n_mut_list = [list(np.arange(0, 2, 1)), list(np.arange(0, 2, 1)), list(np.arange(0, 2, 1))] #starting nmuts, final n_muts, interval

N_adapts = [10**3, 10**4] #number of haploid individuals (positive integer)
alpha_adapt = alpha #mutational sd (positive real number)
u_adapt = 0.001 #mutation probability per generation per genome (0<u<1)
sigma_adapts = [10**(-1), 10**0] #selection strengths

maxgen = 20 #total number of generations population adapts for
gen_rec = 5 #print and save after this many generations

remove_lost = True #If true, remove mutations that are lost (0 for all individuals)
remove = 'derived' #.. any derived (not from ancestor) mutation that is lost 

######################################################################
##FUNCTION FOR POPULATION TO ADAPT##
######################################################################

def main():

	#loop over population size
	i_N = 0
	while i_N < len(N_adapts):
		N_adapt = N_adapts[i_N]

		#loop over selection strength
		i_sigma = 0
		while i_sigma < len(sigma_adapts):
			sigma_adapt = sigma_adapts[i_sigma]

			#loop over dimensions
			l = 0
			while l < len(ns):
				n = ns[l]

				#set n-dependent parameters
				theta_list = [np.append(1,[0]*(n-1))] #optimum

				# open output files
				fileHandles = open_output_files(n, N_adapt, alpha_adapt, u_adapt, sigma_adapt, data_dir) 

				#loop over optima
				j = 0
				while j < len(theta_list):
					
					#set optimum
					theta = theta_list[j]
						
					#loop over all n_muts values
					i = 0
					while i < len(n_mut_list[l]):

						n_muts = n_mut_list[l][i] #set number of mutations in ancestor (ie how much SGV)

						#loop over all replicates
						rep = 0
						while rep < nreps:

							#load ancestor
							burn_id = 'n%d_N%d_alpha%.4f_u%.4f_sigma%.4f_rep%d' %(n, N, alpha, u, sigma, rrep[rep]+1)

							filename = "%s/Muts_%s.npy" %(burn_dir, burn_id)
							ancestor_muts = np.load(filename) #load mutations

							filename = "%s/Freqs_%s.npy" %(burn_dir, burn_id)
							ancestor_freqs = np.load(filename) #load frequencies

							#found adapting population
							[popfound, mutfound] = found(n_muts, ancestor_muts, ancestor_freqs, N_adapt, n)

							#initialize adapting population
							[pop, mut] = [popfound, mutfound]

							#intitialize generation counter
							gen = 0

							#run until maxgen
							while gen < maxgen + 1:

								# genotype to phenotype
								phenos = np.dot(pop, mut) #sum mutations held by each individual

								# phenotype to fitness
								w = fitness(phenos, theta, sigma_adapt)

								# wright-fisher (multinomial) sampling
								parents = np.random.multinomial(N_adapt, w/sum(w)) #number of times each parent chosen
								off = np.repeat(pop, parents, axis=0) #offspring genotypes

								# mating and recombination
								off = recomb(off)

								# mutation and population update
								[pop, mut] = mutate(off, u_adapt, alpha_adapt, n, mut)

								# remove lost mutations (all zero columns in pop)
								[pop, mut] = remove_muts(remove, remove_lost, pop, mut, mutfound)

								if gen % gen_rec == 0:
									
									#calculate number of segregating sites
									notlost = pop.any(axis=0) #sites not lost
									fixed = -1*(pop-1)
									notfixed = fixed.any(axis=0) #sites not fixed
									segregating = notlost*notfixed #sites not lost or fixed (segregating)
									numseg = sum(segregating) #number of segregating sites

									mean_dist = np.linalg.norm(np.mean(phenos, axis=0) - theta, axis=0) #mean euclidean distance to optimum

									#parent fitness and load (use parent 1, but could be either)	
									parents = np.random.randint(len(pop), size = N_adapt)
									parent_phenos = np.dot(pop[parents], mut)	
									# mean_parent_pheno = np.mean(parent_phenos, axis=0)
									# parent_fitnesses = fitness(parent_phenos, mean_parent_pheno, sigma_adapt) #parent fitnesses
									psegvar = np.mean(np.var(parent_phenos, axis = 0)) # segregation variance (mean of individual trait variances)

									#print update
									print('N=%d, sigma=%.2f, n=%d, opt=%r, n_muts=%d, rep=%d, gen=%d, segregating sites=%d, mean distance to opt = %.3f, segregation variance =%.3f' %(N_adapt, sigma_adapt, n, [round(x,2) for x in theta], n_muts,  rep+1, gen, numseg, mean_dist, psegvar)) 
									
									#save data
									write_data_to_output(fileHandles, [theta, n_muts, rep+1,  gen, numseg, mean_dist, psegvar])
									
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

				#next dimension
				l += 1

			#go to next sigma value
			i_sigma += 1

		#go to next N value
		i_N += 1

######################################################################
##RUNNING ADAPTATION FUNCTION##
######################################################################    
	
start = time.time()
main()
end = time.time()
print('this took %.2f seconds to complete' %(end-start))
