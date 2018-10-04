#Authors: Matthew Osmond <mmosmond@zoology.ubc.ca> & Ken A. Thompson <ken.thompson@zoology.ubc.ca>
#Description: The role of standing genetic variance in speciation

import numpy as np
import time
# import matplotlib.pyplot as plt
import csv
import math

######################################################################
##FUNCTIONS##
######################################################################

def open_output_files(n, N, alpha, u, sigma, p_mut, data_dir):
	"""
	This function opens the output files and returns file
	handles to each.
	"""
	sim_id = 'n%d_N%d_alpha%.4f_u%.4f_sigma%.4f_pmut%.4f' %(n, N, alpha, u, sigma, p_mut)
	outfile_A = open("%s/FigS4_%s.csv" %(data_dir, sim_id), "w")
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

def found(N, N_adapt, pop, mut, remove_lost, remove):
	"""
	This function creates a founding population from an ancestral one
	"""
	whofounds = np.random.choice(N, size = N_adapt, replace = False) #random choice of N_adapt founders from ancestral population 
	popfound = pop[whofounds] #list of mutations held by each founding individual
	if remove_lost and remove == 'any': #if removing ancestral mutations when lost
		keep = popfound.any(axis = 0) #loci that have mutation in at least one founding individual
		mutfound = mut[keep] #keep those mutations
		popfound = pop[:, keep] #keep those loci
	else:
		mutfound = mut #else just keep all mutations (and all loci)
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

nreps = 1 #number of replicates for each set of parameters
ns = [2] #phenotypic dimensions (positive integer >=1)
data_dir = 'data'

######################################################################
##PARAMETERS TO MAKE ANCESTOR##
######################################################################

N = 10**3 #number of individuals (positive integer >=1)
n_mut_list = [[0, 10]] #number of mutations (positive integer >=1)
p_mut = 10**(-2) #probability of having mutation at any one locus (0<=p<=1) #set this to zero for de novo only
alpha = 10**(-2) #mutational sd (positive real number)

######################################################################
##PARAMETERS FOR ADAPTING POPULATIONS##
######################################################################

N_adapts = [10**3] #number of haploid individuals (positive integer)
alpha_adapt = alpha #mutational sd (positive real number)
u_adapt = 10**(-3) #mutation probability per generation per genome (0<u<1)
sigma_adapts = [10**0] #selection strengths

opt_dist = 1 #distance to optima

n_angles = 2 #number of angles between optima to simulate (including 0 and 180)
angles = [math.pi*x/(n_angles-1) for x in range(n_angles)] #angles to use (in radians)

maxgen = 10**3 #total number of generations populations adapt for

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
				theta1 = np.append(opt_dist,[0]*(n-1)) #set one optima to be fixed
				if n == 2:
					theta2_list = np.array([[opt_dist*math.cos(x), opt_dist*math.sin(x)] for x in angles]) #optima to use
				elif n > 2:
					theta2_list = np.array([np.append([opt_dist*math.cos(x), opt_dist*math.sin(x)], [0]*(n-2)) for x in angles]) #optima to use
				
				# open output files
				fileHandles = open_output_files(n, N, alpha, u_adapt, sigma_adapt, p_mut, data_dir) 

				#loop over optima
				j = 0
				while j < len(theta2_list):
					
					#set optima
					theta2 = theta2_list[j]
						
					# #set up plot of hybrid load versus number of ancestral mutations (n_muts)
					# plt.axis([0, max(n_mut_list)+1, 0, 0.1])
					# plt.ylabel('hybrid load at generation %d (mean $\pm$ SD of %d replicates)' %(maxgen,nreps))
					# plt.xlabel('number of ancestral mutations')
					# plt.ion()

					#loop over all n_muts values
					i = 0
					while i < len(n_mut_list[l]):

						n_muts = n_mut_list[l][i] #set number of mutations in ancestor (ie how much SGV)

						# hyloads = [0] * nreps #initialize vector to store hybrid loads in from each replicate

						#loop over all replicates
						rep = 0
						while rep < nreps:

							#make ancestor
							if n_muts > 0:
								pop = np.random.binomial(1, p_mut, (N, n_muts)) #p_mut chance of having each of n_muts mutations, for all K individuals
								mut = np.random.normal(0, alpha, (n_muts, n)) #create n_muts mutations, each with a random normal phenotypic effect in each n dimension with mean 0 and sd alpha
							else: #de novo only, even if p_mut>0
								pop = np.array([[1]] * N)
								mut = np.array([[0] * n])

							#found identical populations
							[popfound, mutfound] = found(N, N_adapt, pop, mut, remove_lost, remove)
							[pop1, mut1] = [popfound, mutfound]
							[pop2, mut2] = [popfound, mutfound]

							#intitialize generation counter
							gen = 0

							#run until maxgen
							while gen < maxgen + 1:

								# genotype to phenotype
								phenos1 = np.dot(pop1, mut1) #sum mutations held by each individual
								phenos2 = np.dot(pop2, mut2) #sum mutations held by each individual

								# phenotype to fitness
								w1 = fitness(phenos1, theta1, sigma_adapt)
								w2 = fitness(phenos2, theta2, sigma_adapt)

								# wright-fisher (multinomial) sampling
								parents1 = np.random.multinomial(N_adapt, w1/sum(w1)) #number of times each parent chosen
								off1 = np.repeat(pop1, parents1, axis=0) #offspring genotypes
								parents2 = np.random.multinomial(N_adapt, w2/sum(w2)) #number of times each parent chosen
								off2 = np.repeat(pop2, parents2, axis=0) #offspring genotypes

								# mating and recombination
								off1 = recomb(off1)
								off2 = recomb(off2)

								# mutation and population update
								[pop1, mut1] = mutate(off1, u_adapt, alpha_adapt, n, mut1)
								[pop2, mut2] = mutate(off2, u_adapt, alpha_adapt, n, mut2)

								# remove lost mutations (all zero columns in pop)
								[pop1, mut1] = remove_muts(remove, remove_lost, pop1, mut1, mutfound)
								[pop2, mut2] = remove_muts(remove, remove_lost, pop2, mut2, mutfound)

								# go to next generation
								gen += 1

							#make variables to hold offspring phenotypes
							offphenos = dict()
							offpheno = []

							#make each of nHybrids hybrids
							for k in range(nHybrids):
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
							mean_offpheno = np.mean(offpheno, axis=0)
							offfit  = fitness(offpheno, mean_offpheno, sigma_adapt)
							hyload = -np.mean(np.log(offfit)) #hybrid load as defined by Chevin et al 2014

							#calculate genetic parallelism across ancestrally-segregating loci that have been segregating in adapting populations since divergence
							p = sum(pop1[:, len(mutfound)-n_muts:n_muts]) / N_adapt #frequency of derived alleles in pop1
							q = sum(pop2[:, len(mutfound)-n_muts:n_muts]) / N_adapt #frequency of derived alleles in pop2
							EH = np.mean(p*(1-q)+(1-p)*q) #expected heterozygosity in hybrids
							P = 1 - EH #our measure of parallelism (expected homozygosity in F1 diploid hybrids; will be 0 if fixed all different alleles, 1 if fixed all same alleles)

							#calculate genetic parallelism across ancestrally-shared segregating that have been segregating in adapting populations since divergence plus those loci that have mutations unique to one adapting population
							p = sum(pop1[:, len(mutfound)-n_muts:n_muts]) / N_adapt #frequency of derived alleles in pop1 
							q = sum(pop2[:, len(mutfound)-n_muts:n_muts]) / N_adapt #frequency of derived alleles in pop2
							EH_1 = p*(1-q)+(1-p)*q #expected heterozygosities at those loci
							p = sum(pop1[:, len(mutfound):]) / N_adapt #frequency of unique derived alleles in pop1 = expected heterozygosity at loci with mutations unique to pop1
							q = sum(pop2[:, len(mutfound):]) / N_adapt #frequency of unique derived alleles in pop2 = expected heterozygosity at loci with mutations unique to pop2
							EH_2 = np.append(p,q) #list of expected heterozygosities at unique loci
							EH = np.mean(np.append(EH_1,EH_2)) #expected heterozygosity across all loci considered
							P_all = 1 - EH #our measure of parallelism

							#print an update
							print('N=%d, sigma=%.2f, n=%d, angle=%r, rep=%d, n_muts=%d, distance between optima=%.3f,  hybrid load=%.3f, parallelism=%.4f, parallelism_all=%.4f' %(N_adapt, sigma_adapt, n, round(angles[j]*180/math.pi,2), rep+1, n_muts, opt_dist * (2*(1-math.cos(angles[j])))**(0.5), hyload, P, P_all)) 
							
							#save data
							write_data_to_output(fileHandles, [round(angles[j]*180/math.pi,2), rep+1, n_muts, opt_dist * (2*(1-math.cos(angles[j])))**(0.5), hyload, P, P_all])

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
