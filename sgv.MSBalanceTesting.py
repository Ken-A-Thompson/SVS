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
	sim_id = 'n%d_K%d_alpha%.2f_B%d_u%.3f' %(n, K, alpha, B, u)
	outfile_A = open("%s/MSBalance_%s.csv" %(data_dir, sim_id), "w")
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

def remove_muts(remove_lost, pop, mut):
	"""
	This function creates mutations and updates population
	"""
	if remove_lost:
		keep = pop.any(axis=0)
		mut = mut[keep]
		pop = pop[:, keep]
	return [pop, mut]

######################################################################
##PARAMETERS##
######################################################################

nreps = 1 #number of replicates for each set of parameters
n = 2 #phenotypic dimensions (positive integer >=1)
data_dir = 'data'

K = 1000 #number of individuals (positive integer >=1)
Es = 0.05 #expected selective effect
alpha = 2*Es/n #mutational sd (positive real number)

theta = [0] * n #optimum

B = 2 #number of offspring per generation per parent (positive integer)
u = 0.01 #mutation probability per generation per genome (0<u<1)

maxgen = 10000 #total number of generations population adapts for
gen_rec = 1000 #print and save after this many generations

remove_lost = True #If true, remove mutations that are lost (0 for all individuals)

######################################################################
##FUNCTION FOR POPULATION TO ADAPT##
######################################################################

def main():

	#plot hybrid load
	# plt.axis([-0.1, 0, 0, 10])
	# plt.ion()

	# open output files
	fileHandles = open_output_files(n, K, alpha, B, u, data_dir) 

	#loop over all replicates
	rep = 0
	while rep < nreps:

		pop = np.array([[1]] * K) #start all individuals with one "mutation"
		mut = np.array([[0] * n]) #but let the "mutation" do nothing (ie stay at origin)

		#intitialize generation counter
		gen = 0

		#run until maxgen
		while gen < maxgen + 1:

			# genotype to phenotype
			phenos = np.dot(pop, mut) #sum mutations held by each individual

			# viability selection
			surv = viability(phenos, theta, pop, K)

			# end simulation if population extinct (or unable to produce offspring)        
			if len(surv) < 2: 
				print("Extinct")              
				break 
				
			# meiosis
			off = recomb(surv, B)

			# mutation and population update
			[pop, mut] = mutate(off, u, alpha, n, mut)

			# remove lost mutations (all zero columns in pop)
			[pop, mut] = remove_muts(remove_lost, pop, mut)

			if gen % gen_rec == 0:
				
				w = survival(np.linalg.norm(mut - theta, axis=1)) #fitness of each mutation alone
				# w = survival(np.linalg.norm(np.dot(pop,mut) - theta, axis=1)) #fitness of each individual
				s = np.log(w)

				fig, (ax1, ax2, ax3) = plt.subplots(3,1)

				ax1.hist(s) #histogram of s of individual mutations
				ax2.hist(sum(pop)[1:]/len(pop)) #histogram of mutant frequencies (site frequency spectrum)
				ax3.hist(np.sum(pop[:,1:],axis=1)) #histogram of number of mutations in a single individual
				plt.pause(0.01)

				#print update
				print('gen=%d n=%d' %(gen, len(pop))) 
				
				#save data
				for x in range(len(w)):
					write_data_to_output(fileHandles, [gen] + [s[x]]) #gen, s

			# go to next generation
			gen += 1

		# go to next rep
		rep += 1

	# cleanup
	close_output_files(fileHandles)

######################################################################
##RUNNING ADAPTATION FUNCTION##
######################################################################    
	
start = time.time()
main()
end = time.time()
print('this took %.2f seconds to complete' %(end-start))